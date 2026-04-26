"""Kalman filter: 4-state with bounded transient buffer, free TDEE drift.

States:
  FM (lbs)     — fat mass
  LM (lbs)     — lean mass
  BUF (lbs)    — non-glycogen, non-sodium transient buffer (gut contents, hydration)
  TDEE (cal/d) — total daily expenditure, free random walk

Design notes:
  - LM is anchored to composition scans + the AA training term; otherwise evolves
    via the lean fraction of energy balance with workout-modulated Forbes (CB).
  - BUF is a generalized water/gut buffer that decays toward zero with HL ~5 days
    and is capped at ±3 lb. It absorbs daily weight noise that the P1 glycogen+sodium
    correction does not capture.
  - FM has hard rate caps from physiology: gain at +0.5 lb/day (Levine NEAT),
    loss at FM × 31/3500 = FM × 0.0089 lb/day (Alpert mobilization ceiling).
  - TDEE is a free random walk with daily-step rate cap (no anchor, no mean reversion).
  - Calorimetry RMR observations are NOT fed to the filter — they are held out
    for posterior validation.

Sweep: σ_TDEE ∈ {15, 30, 50} × max_dTDEE ∈ {30, 60, 100}.
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

# ── Physical constants ────────────────────────────────────────────────
# Hall's body-composition-model energy densities (doi:10.1152/ajpendo.00559.2006).
# These plus Forbes-coupled partition give correct energy conservation:
#   intake - TDEE = dFM × E_FAT + dLM × E_LEAN,
#   dFM / dW = ff(FM)   [Forbes]
# A naive partition (dFM = deficit × ff / 3500, dLM = deficit × (1-ff) / E_LEAN)
# *over-estimates* lean loss in deficit because the deficit energy is split
# independently between two tissues at different densities, violating Forbes.
E_FAT = 3500.0                 # cal/lb fat (textbook)
E_LEAN = 1100.0                # cal/lb lean (Hall: ~1100 cal/lb, mostly water+protein)
FORBES_C_KG = 10.4
LBS_TO_KG = 0.453592

ALPERT_RATE = 31.0 / E_FAT             # 0.00886 lb/lb-FM/day max FM mobilization
DFM_GAIN_CAP = 0.5             # lb/day, peak overfeeding
DLM_LOSS_CAP = 0.5             # lb/day, extreme protein turnover (rarely hit)
DLM_GAIN_CAP = 0.15            # lb/day, AA per-session ~0.09 + small headroom

BUF_HL = 5.0                   # days
BUF_CAP = 3.0                  # lbs
BUF_DECAY = np.exp(-np.log(2) / BUF_HL)

# ── Process noise (Q) ─────────────────────────────────────────────────
Q_FM = 0.005                   # lb²/day — daily FM std ≈ 0.07 lb
Q_LM = 0.005                   # lb²/day — daily LM std ≈ 0.07 lb
Q_BUF = 0.16                   # lb²/day — daily BUF std ≈ 0.4 lb

# ── Observation noise (R) ─────────────────────────────────────────────
R_WEIGHT = 0.97                # lb² — measured from consecutive-day weight variance
                               #       after P1 glycogen+sodium correction
R_FM_SCAN = 2.25               # lb² — BOD POD/InBody at ~1.5 lb std
R_LM_SCAN = 2.25               # lb² — same
R_BOOT_FM = 0.50               # tighter scan noise during burn-in for stable init
R_BOOT_LM = 0.50

# ── Training term (finding AA) ────────────────────────────────────────
TRAINING_DELTA = 0.091         # lb LM per session
TRAINING_HL = 275.0            # days
TRAINING_DECAY = np.log(2) / TRAINING_HL

# ── Workout-modulated Forbes partition (finding CB) ──────────────────
# CB regression on cuts: each workout/wk reduces lean fraction of weight loss
# by ~8.4 percentage points beyond Forbes baseline. Implementation: shift ff
# (fat fraction of weight change) by alpha × workouts_per_week_30d.
# Set FORBES_WORKOUT_ALPHA = 0 to recover the base Forbes partition.
FORBES_WORKOUT_ALPHA = 0.084   # per workout/wk, additive to ff

# ── Sweep grid ────────────────────────────────────────────────────────
SIGMA_TDEE_GRID = [15, 30, 50]    # cal/day daily noise
MAX_DTDEE_GRID = [30, 60, 100]    # cal/day hard cap


def load_data():
    smooth = pd.read_csv(ROOT / "analysis" / "P1_smoothed_weight.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    strength = pd.read_csv(ROOT / "workout" / "strength.csv", parse_dates=["date"])
    comp_scans = pd.read_csv(ROOT / "composition" / "composition.csv", parse_dates=["date"])
    rmr = pd.read_csv(ROOT / "RMR" / "rmr.csv", parse_dates=["date"])
    p2 = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    temp = pd.read_csv(ROOT / "analysis" / "AS_temperature_baseline.csv", parse_dates=["date"])
    return smooth, intake, strength, comp_scans, rmr, p2, temp


def forbes_fat_fraction(fm_lbs, workouts_per_wk=0.0):
    """Forbes fat-fraction-of-weight-change, modulated by recent workout density.

    Base Forbes: ff_base = 1 - C²/(C + FM_kg)² (free-living calibration).
    Workout modifier: ff = clip(ff_base + alpha × workouts/wk, 0, 1).
    Per CB, alpha = +0.084 — a 1×/wk training cadence shifts ~8 pp from lean
    catabolism to fat catabolism in deficit (and the symmetric direction in surplus).
    """
    fm_kg = max(float(fm_lbs), 5.0) * LBS_TO_KG
    ff_base = 1.0 - (FORBES_C_KG ** 2) / (FORBES_C_KG + fm_kg) ** 2
    ff_modified = ff_base + FORBES_WORKOUT_ALPHA * workouts_per_wk
    return float(np.clip(ff_modified, 0.0, 1.0))


def compute_training_increment(dates, strength_dates):
    """Daily LM increment from the AA training model.

    On a workout day, add TRAINING_DELTA. The exponential decay of past sessions
    is folded into LM dynamics naturally — we don't add a global stock; the
    accumulated effect emerges from the daily increments minus the decay-implied
    bleed. To avoid double-counting decay (LM doesn't decay on its own here), we
    add the *net* daily increment: delta on workout days minus a small bleed
    on non-workout days proportional to time-since-last-session.
    """
    n = len(dates)
    workout_set = set(np.datetime64(d, "D") for d in strength_dates)
    out = np.zeros(n)
    accumulated = 0.0
    for i in range(n):
        d = np.datetime64(dates[i], "D")
        if d in workout_set:
            out[i] += TRAINING_DELTA
            accumulated += TRAINING_DELTA
        # bleed: accumulated × decay rate
        bleed = accumulated * TRAINING_DECAY
        out[i] -= bleed
        accumulated -= bleed
    return out


def build_inputs(smooth, intake, strength, comp_scans, p2):
    dates = intake["date"].values
    calories = intake["calories"].values
    n = len(dates)

    expected_rmr = p2.set_index("date")["expected_rmr"].reindex(intake["date"]).values

    # Sparse weight observations from P1 (glycogen+sodium-corrected)
    obs_map = smooth.set_index("date")["smoothed_weight_lbs"]
    weight_obs = np.full(n, np.nan)
    for i, d in enumerate(dates):
        if d in obs_map.index:
            v = obs_map[d]
            if not np.isnan(v):
                weight_obs[i] = v

    # Sparse composition scans (FM and LM at ~70 dates)
    fm_scan = np.full(n, np.nan)
    lm_scan = np.full(n, np.nan)
    scan_idx = {pd.Timestamp(d): i for i, d in enumerate(dates)}
    for _, row in comp_scans.iterrows():
        d = pd.Timestamp(row["date"])
        if d in scan_idx:
            i = scan_idx[d]
            fm_scan[i] = row["fat_mass_lbs"]
            # CLAUDE.md: derive LM as weight - fat_mass (raw "Lean Mass" col is dry lean)
            lm_scan[i] = row["weight_lbs"] - row["fat_mass_lbs"]

    # Daily training increment for LM dynamics
    training_inc = compute_training_increment(dates, strength["date"].values)

    # Workouts per week, 30-day trailing — used to modulate Forbes partition
    workout_set = set(pd.to_datetime(strength["date"]).dt.normalize())
    workout_indicator = np.array([1 if pd.Timestamp(d).normalize() in workout_set else 0
                                  for d in dates], dtype=float)
    workouts_per_wk = (
        pd.Series(workout_indicator).rolling(30, min_periods=1).sum().values * 7.0 / 30.0
    )

    return dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr, training_inc, workouts_per_wk


def kalman_v4_forward(calories, weight_obs, fm_scan, lm_scan, expected_rmr,
                      training_inc, workouts_per_wk, fm_start, lm_start,
                      sigma_tdee, max_dtdee):
    """Forward EKF on 4-state (FM, LM, BUF, TDEE)."""
    n = len(calories)
    q_tdee = sigma_tdee ** 2

    # Initial state — anchored to first scan
    first_scan = np.where(~np.isnan(fm_scan))[0][0]
    init_tdee = (expected_rmr[first_scan] * 1.15
                 if not np.isnan(expected_rmr[first_scan]) else 2100.0)

    x = np.array([fm_start, lm_start, 0.0, init_tdee])
    P = np.diag([4.0, 4.0, 1.0, 200.0 ** 2])  # wide TDEE prior

    x_filt = np.full((n, 4), np.nan)
    P_filt = np.full((n, 4, 4), np.nan)
    x_pred_arr = np.full((n, 4), np.nan)
    P_pred_arr = np.full((n, 4, 4), np.nan)
    innovations = np.full(n, np.nan)  # weight innovations

    x_filt[first_scan] = x
    P_filt[first_scan] = P
    x_pred_arr[first_scan] = x
    P_pred_arr[first_scan] = P

    Q = np.diag([Q_FM, Q_LM, Q_BUF, q_tdee])

    for i in range(first_scan, n - 1):
        fm, lm, buf, tdee = x
        ff = forbes_fat_fraction(fm, workouts_per_wk[i])
        deficit = (calories[i] if not np.isnan(calories[i]) else 0.0) - tdee

        # Proper Forbes-coupled energy partition (energy-conserving):
        #   deficit = dFM × E_FAT + dLM × E_LEAN
        #   dFM / dW = ff
        #   ⇒ dW = deficit / (ff × E_FAT + (1 - ff) × E_LEAN)
        e_per_lb = ff * E_FAT + (1 - ff) * E_LEAN
        dW_raw = deficit / e_per_lb
        dFM_raw = ff * dW_raw
        dLM_raw = (1 - ff) * dW_raw + training_inc[i]

        # Apply hard rate caps
        dFM = np.clip(dFM_raw, -fm * ALPERT_RATE, DFM_GAIN_CAP)
        dLM = np.clip(dLM_raw, -DLM_LOSS_CAP, DLM_GAIN_CAP)

        x_pred = np.array([
            max(fm + dFM, 5.0),
            max(lm + dLM, 50.0),
            BUF_DECAY * buf,
            tdee,  # random walk; η added via Q
        ])

        # Jacobian (proper coupled-Forbes partition; ff held locally constant)
        F = np.array([
            [1.0, 0.0, 0.0, -ff / e_per_lb],
            [0.0, 1.0, 0.0, -(1 - ff) / e_per_lb],
            [0.0, 0.0, BUF_DECAY, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
        P_pred = F @ P @ F.T + Q

        x_pred_arr[i + 1] = x_pred
        P_pred_arr[i + 1] = P_pred

        x = x_pred.copy()
        P = P_pred.copy()

        # Sequential observation updates
        # 1. Weight: H = [1, 1, 1, 0]
        if not np.isnan(weight_obs[i + 1]):
            H = np.array([[1.0, 1.0, 1.0, 0.0]])
            y = weight_obs[i + 1]
            y_pred = (H @ x)[0]
            innovation = y - y_pred
            S = (H @ P @ H.T + R_WEIGHT)[0, 0]
            K = (P @ H.T) / S
            x = x + (K * innovation).flatten()
            P = P - np.outer(K.flatten(), H @ P)
            innovations[i + 1] = innovation

        # 2. FM scan: H = [1, 0, 0, 0]
        if not np.isnan(fm_scan[i + 1]):
            H = np.array([[1.0, 0.0, 0.0, 0.0]])
            y = fm_scan[i + 1]
            r_use = R_BOOT_FM if (i + 1 - first_scan) < 30 else R_FM_SCAN
            S = (H @ P @ H.T + r_use)[0, 0]
            K = (P @ H.T) / S
            x = x + (K * (y - (H @ x)[0])).flatten()
            P = P - np.outer(K.flatten(), H @ P)

        # 3. LM scan: H = [0, 1, 0, 0]
        if not np.isnan(lm_scan[i + 1]):
            H = np.array([[0.0, 1.0, 0.0, 0.0]])
            y = lm_scan[i + 1]
            r_use = R_BOOT_LM if (i + 1 - first_scan) < 30 else R_LM_SCAN
            S = (H @ P @ H.T + r_use)[0, 0]
            K = (P @ H.T) / S
            x = x + (K * (y - (H @ x)[0])).flatten()
            P = P - np.outer(K.flatten(), H @ P)

        # Apply hard bounds post-update
        x[0] = max(x[0], 5.0)
        x[1] = max(x[1], 50.0)
        x[2] = np.clip(x[2], -BUF_CAP, BUF_CAP)
        # Rate cap on TDEE between consecutive filtered values
        prev_tdee = x_filt[i, 3] if not np.isnan(x_filt[i, 3]) else tdee
        x[3] = prev_tdee + np.clip(x[3] - prev_tdee, -max_dtdee, max_dtdee)

        x_filt[i + 1] = x
        P_filt[i + 1] = P

    return x_filt, P_filt, x_pred_arr, P_pred_arr, innovations, first_scan


def rts_smoother_v4(x_filt, P_filt, x_pred, P_pred, first_scan, workouts_per_wk):
    n = len(x_filt)
    x_smooth = np.copy(x_filt)
    P_smooth = np.copy(P_filt)

    valid = ~np.isnan(x_filt[:, 0])
    last_valid = np.where(valid)[0][-1]

    for i in range(last_valid - 1, first_scan - 1, -1):
        if np.isnan(x_filt[i, 0]) or np.isnan(P_pred[i + 1, 0, 0]):
            continue
        fm = x_filt[i, 0]
        ff = forbes_fat_fraction(fm, workouts_per_wk[i])
        e_per_lb = ff * E_FAT + (1 - ff) * E_LEAN
        F = np.array([
            [1.0, 0.0, 0.0, -ff / e_per_lb],
            [0.0, 1.0, 0.0, -(1 - ff) / e_per_lb],
            [0.0, 0.0, BUF_DECAY, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
        try:
            P_pred_inv = np.linalg.inv(P_pred[i + 1])
        except np.linalg.LinAlgError:
            continue
        G = P_filt[i] @ F.T @ P_pred_inv
        x_smooth[i] = x_filt[i] + G @ (x_smooth[i + 1] - x_pred[i + 1])
        P_smooth[i] = P_filt[i] + G @ (P_smooth[i + 1] - P_pred[i + 1]) @ G.T

    return x_smooth, P_smooth


def diagnostics(dates, x_smooth, innovations, fm_scan, lm_scan, calories, rmr_df,
                weight_obs, temp_df=None):
    """Return a dict of pre-registered fit metrics for the filter."""
    out = {}

    # 1. Innovation autocorrelation
    valid_innov = innovations[~np.isnan(innovations)]
    centered = valid_innov - valid_innov.mean()
    out["n_innov"] = len(valid_innov)
    out["innov_std"] = float(valid_innov.std())
    out["innov_acf1"] = float(np.corrcoef(centered[:-1], centered[1:])[0, 1])

    # 2. FM scan fit
    fm_errs = []
    for i in range(len(fm_scan)):
        if not np.isnan(fm_scan[i]) and not np.isnan(x_smooth[i, 0]):
            fm_errs.append(x_smooth[i, 0] - fm_scan[i])
    out["fm_scan_mae"] = float(np.mean(np.abs(fm_errs)))
    out["fm_scan_n"] = len(fm_errs)

    # 3. LM scan fit
    lm_errs = []
    for i in range(len(lm_scan)):
        if not np.isnan(lm_scan[i]) and not np.isnan(x_smooth[i, 1]):
            lm_errs.append(x_smooth[i, 1] - lm_scan[i])
    out["lm_scan_mae"] = float(np.mean(np.abs(lm_errs)))
    out["lm_scan_n"] = len(lm_errs)

    # 4. Calorimetry posterior (held out)
    ratios = []
    rmr_dates = pd.to_datetime(rmr_df["date"]).values
    rmr_kcal = rmr_df["rmr_kcal"].values
    for d, r in zip(rmr_dates, rmr_kcal):
        idx = np.searchsorted(dates, np.datetime64(d))
        if idx < len(x_smooth) and not np.isnan(x_smooth[idx, 3]):
            ratios.append(x_smooth[idx, 3] / r)
    out["tdee_rmr_n"] = len(ratios)
    out["tdee_rmr_mean"] = float(np.mean(ratios))
    out["tdee_rmr_std"] = float(np.std(ratios))
    out["tdee_rmr_min"] = float(np.min(ratios))
    out["tdee_rmr_max"] = float(np.max(ratios))
    out["tdee_rmr_in_band"] = float(np.mean([(0.95 <= r <= 1.30) for r in ratios]))

    # 5. Cumulative energy balance closure
    valid_cal = ~np.isnan(calories) & ~np.isnan(x_smooth[:, 3])
    cum_deficit = np.cumsum((calories - x_smooth[:, 3])[valid_cal])
    fm_first = x_smooth[np.where(valid_cal)[0][0], 0]
    fm_last = x_smooth[np.where(valid_cal)[0][-1], 0]
    fm_change_observed = fm_last - fm_first
    fm_change_implied = cum_deficit[-1] / E_FAT
    out["eb_residual_lb"] = float(fm_change_implied - fm_change_observed)

    # 6. BUF behavior
    buf = x_smooth[:, 2]
    valid_buf = buf[~np.isnan(buf)]
    out["buf_mean"] = float(np.mean(valid_buf))
    out["buf_std"] = float(np.std(valid_buf))
    out["buf_max_abs"] = float(np.max(np.abs(valid_buf)))
    out["buf_frac_under_1lb"] = float(np.mean(np.abs(valid_buf) < 1.0))
    out["buf_frac_under_cap"] = float(np.mean(np.abs(valid_buf) < BUF_CAP * 0.99))

    # 7. TDEE smoothness
    tdee = x_smooth[~np.isnan(x_smooth[:, 3]), 3]
    out["tdee_dt_std"] = float(np.std(np.diff(tdee)))
    out["tdee_min"] = float(np.min(tdee))
    out["tdee_max"] = float(np.max(tdee))
    out["tdee_mean"] = float(np.mean(tdee))

    # 8. TDEE vs corrected baseline temperature (held-out posterior check)
    if temp_df is not None:
        date_idx = {pd.Timestamp(d): i for i, d in enumerate(dates)}
        tdee_pairs = []
        temp_pairs = []
        for _, row in temp_df.iterrows():
            d = pd.Timestamp(row["date"])
            if d in date_idx:
                i = date_idx[d]
                if not np.isnan(x_smooth[i, 3]) and not np.isnan(row["baseline_mean"]):
                    tdee_pairs.append(x_smooth[i, 3])
                    temp_pairs.append(row["baseline_mean"])
        if len(tdee_pairs) >= 30:
            r = np.corrcoef(tdee_pairs, temp_pairs)[0, 1]
            # Partial correlation controlling for FM (since both might trend)
            tdee_arr = np.array(tdee_pairs)
            temp_arr = np.array(temp_pairs)
            fm_pairs = []
            for _, row in temp_df.iterrows():
                d = pd.Timestamp(row["date"])
                if d in date_idx:
                    i = date_idx[d]
                    if not np.isnan(x_smooth[i, 3]) and not np.isnan(row["baseline_mean"]):
                        fm_pairs.append(x_smooth[i, 0])
            fm_arr = np.array(fm_pairs)
            # Residualize both vs FM
            def resid(y, x):
                X = np.column_stack([x, np.ones_like(x)])
                beta, *_ = np.linalg.lstsq(X, y, rcond=None)
                return y - X @ beta
            tdee_r = resid(tdee_arr, fm_arr)
            temp_r = resid(temp_arr, fm_arr)
            r_partial = np.corrcoef(tdee_r, temp_r)[0, 1]
            out["tdee_temp_n"] = len(tdee_pairs)
            out["tdee_temp_r"] = float(r)
            out["tdee_temp_r_partial_fm"] = float(r_partial)
        else:
            out["tdee_temp_n"] = len(tdee_pairs)
            out["tdee_temp_r"] = float("nan")
            out["tdee_temp_r_partial_fm"] = float("nan")

    return out


def run_one(dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr, training_inc,
            workouts_per_wk, fm_start, lm_start, rmr_df, sigma_tdee, max_dtdee, temp_df=None):
    x_filt, P_filt, x_pred, P_pred, innov, first_scan = kalman_v4_forward(
        calories, weight_obs, fm_scan, lm_scan, expected_rmr,
        training_inc, workouts_per_wk, fm_start, lm_start, sigma_tdee, max_dtdee,
    )
    x_smooth, P_smooth = rts_smoother_v4(x_filt, P_filt, x_pred, P_pred, first_scan, workouts_per_wk)
    diag = diagnostics(dates, x_smooth, innov, fm_scan, lm_scan, calories, rmr_df,
                       weight_obs, temp_df)
    return x_filt, x_smooth, innov, diag


def main():
    smooth, intake, strength, comp_scans, rmr_df, p2, temp_df = load_data()
    dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr, training_inc, workouts_per_wk = \
        build_inputs(smooth, intake, strength, comp_scans, p2)

    first_scan = np.where(~np.isnan(fm_scan))[0][0]
    fm_start = fm_scan[first_scan]
    lm_start = lm_scan[first_scan]
    print(f"Start: {str(dates[first_scan])[:10]}  FM={fm_start:.1f}  LM={lm_start:.1f}")
    print(f"Sweep: σ_TDEE × max_dTDEE")
    print(f"  σ_TDEE ∈ {SIGMA_TDEE_GRID} cal/day")
    print(f"  max_dTDEE ∈ {MAX_DTDEE_GRID} cal/day")
    print()

    results = []
    best_diag = None
    best_x = None
    best_innov = None
    best_key = None
    best_score = float("inf")

    print(f"{'σ_TDEE':>7} {'max_d':>6} {'FM_MAE':>7} {'LM_MAE':>7} {'innov_acf':>10} "
          f"{'innov_std':>10} {'TDEE/RMR':>9} {'(std)':>6} {'inband':>7} {'EB(lb)':>7} "
          f"{'BUF_max':>8} {'r(T,°F)':>8} {'r|FM':>6}")
    for sigma_tdee in SIGMA_TDEE_GRID:
        for max_dtdee in MAX_DTDEE_GRID:
            x_filt, x_smooth, innov, diag = run_one(
                dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr,
                training_inc, workouts_per_wk, fm_start, lm_start, rmr_df,
                sigma_tdee, max_dtdee, temp_df=temp_df,
            )
            results.append({
                "sigma_tdee": sigma_tdee,
                "max_dtdee": max_dtdee,
                **diag,
            })
            print(
                f"{sigma_tdee:>7} {max_dtdee:>6} "
                f"{diag['fm_scan_mae']:>7.2f} {diag['lm_scan_mae']:>7.2f} "
                f"{diag['innov_acf1']:>10.3f} {diag['innov_std']:>10.3f} "
                f"{diag['tdee_rmr_mean']:>9.3f} {diag['tdee_rmr_std']:>6.3f} "
                f"{diag['tdee_rmr_in_band']:>7.2f} {diag['eb_residual_lb']:>7.1f} "
                f"{diag['buf_max_abs']:>8.2f} "
                f"{diag.get('tdee_temp_r', float('nan')):>8.3f} "
                f"{diag.get('tdee_temp_r_partial_fm', float('nan')):>6.3f}"
            )
            # Composite score: rank by held-out metrics with sane priors
            score = (diag['fm_scan_mae'] + diag['lm_scan_mae']
                     + max(diag['innov_acf1'], 0) * 5
                     + max(0, abs(diag['tdee_rmr_mean'] - 1.15) - 0.10) * 10
                     + max(0, 0.05 - diag['tdee_rmr_std']) * 20
                     - max(0, diag.get('tdee_temp_r_partial_fm', 0)) * 5)
            if score < best_score:
                best_score = score
                best_key = (sigma_tdee, max_dtdee)
                best_x = x_smooth
                best_diag = diag
                best_innov = innov

    # Save best run.
    print(f"\nBest: σ_TDEE={best_key[0]}, max_dTDEE={best_key[1]} (composite score {best_score:.2f})")
    out = pd.DataFrame({
        "date": dates,
        # smoothed states (filtered/std placeholders kept null for downstream interface)
        "fat_mass_lbs_filtered": np.round(best_x[:, 0], 2),
        "fat_mass_std_filtered": np.nan,
        "tdee_filtered": np.round(best_x[:, 3], 0),
        "tdee_std_filtered": np.nan,
        "fat_mass_lbs": np.round(best_x[:, 0], 2),
        "fat_mass_std": np.nan,
        "tdee": np.round(best_x[:, 3], 0),
        "tdee_std": np.nan,
        "innovation": np.round(best_innov, 3),
        # additional state columns
        "lm_lbs": np.round(best_x[:, 1], 2),
        "buf_lbs": np.round(best_x[:, 2], 3),
    })
    out["date"] = pd.to_datetime(out["date"]).dt.strftime("%Y-%m-%d")
    path = ROOT / "analysis" / "P4_kalman_daily.csv"
    out.to_csv(path, index=False)
    print(f"Wrote {len(out)} rows to {path}")

    # Save sweep summary
    pd.DataFrame(results).to_csv(ROOT / "analysis" / "P4_kalman_sweep.csv", index=False)
    print(f"Wrote sweep results to analysis/P4_kalman_sweep.csv")

    # Calorimetry per-point posterior dump for the best config
    print(f"\n=== Calorimetry posterior (best config — held out, never seen by filter) ===")
    rmr_dates = pd.to_datetime(rmr_df["date"]).values
    rmr_kcal = rmr_df["rmr_kcal"].values
    print(f"  {'date':<12} {'RMR':>5} {'TDEE':>5} {'ratio':>6}")
    ratios = []
    for d, r in zip(rmr_dates, rmr_kcal):
        idx = np.searchsorted(dates, np.datetime64(d))
        if idx < len(best_x) and not np.isnan(best_x[idx, 3]):
            t = best_x[idx, 3]
            ratios.append(t / r)
            print(f"  {str(d)[:10]:<12} {r:>5.0f} {t:>5.0f} {t/r:>6.2f}")
    print(f"  Mean {np.mean(ratios):.3f}  std {np.std(ratios):.3f}  "
          f"range [{np.min(ratios):.2f}, {np.max(ratios):.2f}]")

    # Per-scan FM and LM error (best config)
    print(f"\n=== Composition scan errors (best config) ===")
    print(f"  {'date':<12} {'FM_obs':>6} {'FM_filt':>7} {'FM_err':>7} "
          f"{'LM_obs':>6} {'LM_filt':>7} {'LM_err':>7}")
    for i in range(len(fm_scan)):
        if not np.isnan(fm_scan[i]):
            d = pd.Timestamp(dates[i]).strftime("%Y-%m-%d")
            print(f"  {d:<12} {fm_scan[i]:>6.1f} {best_x[i, 0]:>7.1f} "
                  f"{best_x[i, 0] - fm_scan[i]:>+7.1f} "
                  f"{lm_scan[i]:>6.1f} {best_x[i, 1]:>7.1f} "
                  f"{best_x[i, 1] - lm_scan[i]:>+7.1f}")


if __name__ == "__main__":
    main()
