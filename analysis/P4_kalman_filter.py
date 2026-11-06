"""Kalman filter: 4-state with bounded transient buffer, free TDEE drift.

States:
  FM (lbs)     — fat mass
  LM (lbs)     — lean mass
  BUF (lbs)    — non-glycogen, non-sodium transient buffer (gut contents, hydration)
  TDEE (cal/d) — total daily expenditure, free random walk

Design notes:
  - LM is anchored to composition scans; between scans it evolves via the lean
    fraction of energy balance with workout-modulated Forbes (CB).
  - BUF is a generalized water/gut buffer that decays toward zero with HL ~5 days
    and is capped at ±3 lb. It absorbs daily weight noise that the P1 glycogen+sodium
    correction does not capture.
  - FM has hard rate caps from physiology: gain at +0.5 lb/day (Levine NEAT),
    loss at FM × 31/3500 = FM × 0.0089 lb/day (Alpert mobilization ceiling).
  - TDEE is a free random walk with daily-step rate cap (no anchor, no mean reversion).
  - Calorimetry RMR observations are NOT fed to the filter — they are held out
    for posterior validation.

Tuning: σ_TDEE, Q_BUF, R_WEIGHT and the intake model come from analysis/results/
DM_intake_calibration.json (one likelihood over all assimilated observations);
the σ_TDEE × max_dTDEE sweep printed by main() is diagnostic only.
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

BUF_HL = 5.0                   # days (default; DM fits it, see FILTER)
BUF_CAP = 3.0                  # lbs
BUF_DECAY = np.exp(-np.log(2) / BUF_HL)
# Intake-driven transient: a day's intake above/below the record mean moves
# next morning's weight through gut contents and water that P1's glycogen and
# sodium corrections do not capture. Modelled as a deterministic input to BUF,
# BETA_BUF lb per 1000 consumed cal of deviation, decaying with BUF's half-life.
# Without this term the filter reads the transient as stored energy and the
# intake multiplier runs away against Q_BUF (DM, 2026-11). Default 0 = legacy.
BETA_BUF = 0.0

# ── Process noise (Q) ─────────────────────────────────────────────────
Q_FM = 0.005                   # lb²/day — daily FM std ≈ 0.07 lb
Q_LM = 0.005                   # lb²/day — daily LM std ≈ 0.07 lb
Q_BUF = 0.32                   # lb²/day — daily BUF std ≈ 0.57 lb
                               # (raised from 0.16 to absorb daily water/gut wobble that
                               # otherwise leaked into TDEE drift; biological prior says
                               # daily metabolism doesn't actually wiggle that much)

# ── Observation noise (R) ─────────────────────────────────────────────
R_WEIGHT = 0.97                # lb² — measured from consecutive-day weight variance
                               #       after P1 glycogen+sodium correction
R_FM_SCAN = 2.25               # lb² — BOD POD/InBody at ~1.5 lb std
R_LM_SCAN = 2.25               # lb² — same
R_BOOT_FM = 0.50               # tighter scan noise during burn-in for stable init
R_BOOT_LM = 0.50

# ── Intake correction (single source of truth: analysis/intake_model.py) ──
# Logged intake is converted to consumed energy by `true_intake`, whose form
# and parameters are selected by analysis/DM_intake_calibration.py (innovation
# likelihood over candidate correction forms) and stored in
# analysis/results/DM_intake_calibration.json. Nothing is typed here.
#
# Why not energy-balance closure: TDEE is a free random walk, so scaling intake
# by m raises TDEE by (m-1)*intake almost day for day (measured leakage
# dEB/dm = 24.9 lb per unit m against 3330.7 lb/unit with TDEE fixed; the
# filter absorbs 99.25% of any level change). The information about the scale
# lives in the weight innovation sequence instead: a wrong scale mis-predicts
# the weight RESPONSE to intake variation. See CR/DM. CQ estimates the same
# quantity by regression outside the filter; CQ and CR/DM are consistency
# checks on each other, not independent confirmations.
sys.path.insert(0, str(ROOT / "analysis"))
from intake_model import true_intake, INTAKE_MULTIPLIER, CALIBRATION, FILTER, describe as describe_intake_model, describe_filter  # noqa: E402
from results import save as save_result  # noqa: E402

# ── Workout-modulated Forbes partition (finding CB) ──────────────────
# CB regression on cuts: each workout/wk reduces lean fraction of weight loss
# by ~8.4 percentage points beyond Forbes baseline. Implementation: shift ff
# (fat fraction of weight change) by alpha × workouts_per_week_30d.
# Set FORBES_WORKOUT_ALPHA = 0 to recover the base Forbes partition.
FORBES_WORKOUT_ALPHA = 0.084   # per workout/wk, additive to ff

# ── Sweep grid ────────────────────────────────────────────────────────
SIGMA_TDEE_GRID = [15, 25, 35, 50]    # cal/day daily noise
MAX_DTDEE_GRID = [60, 100]            # cal/day hard cap


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


def build_inputs(smooth, intake, strength, comp_scans, p2, calories_true=None):
    """`calories_true` lets calibration scripts (CR, DM) inject their own
    corrected-intake series; default is the canonical intake_model."""
    dates = intake["date"].values
    calories = (np.asarray(calories_true, dtype=float) if calories_true is not None
                else true_intake(intake))
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

    # Workouts per week, 30-day trailing — used to modulate Forbes partition (CB)
    workout_set = set(pd.to_datetime(strength["date"]).dt.normalize())
    workout_indicator = np.array([1 if pd.Timestamp(d).normalize() in workout_set else 0
                                  for d in dates], dtype=float)
    workouts_per_wk = (
        pd.Series(workout_indicator).rolling(30, min_periods=1).sum().values * 7.0 / 30.0
    )

    return dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr, workouts_per_wk


def predict_step(x, P, cal_i, wpw_i, Q, buf_decay, beta_buf, cal_ref):
    """One-day state prediction (shared by the forward filter and by DM's
    k-step-ahead prediction-error criterion)."""
    fm, lm, buf, tdee = x
    ff = forbes_fat_fraction(fm, wpw_i)
    deficit = (cal_i if not np.isnan(cal_i) else 0.0) - tdee

    # Proper Forbes-coupled energy partition (energy-conserving):
    #   deficit = dFM × E_FAT + dLM × E_LEAN
    #   dFM / dW = ff
    #   ⇒ dW = deficit / (ff × E_FAT + (1 - ff) × E_LEAN)
    e_per_lb = ff * E_FAT + (1 - ff) * E_LEAN
    dW_raw = deficit / e_per_lb
    dFM_raw = ff * dW_raw
    dLM_raw = (1 - ff) * dW_raw

    # Apply hard rate caps
    dFM = np.clip(dFM_raw, -fm * ALPERT_RATE, DFM_GAIN_CAP)
    dLM = np.clip(dLM_raw, -DLM_LOSS_CAP, DLM_GAIN_CAP)

    cal_dev = ((cal_i - cal_ref) if not np.isnan(cal_i) else 0.0) / 1000.0
    x_pred = np.array([
        max(fm + dFM, 5.0),
        max(lm + dLM, 50.0),
        buf_decay * buf + beta_buf * cal_dev,
        tdee,  # random walk; η added via Q
    ])

    # Jacobian (proper coupled-Forbes partition; ff held locally constant)
    F = np.array([
        [1.0, 0.0, 0.0, -ff / e_per_lb],
        [0.0, 1.0, 0.0, -(1 - ff) / e_per_lb],
        [0.0, 0.0, buf_decay, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ])
    P_pred = F @ P @ F.T + Q
    return x_pred, P_pred


def kstep_prediction_errors(x_filt, P_filt, calories, weight_obs, workouts_per_wk, k,
                            sigma_tdee, q_buf=None, r_weight=None, beta_buf=None, buf_hl=None):
    """k-day-ahead weight prediction errors from each filtered state, with no
    intermediate observation updates (TDEE held at its filtered value, BUF
    decaying, FM/LM driven by the next k days of intake). Returns
    (innovations, variances) aligned to the start day. Used by DM: at k of a
    couple of weeks the gut/water transient has washed out and TDEE has not
    had time to drift, which is exactly where the intake scale is identified."""
    n = len(calories)
    q_buf = Q_BUF if q_buf is None else float(q_buf)
    r_weight = R_WEIGHT if r_weight is None else float(r_weight)
    beta_buf = BETA_BUF if beta_buf is None else float(beta_buf)
    buf_decay = BUF_DECAY if buf_hl is None else float(np.exp(-np.log(2) / float(buf_hl)))
    cal_ref = float(np.nanmean(calories))
    Q = np.diag([Q_FM, Q_LM, q_buf, sigma_tdee ** 2])
    H = np.array([1.0, 1.0, 1.0, 0.0])
    innov = np.full(n, np.nan); var = np.full(n, np.nan)
    for t in range(n - k):
        if np.isnan(x_filt[t, 0]) or np.isnan(weight_obs[t + k]) or np.isnan(weight_obs[t]):
            continue
        x, P = x_filt[t].copy(), P_filt[t].copy()
        for j in range(k):
            x, P = predict_step(x, P, calories[t + j], workouts_per_wk[t + j], Q, buf_decay, beta_buf, cal_ref)
        innov[t] = weight_obs[t + k] - H @ x
        var[t] = H @ P @ H + r_weight
    return innov, var


def kalman_v4_forward(calories, weight_obs, fm_scan, lm_scan, expected_rmr,
                      workouts_per_wk, fm_start, lm_start,
                      sigma_tdee, max_dtdee, collect=None, q_buf=None, r_weight=None,
                      beta_buf=None, buf_hl=None):
    """Forward EKF on 4-state (FM, LM, BUF, TDEE).

    `collect`: optional dict; if given, filled with the innovation and
    innovation-variance series for every assimilated channel ('innov_var' for
    weight, plus 'fm_innov'/'fm_var' and 'lm_innov'/'lm_var' for the scans) so
    callers can form the exact prediction-error-decomposition log-likelihood.
    `q_buf`, `r_weight`: override the module defaults (used by DM, which fits
    the noise hyperparameters jointly with the intake model).
    """
    n = len(calories)
    q_tdee = sigma_tdee ** 2
    q_buf = Q_BUF if q_buf is None else float(q_buf)
    r_weight = R_WEIGHT if r_weight is None else float(r_weight)
    beta_buf = BETA_BUF if beta_buf is None else float(beta_buf)
    buf_decay = BUF_DECAY if buf_hl is None else float(np.exp(-np.log(2) / float(buf_hl)))
    cal_ref = float(np.nanmean(calories))      # deviation reference for the BUF input

    # Initial state — anchored to first scan
    first_scan = np.where(~np.isnan(fm_scan))[0][0]
    # Initialise TDEE from the data rather than from an assumed activity factor.
    # The old form was expected_rmr * 1.15, which is the same activity-factor
    # prior that was removed from the tuning score. Mean corrected intake over
    # the first 60 logged days is a neutral starting point; the weight
    # observations move TDEE off it within days, and P is wide (200^2) anyway.
    _init_window = calories[first_scan:first_scan + 60]
    _init_window = _init_window[~np.isnan(_init_window)]
    init_tdee = float(np.mean(_init_window)) if len(_init_window) else 2100.0

    x = np.array([fm_start, lm_start, 0.0, init_tdee])
    P = np.diag([4.0, 4.0, 1.0, 200.0 ** 2])  # wide TDEE prior

    x_filt = np.full((n, 4), np.nan)
    P_filt = np.full((n, 4, 4), np.nan)
    x_pred_arr = np.full((n, 4), np.nan)
    P_pred_arr = np.full((n, 4, 4), np.nan)
    innovations = np.full(n, np.nan)  # weight innovations
    innov_var = np.full(n, np.nan)    # weight innovation variance S
    fm_innov, fm_var = np.full(n, np.nan), np.full(n, np.nan)
    lm_innov, lm_var = np.full(n, np.nan), np.full(n, np.nan)

    x_filt[first_scan] = x
    P_filt[first_scan] = P
    x_pred_arr[first_scan] = x
    P_pred_arr[first_scan] = P

    Q = np.diag([Q_FM, Q_LM, q_buf, q_tdee])

    for i in range(first_scan, n - 1):
        x_pred, P_pred = predict_step(x, P, calories[i], workouts_per_wk[i], Q,
                                      buf_decay, beta_buf, cal_ref)

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
            S = (H @ P @ H.T + r_weight)[0, 0]
            K = (P @ H.T) / S
            x = x + (K * innovation).flatten()
            P = P - np.outer(K.flatten(), H @ P)
            innovations[i + 1] = innovation
            innov_var[i + 1] = S

        # 2. FM scan: H = [1, 0, 0, 0]
        if not np.isnan(fm_scan[i + 1]):
            H = np.array([[1.0, 0.0, 0.0, 0.0]])
            y = fm_scan[i + 1]
            r_use = R_BOOT_FM if (i + 1 - first_scan) < 30 else R_FM_SCAN
            S = (H @ P @ H.T + r_use)[0, 0]
            K = (P @ H.T) / S
            fm_innov[i + 1], fm_var[i + 1] = y - (H @ x)[0], S
            x = x + (K * (y - (H @ x)[0])).flatten()
            P = P - np.outer(K.flatten(), H @ P)

        # 3. LM scan: H = [0, 1, 0, 0]
        if not np.isnan(lm_scan[i + 1]):
            H = np.array([[0.0, 1.0, 0.0, 0.0]])
            y = lm_scan[i + 1]
            r_use = R_BOOT_LM if (i + 1 - first_scan) < 30 else R_LM_SCAN
            S = (H @ P @ H.T + r_use)[0, 0]
            K = (P @ H.T) / S
            lm_innov[i + 1], lm_var[i + 1] = y - (H @ x)[0], S
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

    if collect is not None:
        collect["innov_var"] = innov_var
        collect["fm_innov"], collect["fm_var"] = fm_innov, fm_var
        collect["lm_innov"], collect["lm_var"] = lm_innov, lm_var
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

    # 5. Cumulative energy balance closure.
    #    Σ(intake − TDEE) must equal ΔFM·E_FAT + ΔLM·E_LEAN (the filter's own
    #    partition). The residual is the net energy injected by observation
    #    updates (weight / scan innovations pulling FM and LM directly) rather
    #    than through intake − TDEE. Expressed in lb of fat.
    #    Note (DH): this residual is NOT an identification of INTAKE_MULTIPLIER.
    #    TDEE is a free random walk and absorbs ~99% of any uniform intake
    #    rescaling; over 15 years a 1% multiplier change moves this number by
    #    ~0.25 lb while injecting ~115,000 cal. An earlier version omitted the
    #    ΔLM·E_LEAN term (−14.6 lb of lean over the record ≈ +4.6 lb of fat-
    #    equivalent), which is why "EB closes at m = 1.09" appeared exact.
    valid_cal = ~np.isnan(calories) & ~np.isnan(x_smooth[:, 3])
    idx = np.where(valid_cal)[0]
    cum_deficit = np.cumsum((calories - x_smooth[:, 3])[valid_cal])
    fm_change_observed = x_smooth[idx[-1], 0] - x_smooth[idx[0], 0]
    lm_change_observed = x_smooth[idx[-1], 1] - x_smooth[idx[0], 1]
    fm_change_implied = (cum_deficit[-1] - lm_change_observed * E_LEAN) / E_FAT
    out["eb_residual_lb"] = float(fm_change_implied - fm_change_observed)
    out["eb_residual_fm_only_lb"] = float(cum_deficit[-1] / E_FAT - fm_change_observed)

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


def run_one(dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr,
            workouts_per_wk, fm_start, lm_start, rmr_df, sigma_tdee, max_dtdee, temp_df=None,
            q_buf=None, r_weight=None, beta_buf=None, buf_hl=None):
    x_filt, P_filt, x_pred, P_pred, innov, first_scan = kalman_v4_forward(
        calories, weight_obs, fm_scan, lm_scan, expected_rmr,
        workouts_per_wk, fm_start, lm_start, sigma_tdee, max_dtdee,
        q_buf=q_buf, r_weight=r_weight, beta_buf=beta_buf, buf_hl=buf_hl,
    )
    x_smooth, P_smooth = rts_smoother_v4(x_filt, P_filt, x_pred, P_pred, first_scan, workouts_per_wk)
    diag = diagnostics(dates, x_smooth, innov, fm_scan, lm_scan, calories, rmr_df,
                       weight_obs, temp_df)
    return x_filt, x_smooth, innov, diag


def main():
    smooth, intake, strength, comp_scans, rmr_df, p2, temp_df = load_data()
    dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr, workouts_per_wk = \
        build_inputs(smooth, intake, strength, comp_scans, p2)

    first_scan = np.where(~np.isnan(fm_scan))[0][0]
    fm_start = fm_scan[first_scan]
    lm_start = lm_scan[first_scan]
    print(f"Intake model: {describe_intake_model()}")
    print(f"Filter tuning (from DM calibration): {describe_filter()}")
    print(f"Start: {str(dates[first_scan])[:10]}  FM={fm_start:.1f}  LM={lm_start:.1f}")
    print(f"Diagnostic sweep (does NOT select the canonical run): σ_TDEE ∈ {SIGMA_TDEE_GRID} × "
          f"max_dTDEE ∈ {MAX_DTDEE_GRID}, at the calibrated Q_BUF / R_WEIGHT")
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
    canon = (float(FILTER["sigma_tdee"]), int(FILTER["max_dtdee"]))
    grid = [(s_, c_) for s_ in SIGMA_TDEE_GRID for c_ in MAX_DTDEE_GRID]
    if canon not in grid:
        grid.append(canon)
    for sigma_tdee, max_dtdee in grid:
            x_filt, x_smooth, innov, diag = run_one(
                dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr,
                workouts_per_wk, fm_start, lm_start, rmr_df,
                sigma_tdee, max_dtdee, temp_df=temp_df,
                q_buf=FILTER["q_buf"], r_weight=FILTER["r_weight"],
                beta_buf=FILTER["beta_buf"], buf_hl=FILTER["buf_hl"],
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
            # The canonical run is the DM-calibrated tuning, selected by the
            # prediction-error likelihood over all assimilated observations
            # (weigh-ins + FM + LM scans) jointly with the intake model. No
            # composite score here: earlier versions ranked the sweep by scan
            # MAE + innovation ACF, and before that by terms that asserted an
            # activity factor (|TDEE/RMR − 1.15|) or rewarded a hypothesis the
            # filter is used to test (TDEE–temperature correlation).
            score = (diag['fm_scan_mae'] + diag['lm_scan_mae']
                     + max(diag['innov_acf1'], 0) * 5)
            if (float(sigma_tdee), int(max_dtdee)) == canon:
                best_score = score
                best_key = (sigma_tdee, max_dtdee)
                best_x = x_smooth
                best_diag = diag
                best_innov = innov

    print(f"\nCanonical: σ_TDEE={best_key[0]}, max_dTDEE={best_key[1]}, Q_BUF={FILTER['q_buf']:.3f}, "
          f"R_WEIGHT={FILTER['r_weight']:.2f}, BETA_BUF={FILTER['beta_buf']:.3f} lb/1000cal, BUF_HL={FILTER['buf_hl']:.1f} d "
          f"(from DM; scan-MAE+ACF composite {best_score:.2f} for reference)")
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
        # the intake series the filter actually consumed (downstream scripts use
        # intake_true rather than re-applying a multiplier)
        "intake_logged": intake["calories"].values,
        "intake_true": np.round(calories, 1),
    })
    out["date"] = pd.to_datetime(out["date"]).dt.strftime("%Y-%m-%d")
    path = ROOT / "analysis" / "P4_kalman_daily.csv"
    out.to_csv(path, index=False)
    print(f"Wrote {len(out)} rows to {path}")

    # Save sweep summary
    pd.DataFrame(results).to_csv(ROOT / "analysis" / "P4_kalman_sweep.csv", index=False)
    print(f"Wrote sweep results to analysis/P4_kalman_sweep.csv")
    save_result("P4_kalman_filter",
                sigma_tdee=best_key[0], max_dtdee=best_key[1], q_buf=FILTER["q_buf"], r_weight=FILTER["r_weight"],
                beta_buf=FILTER["beta_buf"], buf_hl=FILTER["buf_hl"],
                composite_score=best_score,
                intake_model=CALIBRATION["model"], intake_params=CALIBRATION["params"],
                effective_multiplier=INTAKE_MULTIPLIER,
                diagnostics={k: v for k, v in best_diag.items() if np.isscalar(v)},
                start_date=str(dates[first_scan])[:10], end_date=str(dates[-1])[:10])
    print("Wrote analysis/results/P4_kalman_filter.json")

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
