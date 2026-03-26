"""Kalman filter v2: two-state (fat mass, TDEE) with known lean mass input.

Lean mass is a known input from composition interpolation + training model
(finding AA: 0.091 lbs/session, 275-day half-life). This gives the filter
a direct fat mass observation: weight - known_lean = fat + noise.

Fat mass changes through energy balance. TDEE drifts slowly with
mean-reversion toward composition-aware expected RMR.

Constants at the top for sweep.
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analysis"))

# ── Constants ────────────────────────────────────────────────────────
CAL_PER_LB = 3500
FORBES_C_KG = 10.4
LBS_TO_KG = 0.453592

# Process noise
Q_FAT = 0.005              # lbs² per day — small but nonzero (intake undercount)
Q_TDEE = 500

# Observation noise
R_OBS = 0.97               # lbs² — measured from consecutive-day variance after corrections

# Mean reversion
MEAN_REVERT_RATE = 0.005
ACTIVITY_FACTOR = 1.15

# Training model (finding AA)
# Default is disabled: the sweep favored no added workout-driven lean term.
TRAINING_DELTA_LBS = 0.0
TRAINING_HALF_LIFE = 275
TRAINING_DECAY = np.log(2) / TRAINING_HALF_LIFE


def load_data():
    smooth = pd.read_csv(ROOT / "analysis" / "P1_smoothed_weight.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    strength = pd.read_csv(ROOT / "workout" / "strength.csv", parse_dates=["date"])
    comp_scans = pd.read_csv(ROOT / "composition" / "composition.csv", parse_dates=["date"])
    rmr = pd.read_csv(ROOT / "RMR" / "rmr.csv")
    return smooth, comp, intake, strength, comp_scans, rmr


def compute_known_lean(intake_dates, comp_ffm, strength, comp_scans,
                       training_delta_lbs=TRAINING_DELTA_LBS,
                       training_half_life=TRAINING_HALF_LIFE):
    """Lean mass = composition interpolation + training-effect dynamics."""
    dates = intake_dates
    n = len(dates)
    workout_days = strength["date"].values.astype("datetime64[D]")
    training_decay = np.log(2) / training_half_life

    # Training effect at each day
    training = np.zeros(n)
    for i, d in enumerate(dates):
        day = np.datetime64(d, "D")
        days_since = (day - workout_days).astype(float)
        past = days_since[days_since > 0]
        training[i] = training_delta_lbs * np.sum(np.exp(-training_decay * past))

    # P2 already interpolates scan-to-scan lean mass, which implicitly includes
    # whatever training effect existed at the scan anchors. To avoid shifting
    # the anchors, remove the linearly interpolated scan-era training effect
    # instead of subtracting a single global offset.
    date_days = dates.astype("datetime64[D]").astype(np.int64)
    scan_days = comp_scans["date"].values.astype("datetime64[D]").astype(np.int64)
    scan_effects = np.interp(scan_days, date_days, training)
    baseline_training = np.interp(date_days, scan_days, scan_effects)

    return comp_ffm + (training - baseline_training)


def build_inputs(smooth, comp, intake, strength, comp_scans,
                 training_delta_lbs=TRAINING_DELTA_LBS,
                 training_half_life=TRAINING_HALF_LIFE):
    dates = intake["date"].values
    calories = intake["calories"].values
    n = len(dates)

    expected_rmr = comp.set_index("date")["expected_rmr"].reindex(intake["date"]).values
    comp_ffm = comp.set_index("date")["ffm_lbs"].reindex(intake["date"]).values

    # Known lean mass (composition + training)
    known_lean = compute_known_lean(
        dates,
        comp_ffm,
        strength,
        comp_scans,
        training_delta_lbs=training_delta_lbs,
        training_half_life=training_half_life,
    )

    # Weight observations (glycogen+sodium corrected)
    obs_map = smooth.set_index("date")["smoothed_weight_lbs"]
    weight_obs = np.full(n, np.nan)
    for i, d in enumerate(dates):
        if d in obs_map.index:
            v = obs_map[d]
            if not np.isnan(v):
                weight_obs[i] = v

    # Fat mass observation = weight - known_lean
    fm_obs = weight_obs - known_lean

    return dates, calories, fm_obs, expected_rmr, known_lean


def forbes_fat_fraction(fm_lbs):
    fm_kg = max(fm_lbs, 5.0) * LBS_TO_KG
    return 1.0 - (FORBES_C_KG ** 2) / (FORBES_C_KG + fm_kg) ** 2


def kalman_forward(calories, fm_obs, expected_rmr, fm_start,
                   q_tdee=Q_TDEE, mean_revert=MEAN_REVERT_RATE):
    n = len(calories)

    H = np.array([[1.0, 0.0]])  # observe fat mass directly
    Q = np.diag([Q_FAT, q_tdee])
    R_mat = np.array([[R_OBS]])

    first_obs = np.where(~np.isnan(fm_obs))[0][0]
    init_tdee = expected_rmr[first_obs] * ACTIVITY_FACTOR if not np.isnan(expected_rmr[first_obs]) else 2100.0

    x = np.array([fm_start, init_tdee])
    P = np.diag([25.0, 250000.0])

    x_filt = np.full((n, 2), np.nan)
    P_filt = np.full((n, 2, 2), np.nan)
    x_pred_arr = np.full((n, 2), np.nan)
    P_pred_arr = np.full((n, 2, 2), np.nan)
    innovations = np.full(n, np.nan)

    x_filt[first_obs] = x
    P_filt[first_obs] = P
    x_pred_arr[first_obs] = x
    P_pred_arr[first_obs] = P

    for i in range(first_obs, n - 1):
        fm, tdee = x
        ff = forbes_fat_fraction(fm)
        surplus = calories[i] - tdee
        target_tdee = expected_rmr[i] * ACTIVITY_FACTOR if not np.isnan(expected_rmr[i]) else tdee

        x_pred = np.array([
            fm + surplus * ff / CAL_PER_LB,
            tdee + mean_revert * (target_tdee - tdee),
        ])

        F = np.array([
            [1.0, -ff / CAL_PER_LB],
            [0.0, 1.0 - mean_revert],
        ])
        P_pred = F @ P @ F.T + Q

        x_pred_arr[i + 1] = x_pred
        P_pred_arr[i + 1] = P_pred

        if not np.isnan(fm_obs[i + 1]):
            y = fm_obs[i + 1]
            y_pred = (H @ x_pred)[0]
            innovation = y - y_pred
            S = (H @ P_pred @ H.T + R_mat)[0, 0]
            K = (P_pred @ H.T) / S
            x = x_pred + (K * innovation).flatten()
            P = P_pred - np.outer(K.flatten(), H @ P_pred)
            innovations[i + 1] = innovation
        else:
            x = x_pred
            P = P_pred

        x[0] = max(x[0], 0.0)
        x_filt[i + 1] = x
        P_filt[i + 1] = P

    return x_filt, P_filt, x_pred_arr, P_pred_arr, innovations


def rts_smoother(x_filt, P_filt, x_pred, P_pred, calories, expected_rmr,
                 mean_revert=MEAN_REVERT_RATE):
    n = len(x_filt)
    x_smooth = np.copy(x_filt)
    P_smooth = np.copy(P_filt)

    valid = ~np.isnan(x_filt[:, 0])
    last_valid = np.where(valid)[0][-1]

    for i in range(last_valid - 1, -1, -1):
        if np.isnan(x_filt[i, 0]) or np.isnan(P_pred[i + 1, 0, 0]):
            continue
        fm = x_filt[i, 0]
        ff = forbes_fat_fraction(fm)
        F = np.array([
            [1.0, -ff / CAL_PER_LB],
            [0.0, 1.0 - mean_revert],
        ])
        P_pred_inv = np.linalg.inv(P_pred[i + 1])
        G = P_filt[i] @ F.T @ P_pred_inv
        x_smooth[i] = x_filt[i] + G @ (x_smooth[i + 1] - x_pred[i + 1])
        P_smooth[i] = P_filt[i] + G @ (P_smooth[i + 1] - P_pred[i + 1]) @ G.T

    return x_smooth, P_smooth


def main():
    smooth, comp, intake, strength, comp_scans, rmr = load_data()
    dates, calories, fm_obs, expected_rmr, known_lean = \
        build_inputs(smooth, comp, intake, strength, comp_scans)

    first_obs = np.where(~np.isnan(fm_obs))[0][0]
    fm_start = fm_obs[first_obs]
    print(f"Start: {str(dates[first_obs])[:10]}  FM={fm_start:.0f}  "
          f"lean={known_lean[first_obs]:.0f}")
    print(f"Constants: Q_FAT={Q_FAT}  Q_TDEE={Q_TDEE}  R={R_OBS}  "
          f"MEAN_REVERT={MEAN_REVERT_RATE}  ACTIVITY={ACTIVITY_FACTOR}")

    x_filt, P_filt, x_pred, P_pred, innovations = kalman_forward(
        calories, fm_obs, expected_rmr, fm_start)
    x_smooth, P_smooth = rts_smoother(
        x_filt, P_filt, x_pred, P_pred, calories, expected_rmr)

    # Validate
    valid_innov = innovations[~np.isnan(innovations)]
    centered = valid_innov - valid_innov.mean()
    acf1 = np.corrcoef(centered[:-1], centered[1:])[0, 1]
    print(f"\nInnovations: n={len(valid_innov)}  std={valid_innov.std():.3f}  "
          f"lag1_acf={acf1:.4f}")

    # Composition fit
    print(f"\nComposition fit:")
    errors = []
    for _, row in comp_scans.iterrows():
        idx = np.searchsorted(dates, np.datetime64(row["date"]))
        if idx < len(x_smooth) and not np.isnan(x_smooth[idx, 0]):
            err = x_smooth[idx, 0] - row["fat_mass_lbs"]
            errors.append(err)
            if abs(err) > 5:
                print(f"  {str(row['date'])[:10]}: filter={x_smooth[idx,0]:.1f}  "
                      f"meas={row['fat_mass_lbs']:.1f}  err={err:+.1f}")
    print(f"  MAE: {np.mean(np.abs(errors)):.1f} lbs  n={len(errors)}")

    # Calorimetry
    print(f"\nCalorimetry:")
    for _, row in rmr.iterrows():
        idx = np.searchsorted(dates, np.datetime64(row["date"]))
        if idx < len(x_smooth) and not np.isnan(x_smooth[idx, 1]):
            tdee = x_smooth[idx, 1]
            print(f"  {row['date']}: RMR={row['rmr_kcal']}  TDEE={tdee:.0f}  "
                  f"ratio={tdee/row['rmr_kcal']:.2f}")

    tdee_valid = x_smooth[~np.isnan(x_smooth[:, 1]), 1]
    print(f"\nTDEE day-to-day std: {np.std(np.diff(tdee_valid)):.1f} cal/day")

    # Save
    out = pd.DataFrame({
        "date": dates,
        "fat_mass_lbs_filtered": np.round(x_filt[:, 0], 2),
        "fat_mass_std_filtered": np.round(np.sqrt(np.maximum(P_filt[:, 0, 0], 0)), 2),
        "tdee_filtered": np.round(x_filt[:, 1], 0),
        "tdee_std_filtered": np.round(np.sqrt(np.maximum(P_filt[:, 1, 1], 0)), 0),
        "fat_mass_lbs": np.round(x_smooth[:, 0], 2),
        "fat_mass_std": np.round(np.sqrt(np.maximum(P_smooth[:, 0, 0], 0)), 2),
        "tdee": np.round(x_smooth[:, 1], 0),
        "tdee_std": np.round(np.sqrt(np.maximum(P_smooth[:, 1, 1], 0)), 0),
        "innovation": np.round(innovations, 3),
    })
    out["date"] = pd.to_datetime(out["date"]).dt.strftime("%Y-%m-%d")
    path = ROOT / "analysis" / "P4_kalman_daily.csv"
    out.to_csv(path, index=False)
    print(f"\nWrote {len(out)} rows to {path}")

    print(f"\n=== Summary ===")
    print(f"  Innovation ACF: {acf1:.4f} (target < 0.2)")
    print(f"  FM MAE at scans: {np.mean(np.abs(errors)):.1f} lbs")
    print(f"  TDEE smoothness: {np.std(np.diff(tdee_valid)):.1f} cal/day")


if __name__ == "__main__":
    main()
