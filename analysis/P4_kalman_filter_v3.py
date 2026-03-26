"""Kalman filter v3: three-state (fat mass, TDEE, intake bias).

Lean mass is a known input (composition + training model).
Fat mass changes through energy balance with bias-corrected intake.
TDEE drifts slowly. Intake bias is a slow-drifting multiplier estimated
from the data, anchored by 25 Cosmed RMR measurements that constrain TDEE
independently of intake.

State: [fat_mass_lbs, tdee_cal, intake_bias]
  true_intake = logged_calories × (1 + bias)
  fat_mass += (true_intake - tdee) × fat_fraction / CAL_PER_LB
  tdee += mean_reversion toward expected_rmr × activity_factor
  bias drifts slowly (constrained by Cosmed observations)

Observations:
  Type 1 (weight days): smoothed_weight - known_lean = fat_mass + noise
  Type 2 (Cosmed days): measured_rmr = tdee / activity_factor + noise
    (This pins TDEE, forcing the filter to put residual error into bias)

Constants at the top for sweep.
"""

from pathlib import Path
import importlib
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
Q_FAT = 0.005              # lbs²/day
Q_TDEE = 200               # cal²/day (~14 cal/day std)
Q_BIAS = 0.000001          # bias²/day — very slow drift (~0.001/day = ~0.4%/year)

# Observation noise
R_WEIGHT = 0.97            # lbs² — gut, hydration, scale
R_RMR = 15000              # cal² — Cosmed measurement noise (~122 cal std, RMSE=170)

# Mean reversion
MEAN_REVERT_RATE = 0.005
ACTIVITY_FACTOR = 1.15

# Training model (finding AA)
TRAINING_DELTA_LBS = 0.091
TRAINING_HALF_LIFE = 275
TRAINING_DECAY = np.log(2) / TRAINING_HALF_LIFE


def load_data():
    smooth = pd.read_csv(ROOT / "analysis" / "P1_smoothed_weight.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    strength = pd.read_csv(ROOT / "workout" / "strength.csv", parse_dates=["date"])
    comp_scans = pd.read_csv(ROOT / "composition" / "composition.csv", parse_dates=["date"])
    rmr = pd.read_csv(ROOT / "RMR" / "rmr.csv", parse_dates=["date"])
    return smooth, comp, intake, strength, comp_scans, rmr


def compute_known_lean(intake_dates, comp_ffm, strength, comp_scans):
    dates = intake_dates
    n = len(dates)
    workout_days = strength["date"].values.astype("datetime64[D]")
    training = np.zeros(n)
    for i, d in enumerate(dates):
        day = np.datetime64(d, "D")
        days_since = (day - workout_days).astype(float)
        past = days_since[days_since > 0]
        training[i] = TRAINING_DELTA_LBS * np.sum(np.exp(-TRAINING_DECAY * past))
    scan_effects = []
    for d in comp_scans["date"].values:
        idx = np.searchsorted(dates, np.datetime64(d))
        if idx < n:
            scan_effects.append(training[idx])
    median_scan = np.median(scan_effects) if scan_effects else 0.0
    return comp_ffm + (training - median_scan)


def build_inputs(smooth, comp, intake, strength, comp_scans, rmr):
    dates = intake["date"].values
    calories = intake["calories"].values
    n = len(dates)

    expected_rmr = comp.set_index("date")["expected_rmr"].reindex(intake["date"]).values
    comp_ffm = comp.set_index("date")["ffm_lbs"].reindex(intake["date"]).values
    known_lean = compute_known_lean(dates, comp_ffm, strength, comp_scans)

    # Weight observations
    obs_map = smooth.set_index("date")["smoothed_weight_lbs"]
    weight_obs = np.full(n, np.nan)
    for i, d in enumerate(dates):
        if d in obs_map.index:
            v = obs_map[d]
            if not np.isnan(v):
                weight_obs[i] = v
    fm_obs = weight_obs - known_lean

    # RMR observations (Cosmed measurements)
    rmr_obs = np.full(n, np.nan)
    for _, row in rmr.iterrows():
        idx = np.searchsorted(dates, np.datetime64(row["date"]))
        if idx < n:
            rmr_obs[idx] = row["rmr_kcal"]

    return dates, calories, fm_obs, rmr_obs, expected_rmr, known_lean


def forbes_fat_fraction(fm_lbs):
    fm_kg = max(fm_lbs, 5.0) * LBS_TO_KG
    return 1.0 - (FORBES_C_KG ** 2) / (FORBES_C_KG + fm_kg) ** 2


def kalman_forward(calories, fm_obs, rmr_obs, expected_rmr, fm_start):
    n = len(calories)

    H_weight = np.array([[1.0, 0.0, 0.0]])  # observe fat mass
    H_rmr = np.array([[0.0, 1.0 / ACTIVITY_FACTOR, 0.0]])  # observe TDEE/activity = RMR
    Q = np.diag([Q_FAT, Q_TDEE, Q_BIAS])

    first_obs = np.where(~np.isnan(fm_obs))[0][0]
    init_tdee = expected_rmr[first_obs] * ACTIVITY_FACTOR if not np.isnan(expected_rmr[first_obs]) else 2100.0

    # State: [fat_mass, tdee, intake_bias]
    x = np.array([fm_start, init_tdee, 0.0])  # bias starts at 0 (no correction)
    P = np.diag([25.0, 250000.0, 0.01])  # bias uncertainty ±10%

    x_filt = np.full((n, 3), np.nan)
    P_filt = np.full((n, 3, 3), np.nan)
    x_pred_arr = np.full((n, 3), np.nan)
    P_pred_arr = np.full((n, 3, 3), np.nan)
    innovations = np.full(n, np.nan)

    x_filt[first_obs] = x
    P_filt[first_obs] = P
    x_pred_arr[first_obs] = x
    P_pred_arr[first_obs] = P

    for i in range(first_obs, n - 1):
        fm, tdee, bias = x
        ff = forbes_fat_fraction(fm)
        true_intake = calories[i] * (1.0 + bias)
        surplus = true_intake - tdee
        target_tdee = expected_rmr[i] * ACTIVITY_FACTOR if not np.isnan(expected_rmr[i]) else tdee

        x_pred = np.array([
            fm + surplus * ff / CAL_PER_LB,
            tdee + MEAN_REVERT_RATE * (target_tdee - tdee),
            bias,  # bias random walk (drift only through Q)
        ])

        # Jacobian
        F = np.array([
            [1.0, -ff / CAL_PER_LB, calories[i] * ff / CAL_PER_LB],
            [0.0, 1.0 - MEAN_REVERT_RATE, 0.0],
            [0.0, 0.0, 1.0],
        ])
        P_pred = F @ P @ F.T + Q

        x_pred_arr[i + 1] = x_pred
        P_pred_arr[i + 1] = P_pred

        x = x_pred.copy()
        P = P_pred.copy()

        # Update 1: weight observation
        if not np.isnan(fm_obs[i + 1]):
            y = fm_obs[i + 1]
            y_pred = (H_weight @ x)[0]
            innovation = y - y_pred
            S = (H_weight @ P @ H_weight.T)[0, 0] + R_WEIGHT
            K = (P @ H_weight.T) / S
            x = x + (K * innovation).flatten()
            P = P - np.outer(K.flatten(), H_weight @ P)
            innovations[i + 1] = innovation

        # Update 2: RMR observation (on Cosmed days)
        if not np.isnan(rmr_obs[i + 1]):
            y_rmr = rmr_obs[i + 1]
            y_rmr_pred = (H_rmr @ x)[0]
            innovation_rmr = y_rmr - y_rmr_pred
            S_rmr = (H_rmr @ P @ H_rmr.T)[0, 0] + R_RMR
            K_rmr = (P @ H_rmr.T) / S_rmr
            x = x + (K_rmr * innovation_rmr).flatten()
            P = P - np.outer(K_rmr.flatten(), H_rmr @ P)

        x[0] = max(x[0], 0.0)
        x[2] = np.clip(x[2], -0.3, 0.3)  # bias bounded ±30%
        x_filt[i + 1] = x
        P_filt[i + 1] = P

    return x_filt, P_filt, x_pred_arr, P_pred_arr, innovations


def rts_smoother(x_filt, P_filt, x_pred, P_pred, calories, expected_rmr):
    n = len(x_filt)
    x_smooth = np.copy(x_filt)
    P_smooth = np.copy(P_filt)

    valid = ~np.isnan(x_filt[:, 0])
    last_valid = np.where(valid)[0][-1]

    for i in range(last_valid - 1, -1, -1):
        if np.isnan(x_filt[i, 0]) or np.isnan(P_pred[i + 1, 0, 0]):
            continue
        fm, tdee, bias = x_filt[i]
        ff = forbes_fat_fraction(fm)
        F = np.array([
            [1.0, -ff / CAL_PER_LB, calories[i] * ff / CAL_PER_LB],
            [0.0, 1.0 - MEAN_REVERT_RATE, 0.0],
            [0.0, 0.0, 1.0],
        ])
        P_pred_inv = np.linalg.inv(P_pred[i + 1])
        G = P_filt[i] @ F.T @ P_pred_inv
        x_smooth[i] = x_filt[i] + G @ (x_smooth[i + 1] - x_pred[i + 1])
        P_smooth[i] = P_filt[i] + G @ (P_smooth[i + 1] - P_pred[i + 1]) @ G.T

    return x_smooth, P_smooth


def main():
    smooth, comp, intake, strength, comp_scans, rmr = load_data()
    dates, calories, fm_obs, rmr_obs, expected_rmr, known_lean = \
        build_inputs(smooth, comp, intake, strength, comp_scans, rmr)

    first_obs = np.where(~np.isnan(fm_obs))[0][0]
    fm_start = fm_obs[first_obs]
    n_rmr = (~np.isnan(rmr_obs)).sum()

    print(f"Start: {str(dates[first_obs])[:10]}  FM={fm_start:.0f}")
    print(f"Constants: Q_FAT={Q_FAT}  Q_TDEE={Q_TDEE}  Q_BIAS={Q_BIAS}  "
          f"R_WEIGHT={R_WEIGHT}  R_RMR={R_RMR}")
    print(f"Observations: {(~np.isnan(fm_obs)).sum()} weight, {n_rmr} RMR")

    print("\nRunning forward filter...")
    x_filt, P_filt, x_pred, P_pred, innovations = kalman_forward(
        calories, fm_obs, rmr_obs, expected_rmr, fm_start)

    print("Running RTS smoother...")
    x_smooth, P_smooth = rts_smoother(
        x_filt, P_filt, x_pred, P_pred, calories, expected_rmr)

    # Validate
    valid_innov = innovations[~np.isnan(innovations)]
    centered = valid_innov - valid_innov.mean()
    acf1 = np.corrcoef(centered[:-1], centered[1:])[0, 1]
    print(f"\nInnovations: n={len(valid_innov)}  std={valid_innov.std():.3f}  "
          f"lag1_acf={acf1:.4f}")

    # Composition fit
    errors = []
    for _, row in comp_scans.iterrows():
        idx = np.searchsorted(dates, np.datetime64(row["date"]))
        if idx < len(x_smooth) and not np.isnan(x_smooth[idx, 0]):
            err = x_smooth[idx, 0] - row["fat_mass_lbs"]
            errors.append(err)
            if abs(err) > 5:
                print(f"  {str(row['date'])[:10]}: filter={x_smooth[idx,0]:.1f}  "
                      f"meas={row['fat_mass_lbs']:.1f}  err={err:+.1f}")
    print(f"  FM MAE: {np.mean(np.abs(errors)):.1f} lbs  n={len(errors)}")

    # Calorimetry
    print(f"\nCalorimetry:")
    for _, row in rmr.iterrows():
        idx = np.searchsorted(dates, np.datetime64(row["date"]))
        if idx < len(x_smooth) and not np.isnan(x_smooth[idx, 1]):
            tdee = x_smooth[idx, 1]
            bias = x_smooth[idx, 2]
            print(f"  {str(row['date'])[:10]}: RMR={row['rmr_kcal']}  TDEE={tdee:.0f}  "
                  f"ratio={tdee/row['rmr_kcal']:.2f}  bias={bias:+.3f} ({bias*100:+.1f}%)")

    # Bias trajectory
    print(f"\nIntake bias by year (positive = eating more than logged):")
    intake_df = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    for year in range(2011, 2027):
        mask = intake_df["date"].dt.year == year
        idx_range = np.where(mask.values)[0]
        if len(idx_range) > 30:
            biases = x_smooth[idx_range, 2]
            valid_b = biases[~np.isnan(biases)]
            if len(valid_b) > 0:
                print(f"  {year}: bias={valid_b.mean():+.3f} ({valid_b.mean()*100:+.1f}%)  "
                      f"= {valid_b.mean() * intake_df.loc[mask, 'calories'].mean():.0f} cal/day")

    tdee_valid = x_smooth[~np.isnan(x_smooth[:, 1]), 1]
    print(f"\nTDEE day-to-day std: {np.std(np.diff(tdee_valid)):.1f} cal/day")

    # Save (same format as v2, with bias as extra column)
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
        "intake_bias": np.round(x_smooth[:, 2], 4),
        "innovation": np.round(innovations, 3),
    })
    out["date"] = pd.to_datetime(out["date"]).dt.strftime("%Y-%m-%d")
    path = ROOT / "analysis" / "P4_kalman_daily_v3.csv"
    out.to_csv(path, index=False)
    print(f"\nWrote {len(out)} rows to {path}")

    print(f"\n=== Summary ===")
    print(f"  Innovation ACF: {acf1:.4f}")
    print(f"  FM MAE: {np.mean(np.abs(errors)):.1f} lbs")
    print(f"  TDEE std: {np.std(np.diff(tdee_valid)):.1f} cal/day")
    bias_all = x_smooth[~np.isnan(x_smooth[:, 2]), 2]
    print(f"  Bias range: {bias_all.min():+.3f} to {bias_all.max():+.3f} "
          f"({bias_all.min()*100:+.1f}% to {bias_all.max()*100:+.1f}%)")


if __name__ == "__main__":
    main()
