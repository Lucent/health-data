"""Parameter sweep for Kalman v2.

Evaluates TDEE dynamics and training-effect sensitivity against:
- composition scan fit
- innovation whiteness
- calorimetry consistency
- TDEE smoothness

Writes a stable CSV artifact for comparison.
"""

from pathlib import Path
import itertools
import numpy as np
import pandas as pd

from P4_kalman_filter_v2 import (
    ACTIVITY_FACTOR,
    TRAINING_DELTA_LBS,
    TRAINING_HALF_LIFE,
    build_inputs,
    kalman_forward,
    load_data,
    rts_smoother,
)

ROOT = Path(__file__).resolve().parent.parent


def evaluate_candidate(q_tdee, mean_revert, training_delta_lbs, training_half_life):
    smooth, comp, intake, strength, comp_scans, rmr = load_data()
    dates, calories, fm_obs, expected_rmr, known_lean = build_inputs(
        smooth,
        comp,
        intake,
        strength,
        comp_scans,
        training_delta_lbs=training_delta_lbs,
        training_half_life=training_half_life,
    )

    first_obs = np.where(~np.isnan(fm_obs))[0][0]
    fm_start = fm_obs[first_obs]
    x_filt, P_filt, x_pred, P_pred, innovations = kalman_forward(
        calories,
        fm_obs,
        expected_rmr,
        fm_start,
        q_tdee=q_tdee,
        mean_revert=mean_revert,
    )
    x_smooth, P_smooth = rts_smoother(
        x_filt,
        P_filt,
        x_pred,
        P_pred,
        calories,
        expected_rmr,
        mean_revert=mean_revert,
    )

    valid_innov = innovations[~np.isnan(innovations)]
    centered = valid_innov - valid_innov.mean()
    innov_acf1 = np.corrcoef(centered[:-1], centered[1:])[0, 1]
    innov_std = valid_innov.std()

    errors = []
    for _, row in comp_scans.iterrows():
        idx = np.searchsorted(dates, np.datetime64(row["date"]))
        if idx < len(x_smooth) and not np.isnan(x_smooth[idx, 0]):
            errors.append(x_smooth[idx, 0] - row["fat_mass_lbs"])
    errors = np.array(errors)
    fm_mae = np.mean(np.abs(errors))
    fm_rmse = np.sqrt(np.mean(errors ** 2))
    fm_max_abs = np.max(np.abs(errors))

    cal_errors = []
    for _, row in rmr.iterrows():
        idx = np.searchsorted(dates, np.datetime64(row["date"]))
        if idx < len(x_smooth) and not np.isnan(x_smooth[idx, 1]):
            cal_errors.append(abs(x_smooth[idx, 1] / row["rmr_kcal"] - ACTIVITY_FACTOR))
    cal_ratio_err = np.mean(cal_errors)

    tdee_valid = x_smooth[~np.isnan(x_smooth[:, 1]), 1]
    tdee_dd_std = np.std(np.diff(tdee_valid))

    score = abs(innov_acf1) + cal_ratio_err * 2 + fm_mae * 0.05

    return {
        "q_tdee": q_tdee,
        "mean_revert": mean_revert,
        "training_delta_lbs": training_delta_lbs,
        "training_half_life": training_half_life,
        "fm_mae": round(fm_mae, 4),
        "fm_rmse": round(fm_rmse, 4),
        "fm_max_abs": round(fm_max_abs, 4),
        "innov_std": round(innov_std, 4),
        "innov_acf1": round(innov_acf1, 4),
        "cal_ratio_err": round(cal_ratio_err, 4),
        "tdee_dd_std": round(tdee_dd_std, 4),
        "score": round(score, 4),
    }


def main():
    q_values = [25, 50, 100, 200, 500, 1000]
    mean_revert_values = [0.002, 0.005, 0.01, 0.02, 0.05]
    training_delta_values = [0.0, round(TRAINING_DELTA_LBS / 2, 3), TRAINING_DELTA_LBS]
    half_life_values = [TRAINING_HALF_LIFE]

    rows = []
    for q_tdee, mean_revert, training_delta_lbs, training_half_life in itertools.product(
        q_values,
        mean_revert_values,
        training_delta_values,
        half_life_values,
    ):
        rows.append(
            evaluate_candidate(
                q_tdee=q_tdee,
                mean_revert=mean_revert,
                training_delta_lbs=training_delta_lbs,
                training_half_life=training_half_life,
            )
        )

    out = pd.DataFrame(rows).sort_values(
        ["score", "fm_mae", "innov_acf1", "fm_rmse"],
        ascending=[True, True, True, True],
    )
    path = ROOT / "analysis" / "P4_kalman_v2_sweep.csv"
    out.to_csv(path, index=False)

    print("=== Best overall (score) ===")
    print(out.head(10).to_string(index=False))

    print("\n=== Best composition fit ===")
    print(out.sort_values(["fm_mae", "fm_rmse", "innov_acf1"]).head(10).to_string(index=False))

    print("\n=== Best innovation whiteness ===")
    print(out.sort_values(["innov_acf1", "fm_mae", "fm_rmse"]).head(10).to_string(index=False))

    print(f"\nWrote {len(out)} rows to {path}")


if __name__ == "__main__":
    main()
