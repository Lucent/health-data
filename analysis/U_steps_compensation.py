"""Test compensated and state-dependent step effects.

This revisits the weak daily gravitostat result with richer questions:
1. Do recent steps predict future intake or future TDEE/RMR over 7-14 days?
2. Does the answer depend on rising/stable/falling branch state?
3. Do high recent-step periods look protective or compensated?

Outputs:
    analysis/steps_compensation_regression.csv
    analysis/steps_compensation_phase_thresholds.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
TREND_WINDOW_DAYS = 90
TREND_THRESHOLD_LBS = 3.0


def load_data():
    steps = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    daily = intake.merge(steps[["date", "steps"]], on="date", how="left")
    daily = daily.merge(kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="left")
    daily = daily.merge(comp[["date", "weight_lbs", "expected_rmr"]], on="date", how="left")
    daily = daily.merge(tirz[["date", "effective_level"]], on="date", how="left")
    daily["effective_level"] = daily["effective_level"].fillna(0)
    daily = daily[daily["effective_level"] == 0].sort_values("date").reset_index(drop=True)
    daily["tdee_rmr_ratio"] = daily["tdee"] / daily["expected_rmr"]
    daily["foot_pounds"] = daily["weight_lbs"] * daily["steps"]
    daily["fat_delta_90d"] = daily["fat_mass_lbs"].diff(TREND_WINDOW_DAYS)
    daily["phase"] = np.where(
        daily["fat_delta_90d"] <= -TREND_THRESHOLD_LBS,
        "falling",
        np.where(daily["fat_delta_90d"] >= TREND_THRESHOLD_LBS, "rising", "stable"),
    )
    daily["steps_7d"] = daily["steps"].rolling(7, min_periods=5).mean().shift(1)
    daily["foot_pounds_7d"] = daily["foot_pounds"].rolling(7, min_periods=5).mean().shift(1)
    for horizon in [7, 14]:
        daily[f"future{horizon}_calories"] = daily["calories"].shift(-1).rolling(horizon, min_periods=horizon).mean()
        daily[f"future{horizon}_ratio"] = daily["tdee_rmr_ratio"].shift(-1).rolling(horizon, min_periods=horizon).mean()
    return daily


def standardized_beta(df, y_col, x_col, covariates, subset=None):
    if subset is not None:
        df = df[subset].copy()
    df = df.dropna(subset=[y_col, x_col] + covariates).copy()
    y = df[y_col].to_numpy(float)
    X = np.column_stack([df[col].to_numpy(float) for col in [x_col] + covariates])
    mean = X.mean(axis=0)
    std = X.std(axis=0)
    std[std == 0] = 1.0
    X = (X - mean) / std
    y_std = y.std()
    y = (y - y.mean()) / (y_std if y_std != 0 else 1.0)
    beta = np.linalg.lstsq(np.column_stack([np.ones(len(X)), X]), y, rcond=None)[0]
    return len(df), beta[1]


def regression_summary(daily):
    rows = []
    for metric in ["steps_7d", "foot_pounds_7d"]:
        for target, covars in [
            ("future7_calories", ["calories"]),
            ("future14_calories", ["calories"]),
            ("future7_ratio", ["tdee_rmr_ratio"]),
            ("future14_ratio", ["tdee_rmr_ratio"]),
        ]:
            n_rows, beta = standardized_beta(daily, target, metric, covars)
            rows.append(
                {
                    "subset": "all",
                    "metric": metric,
                    "target": target,
                    "n_rows": n_rows,
                    "std_beta": round(beta, 4),
                }
            )

        for phase in ["falling", "stable", "rising"]:
            n_rows, beta = standardized_beta(
                daily,
                "future14_ratio",
                metric,
                ["tdee_rmr_ratio"],
                subset=daily["phase"] == phase,
            )
            rows.append(
                {
                    "subset": phase,
                    "metric": metric,
                    "target": "future14_ratio",
                    "n_rows": n_rows,
                    "std_beta": round(beta, 4),
                }
            )
    return pd.DataFrame(rows)


def phase_threshold_summary(daily):
    rows = []
    for phase in ["falling", "stable", "rising"]:
        phase_df = daily[daily["phase"] == phase].dropna(subset=["steps_7d", "future14_calories", "future14_ratio"])
        if len(phase_df) < 60:
            continue
        for quantile in [0.5, 0.8]:
            threshold = phase_df["steps_7d"].quantile(quantile)
            hi = phase_df[phase_df["steps_7d"] >= threshold]
            lo = phase_df[phase_df["steps_7d"] < threshold]
            rows.append(
                {
                    "phase": phase,
                    "quantile": quantile,
                    "steps_7d_threshold": round(threshold, 1),
                    "high_n": len(hi),
                    "low_n": len(lo),
                    "high_future14_calories": round(hi["future14_calories"].mean(), 1),
                    "low_future14_calories": round(lo["future14_calories"].mean(), 1),
                    "high_future14_ratio": round(hi["future14_ratio"].mean(), 4),
                    "low_future14_ratio": round(lo["future14_ratio"].mean(), 4),
                    "high_same_day_calories": round(hi["calories"].mean(), 1),
                    "low_same_day_calories": round(lo["calories"].mean(), 1),
                }
            )
    return pd.DataFrame(rows)


def save_outputs(reg_df, threshold_df):
    reg_df.to_csv(ROOT / "analysis" / "U_steps_compensation_regression.csv", index=False)
    threshold_df.to_csv(ROOT / "analysis" / "U_steps_compensation_phase_thresholds.csv", index=False)


def print_report(reg_df, threshold_df):
    print("\n=== Steps Compensation ===")
    for _, row in reg_df[reg_df["subset"] == "all"].iterrows():
        print(
            f"{row['metric']:>16} -> {row['target']:>18}: beta={row['std_beta']:+.4f}  n={int(row['n_rows'])}"
        )

    print("\nFuture14 TDEE/RMR by phase:")
    for _, row in reg_df[(reg_df["target"] == "future14_ratio") & (reg_df["subset"] != "all")].iterrows():
        print(
            f"  {row['subset']:>7}  {row['metric']:>16}: beta={row['std_beta']:+.4f}  n={int(row['n_rows'])}"
        )

    print("\nPhase threshold summary (steps_7d):")
    for _, row in threshold_df.iterrows():
        print(
            f"  {row['phase']:>7} q={row['quantile']:.1f} thr={row['steps_7d_threshold']:.0f}: "
            f"future14 cal {row['high_future14_calories']:.0f} vs {row['low_future14_calories']:.0f} | "
            f"future14 ratio {row['high_future14_ratio']:.4f} vs {row['low_future14_ratio']:.4f}"
        )


def main():
    daily = load_data()
    reg_df = regression_summary(daily)
    threshold_df = phase_threshold_summary(daily)
    save_outputs(reg_df, threshold_df)
    print_report(reg_df, threshold_df)


if __name__ == "__main__":
    main()
