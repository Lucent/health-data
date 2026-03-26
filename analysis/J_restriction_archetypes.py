"""Restriction archetypes by post-run hysteresis penalty.

Detects pre-tirzepatide restriction runs using only numeric variables:
calories, carbs, protein, duration, and steps.

Question:
Which kinds of restriction leave the smallest post-run expenditure penalty,
measured as the TDEE/RMR ratio in the 30 days after the run relative to the
30 days before it?

Outputs:
    analysis/restriction_runs.csv
    analysis/restriction_archetype_summary.csv
    analysis/restriction_archetype_regression.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

RESTRICTION_THRESHOLD = 1800
MIN_RUN_DAYS = 3
WINDOW_DAYS = 30

LONG_RUN_DAYS = 6
DEEP_CALORIES = 1550
LOW_CARB_G = 170
LOW_PROTEIN_G = 58
HIGH_STEPS = 4200


def load_data():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    steps = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    daily = intake.merge(steps[["date", "steps"]], on="date", how="left")
    daily = daily.merge(kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="left")
    daily = daily.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    daily = daily.merge(tirz[["date", "effective_level"]], on="date", how="left")
    daily["effective_level"] = daily["effective_level"].fillna(0)
    daily["steps"] = daily["steps"].fillna(daily["steps"].median())
    daily["tdee_rmr_ratio"] = daily["tdee"] / daily["expected_rmr"]
    return daily.sort_values("date").reset_index(drop=True)


def detect_runs(daily):
    d = daily[daily["effective_level"] == 0].copy()
    d["is_restriction"] = d["calories"] < RESTRICTION_THRESHOLD
    groups = (d["is_restriction"] != d["is_restriction"].shift()).cumsum()

    rows = []
    for _, grp in d.groupby(groups):
        if not bool(grp["is_restriction"].iloc[0]):
            continue
        if len(grp) < MIN_RUN_DAYS:
            continue

        start_idx = grp.index[0]
        end_idx = grp.index[-1]
        if start_idx - WINDOW_DAYS < 0 or end_idx + WINDOW_DAYS >= len(d):
            continue

        pre = d.loc[start_idx - WINDOW_DAYS:start_idx - 1]
        post = d.loc[end_idx + 1:end_idx + WINDOW_DAYS]

        rows.append(
            {
                "start_date": d.loc[start_idx, "date"].strftime("%Y-%m-%d"),
                "end_date": d.loc[end_idx, "date"].strftime("%Y-%m-%d"),
                "days": len(grp),
                "mean_calories": round(grp["calories"].mean(), 1),
                "mean_carbs_g": round(grp["carbs_g"].mean(), 1),
                "mean_protein_g": round(grp["protein_g"].mean(), 1),
                "mean_steps": round(grp["steps"].mean(), 1),
                "pre_ratio": round(pre["tdee_rmr_ratio"].mean(), 4),
                "run_ratio": round(grp["tdee_rmr_ratio"].mean(), 4),
                "post_ratio": round(post["tdee_rmr_ratio"].mean(), 4),
                "run_minus_pre_ratio": round(grp["tdee_rmr_ratio"].mean() - pre["tdee_rmr_ratio"].mean(), 4),
                "post_minus_pre_ratio": round(post["tdee_rmr_ratio"].mean() - pre["tdee_rmr_ratio"].mean(), 4),
                "post_minus_run_ratio": round(post["tdee_rmr_ratio"].mean() - grp["tdee_rmr_ratio"].mean(), 4),
                "pre_tdee": round(pre["tdee"].mean(), 1),
                "run_tdee": round(grp["tdee"].mean(), 1),
                "post_tdee": round(post["tdee"].mean(), 1),
                "post_minus_pre_tdee": round(post["tdee"].mean() - pre["tdee"].mean(), 1),
                "fat_delta_30d": round(post["fat_mass_lbs"].mean() - grp["fat_mass_lbs"].mean(), 2),
            }
        )

    runs = pd.DataFrame(rows)
    if runs.empty:
        return runs

    runs["long_run"] = runs["days"] >= LONG_RUN_DAYS
    runs["deep_run"] = runs["mean_calories"] < DEEP_CALORIES
    runs["low_carb"] = runs["mean_carbs_g"] < LOW_CARB_G
    runs["low_protein"] = runs["mean_protein_g"] < LOW_PROTEIN_G
    runs["high_steps"] = runs["mean_steps"] >= HIGH_STEPS
    return runs


def archetype_summary(runs):
    rows = []
    dimensions = ["long_run", "deep_run", "low_carb", "low_protein", "high_steps"]
    for dim in dimensions:
        for value, grp in runs.groupby(dim):
            rows.append(
                {
                    "dimension": dim,
                    "value": bool(value),
                    "n_runs": len(grp),
                    "mean_run_minus_pre_ratio": round(grp["run_minus_pre_ratio"].mean(), 4),
                    "mean_post_minus_pre_ratio": round(grp["post_minus_pre_ratio"].mean(), 4),
                    "mean_post_minus_pre_tdee": round(grp["post_minus_pre_tdee"].mean(), 1),
                    "mean_fat_delta_30d": round(grp["fat_delta_30d"].mean(), 2),
                }
            )
    return pd.DataFrame(rows)


def regression_summary(runs):
    features = ["long_run", "deep_run", "low_carb", "low_protein", "high_steps"]
    X = np.column_stack([np.ones(len(runs))] + [runs[f].astype(float).values for f in features])
    y_ratio = runs["post_minus_pre_ratio"].values
    y_tdee = runs["post_minus_pre_tdee"].values
    coef_ratio = np.linalg.lstsq(X, y_ratio, rcond=None)[0]
    coef_tdee = np.linalg.lstsq(X, y_tdee, rcond=None)[0]

    rows = []
    labels = ["intercept"] + features
    for label, c_ratio, c_tdee in zip(labels, coef_ratio, coef_tdee):
        rows.append(
            {
                "term": label,
                "coef_post_minus_pre_ratio": round(c_ratio, 5),
                "coef_post_minus_pre_tdee": round(c_tdee, 3),
            }
        )
    return pd.DataFrame(rows)


def print_report(runs, summary_df, reg_df):
    print("\n=== Restriction Archetypes by Hysteresis Penalty (Pre-tirzepatide) ===")
    print(
        f"Runs: {len(runs)}  mean pre/run/post TDEE-RMR ratio = "
        f"{runs['pre_ratio'].mean():.4f} / {runs['run_ratio'].mean():.4f} / {runs['post_ratio'].mean():.4f}"
    )
    print(
        f"Average post-run penalty: post - pre ratio = {runs['post_minus_pre_ratio'].mean():+.4f}  "
        f"({runs['post_minus_pre_tdee'].mean():+.1f} cal/day)"
    )

    print("\nDimension summaries:")
    for _, row in summary_df.iterrows():
        print(
            f"  {row['dimension']:>12}={str(row['value']):<5}  n={int(row['n_runs'])}  "
            f"post-pre ratio {row['mean_post_minus_pre_ratio']:+.4f}  "
            f"post-pre TDEE {row['mean_post_minus_pre_tdee']:+.1f}"
        )

    print("\nRegression:")
    for _, row in reg_df.iterrows():
        print(
            f"  {row['term']:>12}: ratio {row['coef_post_minus_pre_ratio']:+.5f}  "
            f"TDEE {row['coef_post_minus_pre_tdee']:+.3f}"
        )


def save_outputs(runs, summary_df, reg_df):
    runs.to_csv(ROOT / "analysis" / "J_restriction_runs.csv", index=False)
    summary_df.to_csv(ROOT / "analysis" / "J_restriction_archetype_summary.csv", index=False)
    reg_df.to_csv(ROOT / "analysis" / "J_restriction_archetype_regression.csv", index=False)


def main():
    runs = detect_runs(load_data())
    summary_df = archetype_summary(runs)
    reg_df = regression_summary(runs)
    save_outputs(runs, summary_df, reg_df)
    print_report(runs, summary_df, reg_df)


if __name__ == "__main__":
    main()
