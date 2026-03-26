"""Explain restriction archetype penalties via post-run rebound vs expenditure.

Takes the detected restriction runs and adds 7-day rebound metrics:
- next-7-day mean calories
- any binge in next 7 days
- any high day in next 7 days
- next-7-day max calories

The main question is whether the archetypes with the worst post-run TDEE/RMR
penalties are simply the ones with the strongest eating rebound, or whether
the expenditure penalty persists independently.

Outputs:
    analysis/restriction_rebound_summary.csv
    analysis/restriction_rebound_regression.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

BINGE_THRESHOLD = 2800
HIGH_THRESHOLD = 2400


def load_data():
    runs = pd.read_csv(ROOT / "analysis" / "J_restriction_runs.csv", parse_dates=["start_date", "end_date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    return runs, intake


def add_rebound_metrics(runs, intake):
    rows = []
    for _, row in runs.iterrows():
        post = intake[(intake["date"] > row["end_date"]) & (intake["date"] <= row["end_date"] + pd.Timedelta(days=7))]
        rows.append(
            {
                "end_date": row["end_date"],
                "next7_mean_calories": round(post["calories"].mean(), 1),
                "next7_binge_rate": round((post["calories"] > BINGE_THRESHOLD).mean(), 4),
                "next7_any_binge": int((post["calories"] > BINGE_THRESHOLD).any()),
                "next7_any_high": int((post["calories"] >= HIGH_THRESHOLD).any()),
                "next7_max_calories": round(post["calories"].max(), 1),
            }
        )
    return runs.merge(pd.DataFrame(rows), on="end_date", how="left")


def archetype_rebound_summary(runs):
    rows = []
    for dim in ["long_run", "deep_run", "low_carb", "low_protein", "high_steps"]:
        for value, grp in runs.groupby(dim):
            rows.append(
                {
                    "dimension": dim,
                    "value": bool(value),
                    "n_runs": len(grp),
                    "next7_mean_calories": round(grp["next7_mean_calories"].mean(), 1),
                    "next7_any_binge_rate": round(grp["next7_any_binge"].mean(), 4),
                    "next7_any_high_rate": round(grp["next7_any_high"].mean(), 4),
                    "next7_max_calories": round(grp["next7_max_calories"].mean(), 1),
                    "post_minus_pre_ratio": round(grp["post_minus_pre_ratio"].mean(), 4),
                }
            )
    return pd.DataFrame(rows)


def regression_summary(runs):
    """OLS: post-run penalty ~ archetype flags + rebound metrics."""
    features = [
        "long_run",
        "deep_run",
        "low_carb",
        "low_protein",
        "high_steps",
        "next7_mean_calories",
        "next7_any_binge",
    ]
    X = np.column_stack([np.ones(len(runs))] + [runs[f].astype(float).values for f in features])
    coef = np.linalg.lstsq(X, runs["post_minus_pre_ratio"].values, rcond=None)[0]
    rows = []
    for term, value in zip(["intercept"] + features, coef):
        rows.append({"term": term, "coef_post_minus_pre_ratio": round(value, 6)})
    return pd.DataFrame(rows)


def save_outputs(summary_df, reg_df):
    summary_df.to_csv(ROOT / "analysis" / "L_restriction_rebound_summary.csv", index=False)
    reg_df.to_csv(ROOT / "analysis" / "L_restriction_rebound_regression.csv", index=False)


def print_report(summary_df, reg_df):
    print("\n=== Restriction Rebound vs Expenditure Penalty ===")
    for _, row in summary_df.iterrows():
        print(
            f"  {row['dimension']:>12}={str(row['value']):<5}  n={int(row['n_runs'])}  "
            f"next7 mean={row['next7_mean_calories']:.0f}  any binge={row['next7_any_binge_rate']*100:5.1f}%  "
            f"post-pre ratio={row['post_minus_pre_ratio']:+.4f}"
        )

    print("\nRegression on post-pre TDEE/RMR penalty:")
    for _, row in reg_df.iterrows():
        print(f"  {row['term']:>18}: {row['coef_post_minus_pre_ratio']:+.6f}")


def main():
    runs, intake = load_data()
    runs = add_rebound_metrics(runs, intake)
    summary_df = archetype_rebound_summary(runs)
    reg_df = regression_summary(runs)
    save_outputs(summary_df, reg_df)
    print_report(summary_df, reg_df)


if __name__ == "__main__":
    main()
