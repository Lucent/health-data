"""Decompose restriction penalties into during-run branch shift vs post-run recovery.

Uses the previously detected restriction runs and separates two components:
1. `run_minus_pre_ratio`: how far TDEE/RMR moves during the run itself
2. `post_minus_run_ratio`: how much it recovers (or deteriorates further) after
   the run ends, relative to the run period

This distinguishes archetypes that are "bad during the cut" from archetypes
that are mainly "bad at recovery."

Outputs:
    analysis/restriction_branch_shift_summary.csv
    analysis/restriction_branch_shift_regression.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def load_runs():
    return pd.read_csv(ROOT / "analysis" / "J_restriction_runs.csv")


def summary_table(runs):
    rows = []
    for dim in ["long_run", "deep_run", "low_carb", "low_protein", "high_steps"]:
        for value, grp in runs.groupby(dim):
            rows.append(
                {
                    "dimension": dim,
                    "value": bool(value),
                    "n_runs": len(grp),
                    "mean_run_minus_pre_ratio": round(grp["run_minus_pre_ratio"].mean(), 4),
                    "mean_post_minus_pre_ratio": round(grp["post_minus_pre_ratio"].mean(), 4),
                    "mean_post_minus_run_ratio": round(grp["post_minus_run_ratio"].mean(), 4),
                }
            )
    return pd.DataFrame(rows)


def regression_table(runs):
    features = ["long_run", "deep_run", "low_carb", "low_protein", "high_steps"]
    X = np.column_stack([np.ones(len(runs))] + [runs[f].astype(float).values for f in features])
    rows = []
    for target in ["run_minus_pre_ratio", "post_minus_run_ratio"]:
        coef = np.linalg.lstsq(X, runs[target].values, rcond=None)[0]
        for term, value in zip(["intercept"] + features, coef):
            rows.append({"target": target, "term": term, "coef": round(value, 6)})
    return pd.DataFrame(rows)


def save_outputs(summary_df, reg_df):
    summary_df.to_csv(ROOT / "analysis" / "M_restriction_branch_shift_summary.csv", index=False)
    reg_df.to_csv(ROOT / "analysis" / "M_restriction_branch_shift_regression.csv", index=False)


def print_report(summary_df, reg_df):
    print("\n=== Restriction Branch Shift vs Recovery ===")
    for _, row in summary_df.iterrows():
        print(
            f"  {row['dimension']:>12}={str(row['value']):<5}  n={int(row['n_runs'])}  "
            f"run-pre {row['mean_run_minus_pre_ratio']:+.4f}  "
            f"post-pre {row['mean_post_minus_pre_ratio']:+.4f}  "
            f"post-run {row['mean_post_minus_run_ratio']:+.4f}"
        )

    print("\nRegression coefficients:")
    for _, row in reg_df.iterrows():
        print(f"  {row['target']:>20}  {row['term']:>12}: {row['coef']:+.6f}")


def main():
    runs = load_runs()
    summary_df = summary_table(runs)
    reg_df = regression_table(runs)
    save_outputs(summary_df, reg_df)
    print_report(summary_df, reg_df)


if __name__ == "__main__":
    main()
