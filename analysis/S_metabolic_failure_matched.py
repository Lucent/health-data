"""Check whether restriction-run phase effects survive matching on fat mass.

Uses the enriched restriction-run table and compares falling vs rising runs
within starting-fat-mass bands, then fits simple regressions controlling for
fat mass directly.

Outputs:
    analysis/metabolic_failure_matched_bands.csv
    analysis/metabolic_failure_matched_regression.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
FAT_BANDS = [20, 35, 50, 65, 80, 95]
OUTCOMES = ["run_minus_pre_ratio", "post_minus_pre_ratio", "post_minus_run_ratio"]


def load_runs():
    runs = pd.read_csv(ROOT / "analysis" / "R_metabolic_failure_runs.csv")
    runs = runs.dropna(subset=["phase", "fat_mass_start"]).copy()
    return runs[runs["phase"].isin(["falling", "rising"])].copy()


def matched_band_summary(runs):
    rows = []
    runs = runs.copy()
    runs["fat_band"] = pd.cut(runs["fat_mass_start"], bins=FAT_BANDS, right=False)
    for band, grp in runs.groupby("fat_band", observed=False):
        if grp.empty:
            continue
        falling = grp[grp["phase"] == "falling"]
        rising = grp[grp["phase"] == "rising"]
        if len(falling) < 5 or len(rising) < 5:
            continue
        row = {
            "fat_band": str(band),
            "falling_runs": len(falling),
            "rising_runs": len(rising),
            "mean_fat_mass_start": round(grp["fat_mass_start"].mean(), 2),
        }
        for outcome in OUTCOMES:
            row[f"falling_{outcome}"] = round(falling[outcome].mean(), 4)
            row[f"rising_{outcome}"] = round(rising[outcome].mean(), 4)
            row[f"falling_minus_rising_{outcome}"] = round(
                falling[outcome].mean() - rising[outcome].mean(), 4
            )
        rows.append(row)
    return pd.DataFrame(rows)


def matched_regression(runs):
    rows = []
    d = runs.copy()
    d["falling_flag"] = (d["phase"] == "falling").astype(float)
    for outcome in OUTCOMES:
        valid = d.dropna(subset=[outcome, "fat_mass_start"]).copy()
        X = np.column_stack(
            [
                np.ones(len(valid)),
                valid["fat_mass_start"].values,
                valid["falling_flag"].values,
            ]
        )
        coef = np.linalg.lstsq(X, valid[outcome].values, rcond=None)[0]
        for term, value in zip(["intercept", "fat_mass_start", "falling_vs_rising"], coef):
            rows.append(
                {
                    "target": outcome,
                    "term": term,
                    "coef": round(value, 6),
                }
            )
    return pd.DataFrame(rows)


def save_outputs(band_df, reg_df):
    band_df.to_csv(ROOT / "analysis" / "S_metabolic_failure_matched_bands.csv", index=False)
    reg_df.to_csv(ROOT / "analysis" / "S_metabolic_failure_matched_regression.csv", index=False)


def print_report(band_df, reg_df):
    print("\n=== Metabolic Failure Matched by Start Fat Mass ===")
    for _, row in band_df.iterrows():
        print(
            f"{row['fat_band']:>10}: falling n={int(row['falling_runs'])} vs rising n={int(row['rising_runs'])} | "
            f"post-pre {row['falling_post_minus_pre_ratio']:+.4f} vs {row['rising_post_minus_pre_ratio']:+.4f}"
        )

    print("\nFat-mass-controlled regression:")
    for _, row in reg_df.iterrows():
        print(f"  {row['target']:>22}  {row['term']:>18}: {row['coef']:+.6f}")


def main():
    runs = load_runs()
    band_df = matched_band_summary(runs)
    reg_df = matched_regression(runs)
    save_outputs(band_df, reg_df)
    print_report(band_df, reg_df)


if __name__ == "__main__":
    main()
