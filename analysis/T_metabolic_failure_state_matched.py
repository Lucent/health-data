"""Match restriction runs on fat mass and pre-run TDEE/RMR state.

This checks whether falling-phase restriction runs still underperform rising
ones after matching on:
1. starting fat mass
2. pre-run TDEE/RMR ratio

Outputs:
    analysis/metabolic_failure_state_match_summary.csv
    analysis/metabolic_failure_state_match_pairs.csv
    analysis/metabolic_failure_state_match_regression.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUTCOMES = ["run_minus_pre_ratio", "post_minus_pre_ratio", "post_minus_run_ratio"]
STRICT_DISTANCE = 0.75
STRICTEST_DISTANCE = 0.50


def load_runs():
    runs = pd.read_csv(ROOT / "analysis" / "R_metabolic_failure_runs.csv")
    runs = runs.dropna(subset=["phase", "fat_mass_start", "pre_ratio"]).copy()
    return runs[runs["phase"].isin(["falling", "rising"])].copy()


def matched_regression(runs):
    rows = []
    d = runs.copy()
    d["falling_flag"] = (d["phase"] == "falling").astype(float)
    for outcome in OUTCOMES:
        X = np.column_stack(
            [
                np.ones(len(d)),
                d["fat_mass_start"].values,
                d["pre_ratio"].values,
                d["falling_flag"].values,
            ]
        )
        coef = np.linalg.lstsq(X, d[outcome].values, rcond=None)[0]
        for term, value in zip(
            ["intercept", "fat_mass_start", "pre_ratio", "falling_vs_rising"],
            coef,
        ):
            rows.append({"target": outcome, "term": term, "coef": round(value, 6)})
    return pd.DataFrame(rows)


def nearest_neighbor_pairs(runs):
    d = runs.copy()
    falling = d[d["phase"] == "falling"].copy().reset_index(drop=True)
    rising = d[d["phase"] == "rising"].copy().reset_index(drop=True)

    for col in ["fat_mass_start", "pre_ratio"]:
        mean = d[col].mean()
        std = d[col].std()
        falling[f"{col}_z"] = (falling[col] - mean) / std
        rising[f"{col}_z"] = (rising[col] - mean) / std

    rows = []
    for i, row in falling.iterrows():
        dist = np.sqrt(
            (rising["fat_mass_start_z"] - row["fat_mass_start_z"]) ** 2
            + (rising["pre_ratio_z"] - row["pre_ratio_z"]) ** 2
        )
        if dist.empty:
            continue
        j = dist.idxmin()
        pair = {
            "falling_index": i,
            "rising_index": int(j),
            "match_distance": round(float(dist.loc[j]), 4),
            "falling_fat_mass_start": round(row["fat_mass_start"], 2),
            "rising_fat_mass_start": round(rising.loc[j, "fat_mass_start"], 2),
            "falling_pre_ratio": round(row["pre_ratio"], 4),
            "rising_pre_ratio": round(rising.loc[j, "pre_ratio"], 4),
        }
        for outcome in OUTCOMES:
            pair[f"falling_{outcome}"] = round(row[outcome], 4)
            pair[f"rising_{outcome}"] = round(rising.loc[j, outcome], 4)
            pair[f"falling_minus_rising_{outcome}"] = round(row[outcome] - rising.loc[j, outcome], 4)
        rows.append(pair)
    return pd.DataFrame(rows)


def match_summary(pair_df):
    rows = []
    for label, mask in [
        ("all_pairs", pair_df["match_distance"] >= 0),
        ("distance_le_0.75", pair_df["match_distance"] <= STRICT_DISTANCE),
        ("distance_le_0.50", pair_df["match_distance"] <= STRICTEST_DISTANCE),
    ]:
        grp = pair_df[mask]
        row = {
            "subset": label,
            "n_pairs": int(len(grp)),
            "mean_match_distance": round(grp["match_distance"].mean(), 4),
        }
        for outcome in OUTCOMES:
            row[f"mean_falling_minus_rising_{outcome}"] = round(
                grp[f"falling_minus_rising_{outcome}"].mean(), 4
            )
        rows.append(row)
    return pd.DataFrame(rows)


def save_outputs(summary_df, pair_df, reg_df):
    summary_df.to_csv(ROOT / "analysis" / "T_metabolic_failure_state_match_summary.csv", index=False)
    pair_df.to_csv(ROOT / "analysis" / "T_metabolic_failure_state_match_pairs.csv", index=False)
    reg_df.to_csv(ROOT / "analysis" / "T_metabolic_failure_state_match_regression.csv", index=False)


def print_report(summary_df, reg_df):
    print("\n=== Metabolic Failure Matched on Fat Mass + Pre-run TDEE/RMR ===")
    for _, row in summary_df.iterrows():
        print(
            f"{row['subset']:>16}: n={int(row['n_pairs'])}  dist={row['mean_match_distance']:.3f}  "
            f"run-pre {row['mean_falling_minus_rising_run_minus_pre_ratio']:+.4f}  "
            f"post-pre {row['mean_falling_minus_rising_post_minus_pre_ratio']:+.4f}"
        )

    print("\nState-matched regression:")
    for _, row in reg_df.iterrows():
        print(f"  {row['target']:>22}  {row['term']:>18}: {row['coef']:+.6f}")


def main():
    runs = load_runs()
    pair_df = nearest_neighbor_pairs(runs)
    summary_df = match_summary(pair_df)
    reg_df = matched_regression(runs)
    save_outputs(summary_df, pair_df, reg_df)
    print_report(summary_df, reg_df)


if __name__ == "__main__":
    main()
