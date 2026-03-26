"""Predict restriction-run metabolic failure from pre-run and run structure.

This shifts the endpoint away from binge behavior and toward direct
metabolic failure: downward movement in TDEE/RMR during and after
restriction runs.

Primary outcomes:
1. run_minus_pre_ratio: immediate branch shift during the run
2. post_minus_pre_ratio: net metabolic failure after the run
3. post_minus_run_ratio: failure to recover after the run

Outputs:
    analysis/metabolic_failure_feature_search.csv
    analysis/metabolic_failure_phase_summary.csv
    analysis/metabolic_failure_runs.csv
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
TREND_WINDOW_DAYS = 90
TREND_THRESHOLD_LBS = 3.0

OUTCOMES = [
    "run_minus_pre_ratio",
    "post_minus_pre_ratio",
    "post_minus_run_ratio",
]

FEATURES = [
    "days",
    "mean_calories",
    "mean_carbs_g",
    "mean_protein_g",
    "mean_steps",
    "pre_ratio",
    "pre_tdee",
    "protein_pct",
    "carb_pct",
    "deficit_vs_pre_tdee",
    "relative_deficit",
    "run_cal_std",
    "run_travel_frac",
    "fat_mass_start",
    "ffm_start",
    "fm_start",
    "phase_code",
]


def parse_travel() -> pd.DataFrame:
    rows = []
    for line in (ROOT / "travel" / "trips.md").read_text().splitlines():
        line = line.strip()
        if not (line.startswith("| 201") or line.startswith("| 202")):
            continue
        parts = [part.strip() for part in line.strip("|").split("|")]
        rows.append(
            {
                "start": pd.Timestamp(parts[0]),
                "end": pd.Timestamp(parts[1]),
                "location": parts[3],
            }
        )
    return pd.DataFrame(rows)


def classify_phase(delta_lbs: float) -> str | float:
    if pd.isna(delta_lbs):
        return np.nan
    if delta_lbs <= -TREND_THRESHOLD_LBS:
        return "falling"
    if delta_lbs >= TREND_THRESHOLD_LBS:
        return "rising"
    return "stable"


def load_enriched_runs() -> pd.DataFrame:
    runs = pd.read_csv(ROOT / "analysis" / "J_restriction_runs.csv", parse_dates=["start_date", "end_date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    steps = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    travel = parse_travel()

    daily = intake.merge(steps[["date", "steps"]], on="date", how="left")
    daily = daily.merge(kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="left")
    daily = daily.merge(comp[["date", "expected_rmr", "ffm_lbs", "fm_lbs"]], on="date", how="left")
    daily["travel"] = 0
    for _, row in travel.iterrows():
        mask = (daily["date"] >= row["start"]) & (daily["date"] <= row["end"])
        daily.loc[mask, "travel"] = 1

    fat_map = dict(zip(kalman["date"], kalman["fat_mass_lbs"]))

    rows = []
    for _, row in runs.iterrows():
        run = daily[(daily["date"] >= row["start_date"]) & (daily["date"] <= row["end_date"])].copy()
        start = daily[daily["date"] == row["start_date"]].copy()

        prev_date = row["start_date"] - pd.Timedelta(days=TREND_WINDOW_DAYS)
        delta = fat_map.get(row["start_date"], np.nan) - fat_map.get(prev_date, np.nan)
        phase = classify_phase(delta)
        phase_code = {"falling": -1, "stable": 0, "rising": 1}.get(phase, np.nan)

        rows.append(
            {
                **row.to_dict(),
                "protein_pct": 4 * row["mean_protein_g"] / row["mean_calories"],
                "carb_pct": 4 * row["mean_carbs_g"] / row["mean_calories"],
                "deficit_vs_pre_tdee": row["pre_tdee"] - row["mean_calories"],
                "relative_deficit": (row["pre_tdee"] - row["mean_calories"]) / row["pre_tdee"],
                "run_cal_std": run["calories"].std(),
                "run_travel_frac": run["travel"].mean(),
                "fat_mass_start": start["fat_mass_lbs"].iloc[0] if not start.empty else np.nan,
                "ffm_start": start["ffm_lbs"].iloc[0] if not start.empty else np.nan,
                "fm_start": start["fm_lbs"].iloc[0] if not start.empty else np.nan,
                "phase": phase,
                "phase_code": phase_code,
            }
        )

    return pd.DataFrame(rows)


def loocv_stats(df: pd.DataFrame, cols: list[str], target: str) -> tuple[float, float, int]:
    valid = df.dropna(subset=cols + [target]).copy()
    x = valid[cols].to_numpy(float)
    y = valid[target].to_numpy(float)
    preds = []

    for i in range(len(valid)):
        mask = np.ones(len(valid), dtype=bool)
        mask[i] = False
        train_x = x[mask]
        train_y = y[mask]
        mean = train_x.mean(axis=0)
        std = train_x.std(axis=0)
        std[std == 0] = 1.0
        train_x = (train_x - mean) / std
        test_x = (x[i] - mean) / std
        coef = np.linalg.lstsq(np.column_stack([np.ones(len(train_x)), train_x]), train_y, rcond=None)[0]
        preds.append(np.r_[1, test_x] @ coef)

    preds = np.array(preds)
    sse = ((y - preds) ** 2).sum()
    sst = ((y - y.mean()) ** 2).sum()
    r2 = 1 - sse / sst
    corr = np.corrcoef(preds, y)[0, 1]
    return r2, corr, len(valid)


def feature_search(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for target in OUTCOMES:
        single_rows = []
        for feature in FEATURES:
            r2, corr, n_rows = loocv_stats(runs, [feature], target)
            single_rows.append((r2, corr, feature, n_rows))
            rows.append(
                {
                    "target": target,
                    "model_size": 1,
                    "features": feature,
                    "loocv_r2": round(r2, 4),
                    "loocv_corr": round(corr, 4),
                    "n_rows": n_rows,
                }
            )

        top_features = [feature for _, _, feature, _ in sorted(single_rows, reverse=True)[:8]]
        for size in [2, 3, 4]:
            best = None
            for combo in combinations(top_features, size):
                r2, corr, n_rows = loocv_stats(runs, list(combo), target)
                if best is None or r2 > best[0]:
                    best = (r2, corr, combo, n_rows)
            rows.append(
                {
                    "target": target,
                    "model_size": size,
                    "features": ",".join(best[2]),
                    "loocv_r2": round(best[0], 4),
                    "loocv_corr": round(best[1], 4),
                    "n_rows": best[3],
                }
            )

    return pd.DataFrame(rows).sort_values(["target", "model_size", "loocv_r2"], ascending=[True, True, False])


def phase_summary(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (phase, high_steps), grp in runs.groupby(["phase", "high_steps"]):
        rows.append(
            {
                "phase": phase,
                "high_steps": bool(high_steps),
                "n_runs": len(grp),
                "mean_steps": round(grp["mean_steps"].mean(), 1),
                "mean_calories": round(grp["mean_calories"].mean(), 1),
                "mean_run_minus_pre_ratio": round(grp["run_minus_pre_ratio"].mean(), 4),
                "mean_post_minus_pre_ratio": round(grp["post_minus_pre_ratio"].mean(), 4),
                "mean_post_minus_run_ratio": round(grp["post_minus_run_ratio"].mean(), 4),
            }
        )
    return pd.DataFrame(rows).sort_values(["phase", "high_steps"])


def save_outputs(search_df: pd.DataFrame, phase_df: pd.DataFrame, runs: pd.DataFrame) -> None:
    search_df.to_csv(ROOT / "analysis" / "R_metabolic_failure_feature_search.csv", index=False)
    phase_df.to_csv(ROOT / "analysis" / "R_metabolic_failure_phase_summary.csv", index=False)
    runs.to_csv(ROOT / "analysis" / "R_metabolic_failure_runs.csv", index=False)


def print_report(search_df: pd.DataFrame, phase_df: pd.DataFrame) -> None:
    print("\n=== Metabolic Failure Predictors (Restriction Runs) ===")
    for target in OUTCOMES:
        best_single = search_df[(search_df["target"] == target) & (search_df["model_size"] == 1)].sort_values(
            "loocv_r2", ascending=False
        ).iloc[0]
        best_two = search_df[(search_df["target"] == target) & (search_df["model_size"] == 2)].sort_values(
            "loocv_r2", ascending=False
        ).iloc[0]
        print(
            f"{target:>22}: best single {best_single['features']} "
            f"(LOOCV R^2={best_single['loocv_r2']:.3f}) | best pair {best_two['features']} "
            f"(R^2={best_two['loocv_r2']:.3f})"
        )

    print("\nPhase x steps summary:")
    for _, row in phase_df.iterrows():
        print(
            f"  {row['phase']:>7}  high_steps={str(row['high_steps']):<5}  n={int(row['n_runs'])}  "
            f"run-pre {row['mean_run_minus_pre_ratio']:+.4f}  "
            f"post-pre {row['mean_post_minus_pre_ratio']:+.4f}"
        )


def main() -> None:
    runs = load_enriched_runs()
    search_df = feature_search(runs)
    phase_df = phase_summary(runs)
    save_outputs(search_df, phase_df, runs)
    print_report(search_df, phase_df)


if __name__ == "__main__":
    main()
