"""Test: does running suppress appetite more than walking at matched steps?

Exercise-induced anorexia hypothesis: high-intensity exercise (running)
suppresses appetite through GLP-1 and PYY release more than moderate
exercise (walking) at the same step count.

Era-matched to avoid confounding running periods (2014-2016) with
restriction periods.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def main():
    exercises = pd.read_csv(ROOT / "steps-sleep" / "exercises.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kf = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    steps = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])

    daily = intake.merge(kf[["date", "fat_mass_lbs", "tdee"]], on="date", how="left")
    daily = daily.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    daily["tdee_rmr_ratio"] = daily["tdee"] / daily["expected_rmr"]
    daily["next_cal"] = daily["calories"].shift(-1)
    daily["next_3d_cal"] = daily["calories"].rolling(3).mean().shift(-1)

    runs = exercises[exercises["type"] == "running"].groupby("date").agg(
        run_sessions=("duration_min", "count"), run_minutes=("duration_min", "sum")
    ).reset_index()
    walks = exercises[exercises["type"] == "walking"].groupby("date").agg(
        walk_sessions=("duration_min", "count"), walk_minutes=("duration_min", "sum")
    ).reset_index()

    daily = daily.merge(runs, on="date", how="left")
    daily = daily.merge(walks, on="date", how="left")
    daily = daily.merge(steps[["date", "steps"]], on="date", how="left")

    daily["ran"] = daily["run_sessions"].notna() & (daily["run_sessions"] > 0)
    daily["walked"] = daily["walk_sessions"].notna() & (daily["walk_sessions"] > 0)

    # Era-matched: 2014-2016 only
    era = daily[(daily["date"] >= "2014-01-01") & (daily["date"] <= "2016-12-31")].copy()
    era = era.dropna(subset=["next_cal", "steps"])
    era["exercise_type"] = "none"
    era.loc[era["walked"] & ~era["ran"], "exercise_type"] = "walk_only"
    era.loc[era["ran"], "exercise_type"] = "run"

    run_era = era[era["exercise_type"] == "run"]
    walk_era = era[era["exercise_type"] == "walk_only"]

    print(f"=== 2014-2016 era-matched run vs walk ===")
    print(f"Run days: {len(run_era)}, Walk-only days: {len(walk_era)}")

    matches = []
    for _, run_row in run_era.iterrows():
        diffs = (walk_era["steps"] - run_row["steps"]).abs()
        best_idx = diffs.idxmin()
        best = walk_era.loc[best_idx]
        if diffs[best_idx] < 2000:
            matches.append({
                "run_cal": run_row["calories"], "walk_cal": best["calories"],
                "run_next_cal": run_row["next_cal"], "walk_next_cal": best["next_cal"],
                "run_next_3d": run_row["next_3d_cal"], "walk_next_3d": best["next_3d_cal"],
                "run_steps": run_row["steps"], "walk_steps": best["steps"],
                "run_tdee_ratio": run_row["tdee_rmr_ratio"],
                "walk_tdee_ratio": best["tdee_rmr_ratio"],
            })

    mdf = pd.DataFrame(matches)
    print(f"Step-matched pairs (within 2000 steps): {len(mdf)}")
    print(f"  Steps: run {mdf['run_steps'].mean():.0f} vs walk {mdf['walk_steps'].mean():.0f}")
    print(f"  Same-day cal: run {mdf['run_cal'].mean():.0f} vs walk {mdf['walk_cal'].mean():.0f}  "
          f"diff: {mdf['run_cal'].mean() - mdf['walk_cal'].mean():+.0f}")
    print(f"  Next-day cal: run {mdf['run_next_cal'].mean():.0f} vs walk {mdf['walk_next_cal'].mean():.0f}  "
          f"diff: {mdf['run_next_cal'].mean() - mdf['walk_next_cal'].mean():+.0f}")
    print(f"  Next-3d cal: run {mdf['run_next_3d'].mean():.0f} vs walk {mdf['walk_next_3d'].mean():.0f}  "
          f"diff: {mdf['run_next_3d'].mean() - mdf['walk_next_3d'].mean():+.0f}")
    print(f"  TDEE/RMR: run {mdf['run_tdee_ratio'].mean():.3f} vs walk {mdf['walk_tdee_ratio'].mean():.3f}  "
          f"diff: {mdf['run_tdee_ratio'].mean() - mdf['walk_tdee_ratio'].mean():+.4f}")


if __name__ == "__main__":
    main()
