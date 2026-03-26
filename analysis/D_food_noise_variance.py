"""Reproduce: food noise manifests as variance, not mean.

THEORIES claims: CV 0.24-0.28 pre → 0.19-0.20 on tirz (25% reduction),
distance→intake -30 cal/kg, rebound r=-0.065, post-restriction binge 26%→13%.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.linalg import lstsq

ROOT = Path(__file__).resolve().parent.parent

def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kf = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    daily = intake.merge(kf[["date", "fat_mass_lbs", "tdee"]], on="date", how="left")
    daily = daily.merge(tirz[["date", "effective_level"]], on="date", how="left")
    daily["effective_level"] = daily["effective_level"].fillna(0)
    daily["on_tirz"] = (daily["effective_level"] > 0).astype(int)
    daily["sp"] = daily["fat_mass_lbs"].rolling(180, center=True, min_periods=90).mean()
    daily["dist"] = daily["fat_mass_lbs"] - daily["sp"]

    # 1. Distance → continuous intake
    pre = daily[daily["on_tirz"] == 0].dropna(subset=["dist"])
    slope = np.polyfit(pre["dist"], pre["calories"], 1)[0]
    print("=== Distance → continuous intake ===")
    print(f"Slope: {slope:.0f} cal/day per lb = {slope * 2.2046:.0f} cal/day per kg")

    # 2. Restriction duration → rebound
    daily["restricted"] = daily["calories"] < 1800
    daily["run_id"] = (daily["restricted"] != daily["restricted"].shift(1)).cumsum()
    runs = daily[daily["restricted"]].groupby("run_id").agg(
        end=("date", "last"), n_days=("date", "count")).reset_index()
    runs = runs[runs["n_days"] >= 3]
    rebounds = []
    for _, run in runs.iterrows():
        after = daily[(daily["date"] > run["end"]) &
                      (daily["date"] <= run["end"] + pd.Timedelta(days=7))]
        if len(after) >= 3:
            rebounds.append({"run_days": run["n_days"],
                             "rebound_cal": after["calories"].mean()})
    rb = pd.DataFrame(rebounds)
    r = np.corrcoef(rb["run_days"], rb["rebound_cal"])[0, 1]
    print(f"\n=== Restriction duration → rebound ===")
    print(f"r={r:.4f} (n={len(rb)} runs)")

    # 3. Intake variance by distance from set point
    print(f"\n=== Intake variance by distance ===")
    for lo, hi in [(-5, -2), (-2, 0), (0, 2), (2, 5)]:
        m = (pre["dist"] >= lo) & (pre["dist"] < hi)
        if m.sum() > 30:
            std = pre.loc[m, "calories"].std()
            cv = std / pre.loc[m, "calories"].mean()
            print(f"  {lo:+d} to {hi:+d}: std={std:.0f}  CV={cv:.3f}  n={m.sum()}")

    # 4. Tirzepatide CV reduction
    post_tirz = daily[daily["on_tirz"] == 1]
    cv_pre = pre["calories"].std() / pre["calories"].mean()
    cv_post = post_tirz["calories"].std() / post_tirz["calories"].mean()
    print(f"\n=== Tirzepatide CV reduction ===")
    print(f"Pre-tirz CV: {cv_pre:.3f}")
    print(f"On-tirz CV: {cv_post:.3f}")
    print(f"Reduction: {(1 - cv_post / cv_pre) * 100:.0f}%")

    # 5. Post-restriction binge rate
    pre_daily = daily[daily["on_tirz"] == 0].copy()
    post_daily = daily[daily["on_tirz"] == 1].copy()
    for label, subset in [("Pre-tirz", pre_daily), ("On-tirz", post_daily)]:
        subset = subset.copy()
        subset["restricted"] = subset["calories"] < 1800
        subset["run_id"] = (subset["restricted"] != subset["restricted"].shift(1)).cumsum()
        runs = subset[subset["restricted"]].groupby("run_id").agg(
            end=("date", "last"), n_days=("date", "count")).reset_index()
        runs = runs[runs["n_days"] >= 3]
        binge_count = 0
        total = 0
        for _, run in runs.iterrows():
            after = subset[(subset["date"] > run["end"]) &
                           (subset["date"] <= run["end"] + pd.Timedelta(days=7))]
            if len(after) >= 3:
                total += 1
                if (after["calories"] > 2800).any():
                    binge_count += 1
        if total > 0:
            print(f"\n{label}: {binge_count}/{total} restriction runs → binge within 7d ({binge_count/total*100:.0f}%)")

if __name__ == "__main__":
    main()
