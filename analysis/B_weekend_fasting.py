"""Reproduce: weekend fasting microcosm claims.

THEORIES claim: 3300 cal deficit, 0.76 lbs FM loss by Mon, back by Fri+7,
no compensatory overeating (post-fast 2369 vs pre-fast 2447).
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kf = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    daily = intake.merge(kf[["date", "fat_mass_lbs", "tdee"]], on="date", how="left")
    daily["dow"] = daily["date"].dt.dayofweek

    fasts = []
    for i in range(1, len(daily) - 7):
        if (daily.iloc[i]["dow"] == 5 and daily.iloc[i]["calories"] < 200
                and i + 1 < len(daily) and daily.iloc[i + 1]["dow"] == 6
                and daily.iloc[i + 1]["calories"] < 200):
            pre = daily.iloc[max(0, i - 5):i]
            post = daily.iloc[i + 2:i + 7]
            if len(pre) >= 3 and len(post) >= 3:
                fm_fri = daily.iloc[i - 1]["fat_mass_lbs"]
                fm_mon = daily.iloc[i + 2]["fat_mass_lbs"]
                fm_wk = daily.iloc[i + 7]["fat_mass_lbs"] if i + 7 < len(daily) else np.nan
                tdee_sat = daily.iloc[i]["tdee"]
                tdee_sun = daily.iloc[i + 1]["tdee"]
                deficit = (tdee_sat - daily.iloc[i]["calories"]) + (tdee_sun - daily.iloc[i + 1]["calories"])
                fasts.append({
                    "date": daily.iloc[i]["date"].strftime("%Y-%m-%d"),
                    "pre_mean": pre["calories"].mean(),
                    "post_mean": post["calories"].mean(),
                    "deficit": deficit,
                    "fm_weekend": fm_mon - fm_fri,
                    "fm_week": fm_wk - fm_fri if not np.isnan(fm_wk) else np.nan,
                })

    df = pd.DataFrame(fasts)
    print(f"Weekend 36-hour fasts: {len(df)}")
    print(f"Pre-fast baseline: {df['pre_mean'].mean():.0f} cal/day")
    print(f"Post-fast week: {df['post_mean'].mean():.0f} cal/day")
    print(f"Mean deficit per fast: {df['deficit'].mean():.0f} cal")
    print(f"Expected fat loss: {df['deficit'].mean()/3500:.2f} lbs")
    print(f"Actual FM change Fri→Mon: {df['fm_weekend'].mean():.2f} lbs")
    print(f"Actual FM change Fri→Fri+7: {df['fm_week'].mean():.2f} lbs")
    print(f"Rebound excess: {(df['post_mean'].mean() - df['pre_mean'].mean()) * 5:.0f} cal (negative = no rebound)")
    print(f"% deficit recovered: {(df['post_mean'].mean() - df['pre_mean'].mean()) * 5 / df['deficit'].mean() * 100:.0f}%")

if __name__ == "__main__":
    main()
