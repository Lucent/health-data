"""Reproduce: TDEE ~2100-2200 across 60 lb FM range, TDEE/RMR ratios.

README/THEORIES claim: body defends narrow expenditure band.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

def main():
    kf = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])

    daily = intake.merge(kf[["date", "fat_mass_lbs", "tdee"]], on="date", how="left")
    daily = daily.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    daily["tdee_rmr_ratio"] = daily["tdee"] / daily["expected_rmr"]

    print("=== TDEE by year ===")
    print(f"{'Year':>5} {'TDEE':>6} {'Intake':>7} {'Gap':>5} {'FM':>5} {'Ratio':>6}")
    for year in range(2011, 2027):
        mask = daily["date"].dt.year == year
        yr = daily[mask]
        if yr["tdee"].notna().sum() > 30:
            print(f"{year:5d} {yr['tdee'].median():6.0f} {yr['calories'].median():7.0f} "
                  f"{yr['tdee'].median()-yr['calories'].median():+5.0f} "
                  f"{yr['fat_mass_lbs'].median():5.0f} {yr['tdee_rmr_ratio'].median():6.3f}")

if __name__ == "__main__":
    main()
