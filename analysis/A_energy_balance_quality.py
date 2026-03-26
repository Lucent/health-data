"""Reproduce: ±5 lbs cumulative energy balance, ~10-15% undercount uniform.

README claim: cumulative energy balance closes to ±5 lbs over 15 years,
0.08 lbs/year drift. Undercount ~10-15%, uniform across phases.
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

    valid = daily.dropna(subset=["tdee", "fat_mass_lbs"]).copy()
    valid["surplus"] = valid["calories"] - valid["tdee"]
    valid["cum_surplus"] = valid["surplus"].cumsum() / 3500
    valid["actual_change"] = valid["fat_mass_lbs"] - valid["fat_mass_lbs"].iloc[0]
    residual = valid["cum_surplus"] - valid["actual_change"]

    print("=== Cumulative energy balance ===")
    print(f"Max abs residual: {residual.abs().max():.1f} lbs")
    print(f"Per-year drift: {residual.iloc[-1] / (len(valid)/365.25):.2f} lbs/year")

    print("\n=== Undercount by trajectory phase ===")
    daily["wt_trend"] = daily["fat_mass_lbs"].rolling(90, center=True, min_periods=30).mean()
    daily["wt_slope"] = daily["wt_trend"].diff(30) / 30
    v = daily.dropna(subset=["tdee", "expected_rmr", "wt_slope"]).copy()
    v["ratio"] = v["tdee"] / v["expected_rmr"]
    for label, mask in [("Gaining", v["wt_slope"] > 0.01),
                        ("Losing", v["wt_slope"] < -0.01),
                        ("Stable", (v["wt_slope"] >= -0.01) & (v["wt_slope"] <= 0.01))]:
        group = v[mask]
        if len(group) > 30:
            undercount = group["expected_rmr"].median() * 1.2 - group["tdee"].median()
            pct = undercount / group["calories"].median() * 100
            print(f"  {label:8s}: TDEE/RMR={group['ratio'].median():.3f}  "
                  f"undercount={undercount:.0f} cal ({pct:.1f}%)  n={len(group)}")

if __name__ == "__main__":
    main()
