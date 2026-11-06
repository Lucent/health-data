"""Reproduce: cumulative energy balance closure and TDEE/RMR by trajectory phase.

The closure residual is NOT evidence about the intake multiplier (free-drift TDEE
absorbs any level change; see DM / Kalman appendix). Phase uniformity of TDEE/RMR
is the descriptive claim checked here.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

def main():
    # Use the consumed-calorie series the canonical Kalman actually used
    # (P4's intake_true column), so TDEE and intake are in the same units.
    import sys
    sys.path.insert(0, str(ROOT / "analysis"))
    from intake_model import describe as describe_intake_model
    print(f"Intake model: {describe_intake_model()}")

    kf = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])

    intake = intake.drop(columns=["calories"]).merge(kf[["date", "intake_true"]], on="date", how="left")
    intake = intake.rename(columns={"intake_true": "calories"})

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

    print("\n=== TDEE/RMR by trajectory phase ===")
    print("  NOTE: this section previously reported an 'undercount %' computed as")
    print("  expected_rmr * 1.2 - tdee, i.e. the answer was a linear function of a")
    print("  stipulated activity factor of 1.2. It measured nothing about intake.")
    print("  The intake scale factor is estimated in CQ (regression outside the")
    print("  filter) and CR (profile innovation likelihood). What is reported here")
    print("  now is the descriptive ratio, with no assumed activity factor.")
    daily["wt_trend"] = daily["fat_mass_lbs"].rolling(90, center=True, min_periods=30).mean()
    daily["wt_slope"] = daily["wt_trend"].diff(30) / 30
    v = daily.dropna(subset=["tdee", "expected_rmr", "wt_slope"]).copy()
    v["ratio"] = v["tdee"] / v["expected_rmr"]
    for label, mask in [("Gaining", v["wt_slope"] > 0.01),
                        ("Losing", v["wt_slope"] < -0.01),
                        ("Stable", (v["wt_slope"] >= -0.01) & (v["wt_slope"] <= 0.01))]:
        group = v[mask]
        if len(group) > 30:
            print(f"  {label:8s}: TDEE/RMR={group['ratio'].median():.3f}  "
                  f"TDEE={group['tdee'].median():.0f}  "
                  f"fitted RMR={group['expected_rmr'].median():.0f}  n={len(group)}")
    print("  (uniformity across phases is the claim worth checking here — it says")
    print("   whatever the scale error is, it does not vary with weight direction)")

if __name__ == "__main__":
    main()
