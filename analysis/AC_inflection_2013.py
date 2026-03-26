"""Analyze the 2013→2014 inflection: what triggered the decade-long regain?

Fat mass bottomed at 17 lbs (Oct 2013) and rose every year for 10 years.
Was the trigger metabolic (TDEE adaptation forcing regain) or behavioral
(binge clustering that intake never recovered from)?
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
    daily = daily.merge(comp[["date", "expected_rmr", "ffm_lbs"]], on="date", how="left")
    daily["tdee_rmr_ratio"] = daily["tdee"] / daily["expected_rmr"]

    # Half-year periods around the inflection
    periods = [
        ("Loss (2012-H2)", "2012-06-01", "2012-12-31"),
        ("Loss (2013-H1)", "2013-01-01", "2013-06-30"),
        ("Bottom (2013-H2)", "2013-07-01", "2013-12-31"),
        ("Inflection (2014-H1)", "2014-01-01", "2014-06-30"),
        ("Regain (2014-H2)", "2014-07-01", "2014-12-31"),
        ("Regain (2015-H1)", "2015-01-01", "2015-06-30"),
        ("Regain (2015-H2)", "2015-07-01", "2015-12-31"),
    ]

    print("=== Half-year trajectory around inflection ===")
    print(f"{'Period':>25} {'Cal':>5} {'TDEE':>5} {'Gap':>5} {'FM Δ':>6} {'Ratio':>6} {'Binge%':>7}")
    for label, start, end in periods:
        p = daily[(daily["date"] >= start) & (daily["date"] <= end)]
        fm_change = p["fat_mass_lbs"].iloc[-1] - p["fat_mass_lbs"].iloc[0]
        binge_pct = (p["calories"] > 2800).mean() * 100
        print(f"{label:>25} {p['calories'].mean():5.0f} {p['tdee'].mean():5.0f} "
              f"{p['calories'].mean() - p['tdee'].mean():+5.0f} {fm_change:+6.1f} "
              f"{p['tdee_rmr_ratio'].mean():6.3f} {binge_pct:6.1f}%")

    # Monthly binge rate through the inflection
    window = daily[(daily["date"] >= "2013-01-01") & (daily["date"] <= "2014-12-31")]
    print(f"\n=== Monthly binge rate ===")
    for ym, g in window.groupby(window["date"].dt.to_period("M")):
        rate = (g["calories"] > 2800).mean() * 100
        fm = g["fat_mass_lbs"].mean()
        bar = "#" * int(rate / 2)
        print(f"  {str(ym):>8} FM={fm:4.0f}  binge={rate:5.1f}%  {bar}")

    # The key comparison: TDEE was recovering but intake was rising faster
    bottom = daily[(daily["date"] >= "2013-07-01") & (daily["date"] <= "2013-09-30")]
    regain = daily[(daily["date"] >= "2014-07-01") & (daily["date"] <= "2014-12-31")]
    print(f"\n=== Bottom vs regain ===")
    print(f"  Bottom (2013 Jul-Sep): cal={bottom['calories'].mean():.0f}  "
          f"tdee={bottom['tdee'].mean():.0f}  ratio={bottom['tdee_rmr_ratio'].mean():.3f}  "
          f"binge={( bottom['calories'] > 2800).mean()*100:.0f}%")
    print(f"  Regain (2014 Jul-Dec): cal={regain['calories'].mean():.0f}  "
          f"tdee={regain['tdee'].mean():.0f}  ratio={regain['tdee_rmr_ratio'].mean():.3f}  "
          f"binge={(regain['calories'] > 2800).mean()*100:.0f}%")
    print(f"  TDEE recovered: {regain['tdee'].mean() - bottom['tdee'].mean():+.0f} cal")
    print(f"  Intake increased: {regain['calories'].mean() - bottom['calories'].mean():+.0f} cal")
    print(f"  → Intake outran TDEE recovery by "
          f"{(regain['calories'].mean() - bottom['calories'].mean()) - (regain['tdee'].mean() - bottom['tdee'].mean()):.0f} cal")


if __name__ == "__main__":
    main()
