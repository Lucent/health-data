"""Reproduce: tirzepatide blood level r=-0.50, 568 cal swing, tachyphylaxis, lean mass.

README/THEORIES claims about tirzepatide pharmacokinetics.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.linalg import lstsq

ROOT = Path(__file__).resolve().parent.parent

def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "composition" / "composition.csv", parse_dates=["date"])

    daily = intake[["date", "calories"]].merge(tirz, on="date", how="inner")

    # r=-0.50 partial
    days_since = (daily["date"] - daily["date"].min()).dt.days.values
    X = np.column_stack([days_since, np.ones(len(daily))])
    res_bl = daily["blood_level"].values - X @ lstsq(X, daily["blood_level"].values, rcond=None)[0]
    res_cal = daily["calories"].values - X @ lstsq(X, daily["calories"].values, rcond=None)[0]
    r = np.corrcoef(res_bl, res_cal)[0, 1]
    print(f"Blood level → intake (partial, controlling time): r={r:.4f}")

    # Weekly sawtooth
    print("\nWeekly cycle by day of injection:")
    for day_num in range(7):
        d = daily[daily["days_since_injection"] == day_num]
        if len(d) > 5:
            print(f"  Day {day_num}: {d['calories'].mean():.0f} cal  "
                  f"blood={d['blood_level'].mean():.1f}  n={len(d)}")

    # Tachyphylaxis
    print(f"\nTachyphylaxis: effective vs raw blood level → intake")
    r_raw = np.corrcoef(daily["blood_level"], daily["calories"])[0, 1]
    r_eff = np.corrcoef(daily["effective_level"], daily["calories"])[0, 1]
    print(f"  Raw blood level: r={r_raw:.4f}")
    print(f"  Effective (with decay): r={r_eff:.4f}")

    # Pre vs post intake
    pre = intake[(intake["date"] >= "2023-09-17") & (intake["date"] < "2024-09-17")]
    post = intake[intake["date"] >= "2024-09-17"]
    print(f"\nPre-tirz mean: {pre['calories'].mean():.0f} cal/day")
    print(f"Post-tirz mean: {post['calories'].mean():.0f} cal/day")
    print(f"Reduction: {pre['calories'].mean() - post['calories'].mean():.0f} cal/day "
          f"({(1 - post['calories'].mean()/pre['calories'].mean())*100:.1f}%)")

    # Lean mass preservation
    c = comp.sort_values("date")
    tirz_start = c[c["date"] >= "2024-09-17"]
    if len(tirz_start) >= 2:
        first = tirz_start.iloc[0]
        last = c.iloc[-1]
        fm_loss = first["fat_mass_lbs"] - last["fat_mass_lbs"]
        lean_change = last["lean_mass_lbs"] - first["lean_mass_lbs"]
        total = fm_loss - lean_change
        print(f"\nLean mass preservation (composition scans):")
        print(f"  FM: {first['fat_mass_lbs']:.1f} → {last['fat_mass_lbs']:.1f} ({-fm_loss:+.1f})")
        print(f"  Lean: {first['lean_mass_lbs']:.1f} → {last['lean_mass_lbs']:.1f} ({lean_change:+.1f})")
        print(f"  Fat % of total loss: {fm_loss/total*100:.0f}%")

    # Dose response
    print("\nDose-response:")
    for dose in sorted(daily["dose_mg"].unique()):
        d = daily[daily["dose_mg"] == dose]
        if len(d) > 7:
            print(f"  {dose:5.1f}mg: {d['calories'].mean():.0f} cal/day  n={len(d)}")

if __name__ == "__main__":
    main()
