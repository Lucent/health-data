"""Reproduce: 5 negative intake-side set point tests.

THEORIES cross-cutting claims: r=+0.26 weight gain→more eating,
ratio 1.64, autocorrelation dies at 30d, binge peaks at FM 50,
binge rate 12.2%→3.4% on tirz.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kf = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    daily = intake.merge(kf[["date", "fat_mass_lbs", "tdee"]], on="date", how="left")
    daily = daily.merge(tirz[["date", "effective_level"]], on="date", how="left")
    daily["effective_level"] = daily["effective_level"].fillna(0)

    # Test 1: weight change → future intake
    daily["fm_change_30d"] = daily["fat_mass_lbs"].diff(30)
    daily["next_7d_intake"] = daily["calories"].rolling(7).mean().shift(-7)
    v = daily.dropna(subset=["fm_change_30d", "next_7d_intake"])
    r = np.corrcoef(v["fm_change_30d"], v["next_7d_intake"])[0, 1]
    print(f"1. Weight change → future intake: r={r:.4f} (positive = gaining→eat more)")

    # Test 2: weekly anti-compensation
    daily["week"] = daily["date"].dt.strftime("%Y-%W")
    weekly = daily.groupby("week")["calories"].agg(["mean", "count"])
    weekly = weekly[weekly["count"] == 7]
    ratio = (weekly["mean"] * 7).std() / (daily["calories"].std() * np.sqrt(7))
    print(f"2. Weekly variance ratio: {ratio:.2f} (>1 = anti-homeostatic)")

    # Test 3: autocorrelation
    cal = daily["calories"].values
    print(f"3. Intake autocorrelation:")
    for lag in [1, 7, 14, 30]:
        valid = ~(np.isnan(cal[:-lag]) | np.isnan(cal[lag:]))
        r = np.corrcoef(cal[:-lag][valid], cal[lag:][valid])[0, 1]
        print(f"   Lag {lag:2d}d: r={r:.4f}")

    # Test 4: binge rate by FM
    v2 = daily.dropna(subset=["fat_mass_lbs"])
    v2["binge"] = v2["calories"] > 2800
    print(f"4. Binge rate by FM:")
    for fm_bin in [20, 30, 40, 50, 60, 70, 80]:
        m = (v2["fat_mass_lbs"] >= fm_bin - 5) & (v2["fat_mass_lbs"] < fm_bin + 5)
        if m.sum() > 30:
            rate = v2.loc[m, "binge"].mean() * 100
            print(f"   FM {fm_bin}±5: {rate:.1f}%  n={m.sum()}")

    # Test 5: binge rate pre vs post tirz
    pre = daily[daily["effective_level"] == 0]
    post = daily[daily["effective_level"] > 0]
    pre_rate = (pre["calories"] > 2800).mean() * 100
    post_rate = (post["calories"] > 2800).mean() * 100
    print(f"5. Binge rate pre-tirz: {pre_rate:.1f}%  post: {post_rate:.1f}%  "
          f"({pre_rate / post_rate:.1f}× reduction)")

if __name__ == "__main__":
    main()
