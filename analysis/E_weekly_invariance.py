"""Reproduce: week-scale intake invariance REJECTED, ratio 1.64.

THEORIES claim: weekly variance 64% higher than independent draws.
Intake autocorrelation: r=0.40 at 1d, dies at 30d.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])

    # Weekly variance ratio
    intake["week"] = intake["date"].dt.strftime("%Y-%W")
    weekly = intake.groupby("week")["calories"].agg(["mean", "count"])
    weekly = weekly[weekly["count"] == 7]
    daily_std = intake["calories"].std()
    weekly_total_std = (weekly["mean"] * 7).std()
    expected = daily_std * np.sqrt(7)
    ratio = weekly_total_std / expected

    print("=== Week-scale intake invariance ===")
    print(f"Daily std: {daily_std:.0f}")
    print(f"Weekly total std: {weekly_total_std:.0f}")
    print(f"Expected if independent: {expected:.0f}")
    print(f"Ratio: {ratio:.2f} (1.0=independent, <1=homeostatic, >1=anti-homeostatic)")

    # Autocorrelation
    cal = intake["calories"].values
    print(f"\nIntake autocorrelation:")
    for lag in [1, 3, 7, 14, 30, 60, 90]:
        valid = ~(np.isnan(cal[:-lag]) | np.isnan(cal[lag:]))
        if valid.sum() > 100:
            r = np.corrcoef(cal[:-lag][valid], cal[lag:][valid])[0, 1]
            print(f"  Lag {lag:3d}d: r={r:.4f}")

if __name__ == "__main__":
    main()
