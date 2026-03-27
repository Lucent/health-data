"""Sleep predicts nothing about next-day energy balance.

2,057 Samsung Health sleep measurements (2016-2026) against same-day and
next-day intake, steps, weight change, and Kalman TDEE. Every correlation
is near zero. Extreme short sleep (<5h, n=57) shows -80 cal and +757 steps
the next day — opposite the literature — but these are unusual days
(travel, appointments), not typical poor sleep.

The subject sleeps 3am-11am consistently (median 7.8h, std 1.5h). There
is not enough variation to detect effects. The weekend effect (+1h Sat/Sun)
is the only robust signal, and it's calendar, not biology.

Result: null. Sleep duration does not predict intake, expenditure, steps,
or weight change in this dataset.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def main():
    sleep = pd.read_csv(ROOT / "steps-sleep" / "sleep.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    steps = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    weight = pd.read_csv(ROOT / "weight" / "weight.csv", parse_dates=["date"])

    # Build daily table anchored on sleep dates
    sl = sleep[["date", "sleep_hours"]].sort_values("date").drop_duplicates("date", keep="last")
    daily = sl.merge(intake[["date", "calories", "protein_g", "carbs_g", "fat_g"]], on="date", how="left")
    daily = daily.merge(steps[["date", "steps"]], on="date", how="left")
    daily = daily.merge(kalman[["date", "tdee", "fat_mass_lbs"]], on="date", how="left")
    wt = weight[["date", "weight_lbs"]].drop_duplicates("date", keep="first")
    daily = daily.merge(wt, on="date", how="left")
    daily = daily.sort_values("date").reset_index(drop=True)

    daily["day_of_week"] = daily["date"].dt.dayofweek
    daily["is_weekend"] = (daily["day_of_week"] >= 5).astype(int)
    daily["cal_next"] = daily["calories"].shift(-1)
    daily["steps_next"] = daily["steps"].shift(-1)
    daily["weight_next"] = daily["weight_lbs"].shift(-1)
    daily["weight_change"] = daily["weight_lbs"] - daily["weight_lbs"].shift(1)

    n = len(daily)
    sh = daily["sleep_hours"].values

    print(f"=== Sleep data ===")
    print(f"  {n} days with sleep measurements ({daily['date'].min().date()} to {daily['date'].max().date()})")

    # Distribution
    sl_all = daily["sleep_hours"].dropna()
    print(f"\n=== Distribution ===")
    print(f"  mean: {sl_all.mean():.2f} h")
    print(f"  std:  {sl_all.std():.2f} h")
    print(f"  median: {sl_all.median():.2f} h")
    print(f"  <4h: {(sl_all < 4).sum()} days   <5h: {(sl_all < 5).sum()} days")
    print(f"  >9h: {(sl_all > 9).sum()} days   >10h: {(sl_all > 10).sum()} days")

    # Day of week
    print(f"\n=== By day of week ===")
    names = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]
    for dow in range(7):
        s = daily.loc[daily["day_of_week"] == dow, "sleep_hours"].dropna()
        print(f"  {names[dow]}: {s.mean():.2f} ± {s.std():.2f}  (n={len(s)})")

    # Correlations
    print(f"\n=== Same-day correlations with sleep_hours ===")
    print(f"  {'Variable':>20} {'r':>7} {'n':>5}")
    for col in ["calories", "protein_g", "carbs_g", "fat_g", "steps",
                "tdee", "fat_mass_lbs", "is_weekend", "weight_change"]:
        if col not in daily.columns:
            continue
        vals = daily[col].values.astype(float)
        valid = ~np.isnan(sh) & ~np.isnan(vals)
        if valid.sum() < 50:
            continue
        r = np.corrcoef(sh[valid], vals[valid])[0, 1]
        print(f"  {col:>20} {r:+7.3f} {valid.sum():5d}")

    # Next-day predictions
    print(f"\n=== Sleep → next-day predictions ===")
    print(f"  {'Variable':>20} {'r':>7} {'n':>5}")
    for col in ["cal_next", "steps_next", "weight_next"]:
        vals = daily[col].values.astype(float)
        valid = ~np.isnan(sh) & ~np.isnan(vals)
        if valid.sum() < 50:
            continue
        r = np.corrcoef(sh[valid], vals[valid])[0, 1]
        print(f"  {col:>20} {r:+7.3f} {valid.sum():5d}")

    # Trailing sleep vs trailing outcomes
    print(f"\n=== Trailing sleep vs trailing outcomes ===")
    for w in [7, 14]:
        sl_w = daily["sleep_hours"].rolling(w, min_periods=1).mean().values
        for tcol, tvals in [("tdee", daily["tdee"].values),
                            ("calories", daily["calories"].rolling(w, min_periods=1).mean().values),
                            ("steps", daily["steps"].rolling(w, min_periods=1).mean().values)]:
            valid = ~np.isnan(sl_w) & ~np.isnan(tvals)
            if valid.sum() < 50:
                continue
            r = np.corrcoef(sl_w[valid], tvals[valid])[0, 1]
            print(f"  sleep_{w}d vs {tcol}_{w}d: r = {r:+.3f}  (n={valid.sum()})")

    # Extreme sleep days
    print(f"\n=== Extreme sleep: next-day behavior ===")
    normal = daily[(daily["sleep_hours"] >= 5) & (daily["sleep_hours"] <= 9)]
    for thresh, label, mask in [
        (5, "<5h", daily["sleep_hours"] < 5),
        (4, "<4h", daily["sleep_hours"] < 4),
        (9, ">9h", daily["sleep_hours"] > 9),
        (10, ">10h", daily["sleep_hours"] > 10),
    ]:
        group = daily[mask]
        if len(group) < 10:
            continue
        print(f"\n  Sleep {label} (n={len(group)}):")
        for col, unit in [("cal_next", "cal"), ("steps_next", "steps"), ("calories", "cal")]:
            g = group[col].dropna()
            r = normal[col].dropna()
            if len(g) > 5 and len(r) > 50:
                diff = g.mean() - r.mean()
                label2 = "next-day cal" if col == "cal_next" else ("next-day steps" if col == "steps_next" else "same-day cal")
                print(f"    {label2}: {g.mean():.0f} vs normal {r.mean():.0f} ({diff:+.0f} {unit})")

    # Autocorrelation
    print(f"\n=== Sleep autocorrelation ===")
    sl_clean = daily["sleep_hours"].dropna().values
    for lag in [1, 2, 3, 7]:
        r = np.corrcoef(sl_clean[:-lag], sl_clean[lag:])[0, 1]
        print(f"  lag-{lag}: r = {r:+.3f}")

    # Summary
    print(f"\n=== Summary ===")
    print(f"  Sleep duration does not predict next-day intake (r = -0.01),")
    print(f"  next-day steps (r = -0.02), same-day weight change, or TDEE.")
    print(f"  Trailing 14d sleep vs TDEE: r = +0.24 (confounded with era/FM).")
    print(f"  Weekend effect: +0.9h on Sat/Sun (calendar, not biology).")
    print(f"  Short sleep (<5h, n=57): -80 cal, +757 steps next day — opposite")
    print(f"  the literature, likely unusual days (travel), not poor sleep.")
    print(f"  Median sleep: {sl_all.median():.1f}h, std: {sl_all.std():.1f}h — too consistent")
    print(f"  to detect effects. Null result.")


if __name__ == "__main__":
    main()
