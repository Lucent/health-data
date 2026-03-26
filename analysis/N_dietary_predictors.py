"""Reproduce: protein leverage, meal timing, fiber, gravitostat claims.

THEORIES claims: protein r=-0.34 same-day, front-loading r=-0.19/+0.48,
fiber r=-0.094 partial, gravitostat r=0.055 partial.
Also: all predictors of weight change table.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.linalg import lstsq

ROOT = Path(__file__).resolve().parent.parent

def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    foods = pd.read_csv(ROOT / "intake" / "intake_foods.csv", parse_dates=["date"])
    kf = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    steps = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])
    sleep = pd.read_csv(ROOT / "steps-sleep" / "sleep.csv", parse_dates=["date"])
    weight = pd.read_csv(ROOT / "weight" / "weight.csv", parse_dates=["date"])
    weight = weight[["date", "weight_lbs"]].dropna(subset=["weight_lbs"])

    daily = intake.merge(weight, on="date", how="left")
    daily = daily.merge(kf[["date", "fat_mass_lbs"]], on="date", how="left")
    daily = daily.merge(steps[["date", "steps"]], on="date", how="left")
    daily = daily.merge(sleep[["date", "sleep_hours"]], on="date", how="left")

    # === Protein leverage ===
    print("=== Protein leverage ===")
    daily["prot_pct"] = daily["protein_g"] * 4 / daily["calories"] * 100
    v = daily.dropna(subset=["prot_pct"])
    r = np.corrcoef(v["prot_pct"], v["calories"])[0, 1]
    print(f"Same-day protein % vs calories: r={r:.4f}")

    v["next_cal"] = v["calories"].shift(-1)
    v2 = v.dropna(subset=["next_cal"])
    X = np.column_stack([v2["calories"].values, np.ones(len(v2))])
    res_p = v2["prot_pct"].values - X @ lstsq(X, v2["prot_pct"].values, rcond=None)[0]
    res_c = v2["next_cal"].values - X @ lstsq(X, v2["next_cal"].values, rcond=None)[0]
    r_partial = np.corrcoef(res_p, res_c)[0, 1]
    print(f"Next-day partial (controlling today cal): r={r_partial:.4f}")

    # === Meal timing ===
    print("\n=== Meal timing / front-loading ===")
    meals = foods[foods["meal"] != "TOTAL"]
    morning = meals[meals["meal"].isin(["Breakfast", "Lunch"])].groupby("date")["calories"].sum()
    morning = morning.reset_index().rename(columns={"calories": "morning_cal"})
    joined = daily.merge(morning, on="date", how="left").dropna(subset=["morning_cal"])
    joined = joined[joined["calories"] > 200]
    joined["morning_pct"] = joined["morning_cal"] / joined["calories"] * 100
    r_pct = np.corrcoef(joined["morning_pct"], joined["calories"])[0, 1]
    r_abs = np.corrcoef(joined["morning_cal"], joined["calories"])[0, 1]
    print(f"Morning % vs total: r={r_pct:.4f}")
    print(f"Morning abs vs total: r={r_abs:.4f}")

    # Fiber
    fiber_meals = meals[meals["meal"].isin(["Breakfast", "Lunch"])].groupby("date")["fiber_g"].sum()
    fiber_meals = fiber_meals.reset_index().rename(columns={"fiber_g": "morning_fiber"})
    fj = joined.merge(fiber_meals, on="date", how="left").dropna(subset=["morning_fiber"])
    X = np.column_stack([fj["morning_cal"].values, np.ones(len(fj))])
    res_f = fj["morning_fiber"].values - X @ lstsq(X, fj["morning_fiber"].values, rcond=None)[0]
    res_c = fj["calories"].values - X @ lstsq(X, fj["calories"].values, rcond=None)[0]
    r_fiber = np.corrcoef(res_f, res_c)[0, 1]
    print(f"Morning fiber partial (controlling morning cal): r={r_fiber:.4f}")

    # === Gravitostat ===
    print("\n=== Gravitostat ===")
    v = daily.dropna(subset=["steps", "fat_mass_lbs"]).copy()
    v["foot_lbs"] = v["steps"] * v["fat_mass_lbs"]
    v["next_cal"] = v["calories"].shift(-1)
    v2 = v.dropna(subset=["next_cal"])
    X = np.column_stack([v2["calories"].values, np.ones(len(v2))])
    res_fp = v2["foot_lbs"].values - X @ lstsq(X, v2["foot_lbs"].values, rcond=None)[0]
    res_nc = v2["next_cal"].values - X @ lstsq(X, v2["next_cal"].values, rcond=None)[0]
    r_grav = np.corrcoef(res_fp, res_nc)[0, 1]
    print(f"Foot-pounds → next-day intake partial: r={r_grav:.4f}")

    # === All predictors of weight change ===
    print("\n=== All predictors of next-day weight change ===")
    print("(Partial r controlling calories, carbs, sodium)")
    mask = daily["weight_lbs"].notna().values
    dates = daily["date"].values[mask]
    weights = daily["weight_lbs"].values[mask]

    for col, label in [("calories", "Calories"), ("carbs_g", "Carbs"),
                        ("sodium_mg", "Sodium"), ("fat_g", "Fat"),
                        ("protein_g", "Protein"), ("fiber_g", "Fiber"),
                        ("sleep_hours", "Sleep"), ("steps", "Steps")]:
        vals, wts, cals, carbs, na = [], [], [], [], []
        col_arr = daily[col].values
        cal_arr = daily["calories"].values
        carb_arr = daily["carbs_g"].values
        na_arr = daily["sodium_mg"].values
        for i in range(1, len(dates)):
            gap = (dates[i] - dates[i - 1]) / np.timedelta64(1, "D")
            if gap != 1:
                continue
            idx = np.searchsorted(daily["date"].values, dates[i - 1])
            v = col_arr[idx]
            if np.isnan(v):
                continue
            vals.append(v)
            wts.append(weights[i] - weights[i - 1])
            cals.append(cal_arr[idx])
            carbs.append(carb_arr[idx])
            na.append(na_arr[idx])
        if len(vals) < 30:
            continue
        vals = np.array(vals)
        wts = np.array(wts)
        X = np.column_stack([np.array(cals), np.array(carbs), np.array(na), np.ones(len(vals))])
        res_v = vals - X @ lstsq(X, vals, rcond=None)[0]
        res_w = wts - X @ lstsq(X, wts, rcond=None)[0]
        r = np.corrcoef(res_v, res_w)[0, 1]
        print(f"  {label:<12s}: partial r={r:+.4f}  n={len(vals)}")

if __name__ == "__main__":
    main()
