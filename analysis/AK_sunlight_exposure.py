#!/usr/bin/env python3
"""AK. Sunlight exposure and energy balance.

AE showed sleep duration has no signal (r≈0 for everything). But sleep
duration barely varies (std 1.5h). Possible sunlight exposure — computed
from wake time + sunrise/sunset — varies 7-11h seasonally.

Two questions:
1. Does sunlight exposure predict intake, TDEE, weight change, or RMR
   better than sleep duration?
2. AD's walk sessions → RMR has a season confound (r=0.83 between walks
   and summer). Does sunlight exposure explain the walk-RMR relationship,
   or do walks remain independent after controlling for sunlight?
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent


def partial_corr(x, y, z):
    """Partial correlation of x,y controlling for z (can be matrix)."""
    if z.ndim == 1:
        z = z.reshape(-1, 1)
    Z = np.column_stack([z, np.ones(len(z))])
    res_x = x - Z @ np.linalg.lstsq(Z, x, rcond=None)[0]
    res_y = y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
    return np.corrcoef(res_x, res_y)[0, 1]


def main():
    sunlight = pd.read_csv(ROOT / "steps-sleep" / "sunlight.csv", parse_dates=["date"])
    sleep = pd.read_csv(ROOT / "steps-sleep" / "sleep.csv", parse_dates=["date"])
    sleep = sleep[sleep["sleep_hours"] >= 2]  # drop naps
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    steps = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])
    rmr = pd.read_csv(ROOT / "RMR" / "rmr.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    # Merge everything
    df = sunlight.merge(sleep[["date", "sleep_hours"]], on="date", how="left")
    df = df.merge(intake[["date", "calories", "protein_g", "carbs_g", "sodium_mg"]], on="date", how="left")
    df = df.merge(kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="left")
    df = df.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    df = df.merge(steps[["date", "steps"]], on="date", how="left")
    df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df["tdee_ratio"] = df["tdee"] / df["expected_rmr"]
    df["tdee_residual"] = df["tdee"] - df["expected_rmr"]

    # Trailing averages
    df = df.sort_values("date").reset_index(drop=True)
    for window in [7, 14, 30]:
        df[f"sunlight_{window}d"] = df["sunlight_hours"].rolling(window, min_periods=3).mean()
        df[f"sleep_{window}d"] = df["sleep_hours"].rolling(window, min_periods=3).mean()

    # Season variable
    df["month"] = pd.to_datetime(df["date"]).dt.month
    df["is_summer"] = df["month"].isin([5, 6, 7, 8]).astype(float)

    print("=" * 70)
    print("SUNLIGHT EXPOSURE AND ENERGY BALANCE")
    print("=" * 70)
    print(f"Days with sunlight + intake + Kalman: {df.dropna(subset=['sunlight_hours', 'calories', 'tdee']).shape[0]}")

    # --- 1. Sunlight vs sleep: variance comparison ---
    print("\n--- Variance comparison ---")
    valid = df.dropna(subset=["sunlight_hours", "sleep_hours"])
    print(f"Sunlight hours: mean={valid['sunlight_hours'].mean():.1f}, std={valid['sunlight_hours'].std():.1f}, CV={valid['sunlight_hours'].std()/valid['sunlight_hours'].mean():.3f}")
    print(f"Sleep hours:    mean={valid['sleep_hours'].mean():.1f}, std={valid['sleep_hours'].std():.1f}, CV={valid['sleep_hours'].std()/valid['sleep_hours'].mean():.3f}")
    print(f"Correlation (sunlight, sleep): {np.corrcoef(valid['sunlight_hours'], valid['sleep_hours'])[0,1]:.3f}")

    # --- 2. Sunlight vs energy balance variables ---
    print("\n--- Sunlight vs energy balance (raw correlations) ---")
    targets = [
        ("calories", "Same-day intake"),
        ("tdee", "Same-day TDEE"),
        ("tdee_ratio", "Same-day TDEE/RMR"),
        ("tdee_residual", "Same-day TDEE residual"),
        ("steps", "Same-day steps"),
    ]

    # Also next-day
    df["next_cal"] = df["calories"].shift(-1)
    df["next_steps"] = df["steps"].shift(-1)
    targets += [
        ("next_cal", "Next-day intake"),
        ("next_steps", "Next-day steps"),
    ]

    print(f"\n{'Target':>25}  {'r(sunlight)':>12}  {'r(sleep)':>10}  {'r(sun|FM)':>10}  {'r(slp|FM)':>10}  {'n':>6}")
    for col, label in targets:
        valid = df.dropna(subset=["sunlight_hours", "sleep_hours", col, "fat_mass_lbs"])
        if len(valid) < 50:
            continue
        sun_vals = valid["sunlight_hours"].values
        slp_vals = valid["sleep_hours"].values
        y = valid[col].values
        fm = valid["fat_mass_lbs"].values

        r_sun = np.corrcoef(sun_vals, y)[0, 1]
        r_slp = np.corrcoef(slp_vals, y)[0, 1]
        r_sun_fm = partial_corr(sun_vals, y, fm)
        r_slp_fm = partial_corr(slp_vals, y, fm)

        print(f"{label:>25}  {r_sun:12.4f}  {r_slp:10.4f}  {r_sun_fm:10.4f}  {r_slp_fm:10.4f}  {len(valid):6d}")

    # --- 3. Trailing sunlight windows ---
    print("\n--- Trailing sunlight windows vs TDEE ---")
    print(f"{'Window':>10}  {'r(sun)':>8}  {'r(sun|FM)':>10}  {'r(slp)':>8}  {'r(slp|FM)':>10}")
    for window in [7, 14, 30]:
        sun_col = f"sunlight_{window}d"
        slp_col = f"sleep_{window}d"
        valid = df.dropna(subset=[sun_col, slp_col, "tdee_ratio", "fat_mass_lbs"])
        if len(valid) < 50:
            continue
        r_sun = np.corrcoef(valid[sun_col], valid["tdee_ratio"])[0, 1]
        r_slp = np.corrcoef(valid[slp_col], valid["tdee_ratio"])[0, 1]
        r_sun_fm = partial_corr(valid[sun_col].values, valid["tdee_ratio"].values, valid["fat_mass_lbs"].values)
        r_slp_fm = partial_corr(valid[slp_col].values, valid["tdee_ratio"].values, valid["fat_mass_lbs"].values)
        print(f"{window:>8}d  {r_sun:8.4f}  {r_sun_fm:10.4f}  {r_slp:8.4f}  {r_slp_fm:10.4f}")

    # --- 4. AD's season confound: walk sessions vs RMR with sunlight control ---
    print("\n" + "=" * 70)
    print("AD SEASON CONFOUND TEST: WALK SESSIONS vs RMR")
    print("=" * 70)

    # Load exercise data to count walk sessions
    exercises = pd.read_csv(ROOT / "steps-sleep" / "exercises.csv", parse_dates=["date"])
    walks = exercises[exercises["type"] == "walking"].copy()

    # Merge RMR measurements with sunlight and walk count
    rmr_df = rmr.merge(comp[["date", "expected_rmr"]], on="date", how="inner")
    rmr_df = rmr_df.merge(kalman[["date", "fat_mass_lbs"]], on="date", how="left")
    rmr_df = rmr_df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    rmr_df["effective_level"] = rmr_df["effective_level"].fillna(0)

    # Count walk sessions in trailing 30 days
    walk_counts = []
    sunlight_30d = []
    sleep_30d_vals = []
    for _, row in rmr_df.iterrows():
        d = row["date"]
        window_start = d - pd.Timedelta(days=30)
        n_walks = len(walks[(walks["date"] > window_start) & (walks["date"] <= d)])
        walk_counts.append(n_walks)

        sun_window = sunlight[(sunlight["date"] > window_start) & (sunlight["date"] <= d)]
        sunlight_30d.append(sun_window["sunlight_hours"].mean() if len(sun_window) > 5 else np.nan)

        slp_window = sleep[(sleep["date"] > window_start) & (sleep["date"] <= d)]
        slp_window = slp_window[slp_window["sleep_hours"] >= 2]
        sleep_30d_vals.append(slp_window["sleep_hours"].mean() if len(slp_window) > 5 else np.nan)

    rmr_df["walk_sessions_30d"] = walk_counts
    rmr_df["sunlight_30d"] = sunlight_30d
    rmr_df["sleep_30d"] = sleep_30d_vals

    # Filter to Samsung era (2016+) with real data
    rmr_df = rmr_df[rmr_df["date"] >= "2016-01-01"].copy()
    rmr_df = rmr_df.dropna(subset=["sunlight_30d", "expected_rmr", "walk_sessions_30d"])

    print(f"\nRMR measurements with sunlight data: {len(rmr_df)}")

    # Correlation matrix
    print("\n--- Correlation matrix ---")
    vars_of_interest = ["rmr_kcal", "expected_rmr", "walk_sessions_30d", "sunlight_30d", "sleep_30d", "fat_mass_lbs"]
    valid_rmr = rmr_df.dropna(subset=vars_of_interest)
    print(f"n = {len(valid_rmr)}")
    mat = valid_rmr[vars_of_interest].corr()
    print(mat.round(3).to_string())

    # Key comparisons
    print("\n--- Key correlations with measured RMR ---")
    for var in ["walk_sessions_30d", "sunlight_30d", "sleep_30d", "expected_rmr", "fat_mass_lbs"]:
        v = valid_rmr.dropna(subset=[var])
        r = np.corrcoef(v[var], v["rmr_kcal"])[0, 1]
        print(f"  {var:>20}: r = {r:.3f} (n={len(v)})")

    # Partial correlations: walks controlling for sunlight, sunlight controlling for walks
    print("\n--- Partial correlations with measured RMR ---")
    v = valid_rmr.dropna(subset=["walk_sessions_30d", "sunlight_30d", "rmr_kcal", "expected_rmr"])

    r_walks_sun = partial_corr(v["walk_sessions_30d"].values, v["rmr_kcal"].values, v["sunlight_30d"].values)
    r_sun_walks = partial_corr(v["sunlight_30d"].values, v["rmr_kcal"].values, v["walk_sessions_30d"].values)

    print(f"  walks | sunlight:    r = {r_walks_sun:.3f}")
    print(f"  sunlight | walks:    r = {r_sun_walks:.3f}")

    # Control for expected_rmr too
    controls = np.column_stack([v["sunlight_30d"].values, v["expected_rmr"].values])
    r_walks_full = partial_corr(v["walk_sessions_30d"].values, v["rmr_kcal"].values, controls)
    controls2 = np.column_stack([v["walk_sessions_30d"].values, v["expected_rmr"].values])
    r_sun_full = partial_corr(v["sunlight_30d"].values, v["rmr_kcal"].values, controls2)

    print(f"  walks | sunlight + expected_rmr:    r = {r_walks_full:.3f}")
    print(f"  sunlight | walks + expected_rmr:    r = {r_sun_full:.3f}")

    # Control for expected_rmr and sleep too
    v2 = valid_rmr.dropna(subset=["walk_sessions_30d", "sunlight_30d", "sleep_30d", "rmr_kcal", "expected_rmr"])
    if len(v2) > 10:
        controls3 = np.column_stack([v2["sunlight_30d"].values, v2["expected_rmr"].values, v2["sleep_30d"].values])
        r_walks_all = partial_corr(v2["walk_sessions_30d"].values, v2["rmr_kcal"].values, controls3)
        controls4 = np.column_stack([v2["walk_sessions_30d"].values, v2["expected_rmr"].values, v2["sleep_30d"].values])
        r_sun_all = partial_corr(v2["sunlight_30d"].values, v2["rmr_kcal"].values, controls4)
        print(f"  walks | sunlight + expected_rmr + sleep:    r = {r_walks_all:.3f}")
        print(f"  sunlight | walks + expected_rmr + sleep:    r = {r_sun_all:.3f}")

    # --- Ridge CV regression: does sunlight improve on walks? ---
    print("\n--- Leave-one-out CV: RMR prediction models ---")
    from sklearn.linear_model import Ridge

    models = {
        "expected_rmr only": ["expected_rmr"],
        "expected_rmr + sunlight_30d": ["expected_rmr", "sunlight_30d"],
        "expected_rmr + sleep_30d": ["expected_rmr", "sleep_30d"],
        "expected_rmr + walk_sessions_30d": ["expected_rmr", "walk_sessions_30d"],
        "expected_rmr + walks + sunlight": ["expected_rmr", "walk_sessions_30d", "sunlight_30d"],
        "expected_rmr + walks + sunlight + sleep": ["expected_rmr", "walk_sessions_30d", "sunlight_30d", "sleep_30d"],
        "sunlight_30d only": ["sunlight_30d"],
        "walk_sessions_30d only": ["walk_sessions_30d"],
    }

    v_cv = valid_rmr.dropna(subset=["rmr_kcal", "expected_rmr", "walk_sessions_30d", "sunlight_30d", "sleep_30d"])
    n_cv = len(v_cv)
    y_cv = v_cv["rmr_kcal"].values

    print(f"\nn = {n_cv}")
    print(f"{'Model':>50}  {'CV RMSE':>8}  {'CV R²':>7}")

    for name, features in models.items():
        X = v_cv[features].values
        # Standardize
        X_mean = X.mean(axis=0)
        X_std = X.std(axis=0)
        X_std[X_std == 0] = 1
        X_norm = (X - X_mean) / X_std

        errors = []
        for i in range(n_cv):
            X_train = np.delete(X_norm, i, axis=0)
            y_train = np.delete(y_cv, i)
            X_test = X_norm[i:i+1]
            y_test = y_cv[i]

            model = Ridge(alpha=1.0)
            model.fit(X_train, y_train)
            pred = model.predict(X_test)[0]
            errors.append((pred - y_test) ** 2)

        rmse = np.sqrt(np.mean(errors))
        ss_res = np.sum(errors)
        ss_tot = np.sum((y_cv - y_cv.mean()) ** 2)
        r2 = 1 - ss_res / ss_tot
        print(f"{name:>50}  {rmse:8.0f}  {r2:7.3f}")

    # --- Show the data points ---
    print("\n--- RMR measurement details ---")
    cols = ["date", "rmr_kcal", "expected_rmr", "walk_sessions_30d", "sunlight_30d", "sleep_30d", "fat_mass_lbs"]
    avail = [c for c in cols if c in valid_rmr.columns]
    print(valid_rmr[avail].to_string(index=False))

    # --- 5. Sunlight vs intake/binge at daily level ---
    print("\n" + "=" * 70)
    print("SUNLIGHT vs INTAKE AND BINGES")
    print("=" * 70)

    daily = df.dropna(subset=["sunlight_hours", "calories", "tdee", "fat_mass_lbs"]).copy()
    daily["binge"] = (daily["calories"] > daily["tdee"] + 1000).astype(float)
    daily["surplus"] = daily["calories"] - daily["tdee"]

    # Bin by sunlight hours
    bins = [(5, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 14)]
    print(f"\n{'Sunlight bin':>14}  {'n':>5}  {'intake':>7}  {'surplus':>8}  {'binge%':>7}  {'TDEE/RMR':>9}")
    for lo, hi in bins:
        sub = daily[(daily["sunlight_hours"] >= lo) & (daily["sunlight_hours"] < hi)]
        if len(sub) < 30:
            continue
        print(f"{lo:>5}-{hi:<4}h  {len(sub):5d}  {sub['calories'].mean():7.0f}  "
              f"{sub['surplus'].mean():8.0f}  {100*sub['binge'].mean():7.1f}  "
              f"{sub['tdee_ratio'].mean():9.4f}")


if __name__ == "__main__":
    main()
