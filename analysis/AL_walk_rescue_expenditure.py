#!/usr/bin/env python3
"""AL. Does walking rescue the expenditure defense that tirzepatide suppresses?

AJ shows tirz eliminates 63% of the falling-phase TDEE bonus.
AD shows walk sessions raise RMR at +14 cal/session (30d trailing).
AH#5 says exercise and set point are independent.

If walking is additive, tirz-era days with high walk frequency should show
higher TDEE/RMR than low-walk tirz days, partially offsetting the lost
expenditure defense. Test with Kalman TDEE (all days) and calorimetry
(where available).
"""

import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def partial_corr(x, y, z):
    if z.ndim == 1:
        z = z.reshape(-1, 1)
    Z = np.column_stack([z, np.ones(len(z))])
    res_x = x - Z @ np.linalg.lstsq(Z, x, rcond=None)[0]
    res_y = y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
    return np.corrcoef(res_x, res_y)[0, 1]


def main():
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    exercises = pd.read_csv(ROOT / "steps-sleep" / "exercises.csv", parse_dates=["date"])
    steps = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    rmr = pd.read_csv(ROOT / "RMR" / "rmr.csv", parse_dates=["date"])

    walks = exercises[exercises["type"] == "walking"].copy()

    # Build daily table
    df = kalman.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    df = df.merge(tirz[["date", "effective_level", "blood_level"]], on="date", how="left")
    df = df.merge(steps[["date", "steps"]], on="date", how="left")
    df = df.merge(intake[["date", "calories"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df["blood_level"] = df["blood_level"].fillna(0)
    df["on_tirz"] = df["effective_level"] > 0

    df["tdee_ratio"] = df["tdee"] / df["expected_rmr"]
    df["tdee_residual"] = df["tdee"] - df["expected_rmr"]

    # Count walk sessions in trailing windows
    df = df.sort_values("date").reset_index(drop=True)
    for window in [14, 30, 60]:
        counts = []
        for d in df["date"]:
            start = d - pd.Timedelta(days=window)
            n = len(walks[(walks["date"] > start) & (walks["date"] <= d)])
            counts.append(n)
        df[f"walks_{window}d"] = counts

    # Phase classification
    TREND_WINDOW = 90
    TREND_THRESHOLD = 3.0
    df["fat_delta"] = df["fat_mass_lbs"].diff(TREND_WINDOW)
    df["phase"] = df["fat_delta"].apply(
        lambda x: "falling" if x <= -TREND_THRESHOLD else ("rising" if x >= TREND_THRESHOLD else "stable") if pd.notna(x) else np.nan
    )

    df = df.dropna(subset=["fat_mass_lbs", "tdee", "expected_rmr"]).copy()

    print("=" * 70)
    print("DOES WALKING RESCUE THE EXPENDITURE DEFENSE ON TIRZEPATIDE?")
    print("=" * 70)

    # --- 1. Tirz era: walk frequency vs TDEE/RMR ---
    tirz_days = df[df["on_tirz"]].copy()
    pre_days = df[~df["on_tirz"]].copy()

    print(f"\nTirz days: {len(tirz_days)}")
    print(f"Pre-tirz days: {len(pre_days)}")

    print("\n--- Walk frequency distribution (tirz era) ---")
    print(f"  walks_30d: mean={tirz_days['walks_30d'].mean():.1f}, "
          f"std={tirz_days['walks_30d'].std():.1f}, "
          f"min={tirz_days['walks_30d'].min()}, max={tirz_days['walks_30d'].max()}")

    # --- 2. Tirz era: walk bins vs TDEE ---
    print("\n--- Tirz era: TDEE/RMR by walk frequency (30d) ---")
    walk_bins = [(0, 5), (5, 10), (10, 15), (15, 20), (20, 30), (30, 50)]
    print(f"{'Walks/30d':>12}  {'n':>5}  {'TDEE/RMR':>9}  {'TDEE-RMR':>9}  {'FM':>6}  {'eff_lvl':>8}  {'intake':>7}")
    for lo, hi in walk_bins:
        sub = tirz_days[(tirz_days["walks_30d"] >= lo) & (tirz_days["walks_30d"] < hi)]
        if len(sub) < 10:
            continue
        print(f"{lo:>5}-{hi:<4}  {len(sub):5d}  {sub['tdee_ratio'].mean():9.4f}  "
              f"{sub['tdee_residual'].mean():9.1f}  {sub['fat_mass_lbs'].mean():6.1f}  "
              f"{sub['effective_level'].mean():8.2f}  {sub['calories'].mean():7.0f}")

    # --- 3. Same for pre-tirz era ---
    print("\n--- Pre-tirz era: TDEE/RMR by walk frequency (30d) ---")
    print(f"{'Walks/30d':>12}  {'n':>5}  {'TDEE/RMR':>9}  {'TDEE-RMR':>9}  {'FM':>6}  {'intake':>7}")
    for lo, hi in walk_bins:
        sub = pre_days[(pre_days["walks_30d"] >= lo) & (pre_days["walks_30d"] < hi)]
        if len(sub) < 10:
            continue
        print(f"{lo:>5}-{hi:<4}  {len(sub):5d}  {sub['tdee_ratio'].mean():9.4f}  "
              f"{sub['tdee_residual'].mean():9.1f}  {sub['fat_mass_lbs'].mean():6.1f}  "
              f"{sub['calories'].mean():7.0f}")

    # --- 4. Partial correlations within tirz era ---
    print("\n--- Tirz era: correlations of walks_30d with TDEE ---")
    t = tirz_days.dropna(subset=["walks_30d", "tdee_ratio", "fat_mass_lbs", "effective_level"])
    r_raw = np.corrcoef(t["walks_30d"], t["tdee_ratio"])[0, 1]
    r_fm = partial_corr(t["walks_30d"].values, t["tdee_ratio"].values, t["fat_mass_lbs"].values)
    controls = np.column_stack([t["fat_mass_lbs"].values, t["effective_level"].values])
    r_fm_lvl = partial_corr(t["walks_30d"].values, t["tdee_ratio"].values, controls)
    print(f"  Raw: r = {r_raw:.3f}")
    print(f"  Controlling for FM: r = {r_fm:.3f}")
    print(f"  Controlling for FM + drug level: r = {r_fm_lvl:.3f}")

    # --- Same for pre-tirz ---
    print("\n--- Pre-tirz era: correlations of walks_30d with TDEE ---")
    p = pre_days.dropna(subset=["walks_30d", "tdee_ratio", "fat_mass_lbs"])
    r_raw_pre = np.corrcoef(p["walks_30d"], p["tdee_ratio"])[0, 1]
    r_fm_pre = partial_corr(p["walks_30d"].values, p["tdee_ratio"].values, p["fat_mass_lbs"].values)
    print(f"  Raw: r = {r_raw_pre:.3f}")
    print(f"  Controlling for FM: r = {r_fm_pre:.3f}")

    # --- 5. Interaction: walks × tirz on TDEE ---
    print("\n--- Regression: TDEE/RMR ~ FM + on_tirz + walks_30d + walks×tirz ---")
    valid = df.dropna(subset=["tdee_ratio", "fat_mass_lbs", "walks_30d"]).copy()
    # Only use days with Samsung walk data (post-2016)
    valid = valid[valid["date"] >= "2016-01-01"]

    X = np.column_stack([
        np.ones(len(valid)),
        valid["fat_mass_lbs"].values,
        valid["on_tirz"].astype(float).values,
        valid["walks_30d"].values,
        (valid["on_tirz"].astype(float) * valid["walks_30d"]).values,
    ])
    y = valid["tdee_ratio"].values
    coef = np.linalg.lstsq(X, y, rcond=None)[0]
    labels = ["intercept", "FM", "on_tirz", "walks_30d", "walks×tirz"]
    for name, c in zip(labels, coef):
        print(f"  {name:>15}: {c:.6f}")

    walk_effect_pre = coef[3]
    walk_effect_tirz = coef[3] + coef[4]
    print(f"\n  Walk effect (per session/30d) pre-tirz:  {walk_effect_pre:.6f} on ratio")
    print(f"  Walk effect (per session/30d) on tirz:    {walk_effect_tirz:.6f} on ratio")
    print(f"  Interaction (walks×tirz):                 {coef[4]:.6f}")

    # Convert to cal/day
    mean_rmr = valid["expected_rmr"].mean()
    print(f"\n  At mean RMR ({mean_rmr:.0f}):")
    print(f"  Walk effect pre-tirz:  {walk_effect_pre * mean_rmr:.1f} cal/day per walk-session/30d")
    print(f"  Walk effect on tirz:   {walk_effect_tirz * mean_rmr:.1f} cal/day per walk-session/30d")

    # --- 6. Phase-specific: does walking help during tirz falling? ---
    print("\n--- Tirz falling: walk frequency vs TDEE ---")
    tirz_falling = tirz_days[tirz_days["phase"] == "falling"].copy()
    print(f"Tirz falling days: {len(tirz_falling)}")

    if len(tirz_falling) > 30:
        # Median split
        median_walks = tirz_falling["walks_30d"].median()
        hi_walk = tirz_falling[tirz_falling["walks_30d"] >= median_walks]
        lo_walk = tirz_falling[tirz_falling["walks_30d"] < median_walks]
        print(f"  Median walks_30d: {median_walks:.0f}")
        print(f"  High-walk (n={len(hi_walk)}): TDEE/RMR={hi_walk['tdee_ratio'].mean():.4f}, "
              f"TDEE-RMR={hi_walk['tdee_residual'].mean():.1f}, FM={hi_walk['fat_mass_lbs'].mean():.1f}, "
              f"walks={hi_walk['walks_30d'].mean():.1f}")
        print(f"  Low-walk  (n={len(lo_walk)}): TDEE/RMR={lo_walk['tdee_ratio'].mean():.4f}, "
              f"TDEE-RMR={lo_walk['tdee_residual'].mean():.1f}, FM={lo_walk['fat_mass_lbs'].mean():.1f}, "
              f"walks={lo_walk['walks_30d'].mean():.1f}")
        print(f"  Δ TDEE/RMR: {hi_walk['tdee_ratio'].mean() - lo_walk['tdee_ratio'].mean():.4f}")
        print(f"  Δ TDEE-RMR: {hi_walk['tdee_residual'].mean() - lo_walk['tdee_residual'].mean():.1f} cal/day")

    # --- 7. Pre-tirz falling comparison ---
    print("\n--- Pre-tirz falling: walk frequency vs TDEE ---")
    pre_falling = pre_days[pre_days["phase"] == "falling"].copy()
    pre_falling = pre_falling[pre_falling["date"] >= "2016-01-01"]  # Samsung era only
    print(f"Pre-tirz falling days (Samsung era): {len(pre_falling)}")

    if len(pre_falling) > 30:
        median_walks_pre = pre_falling["walks_30d"].median()
        hi_walk_pre = pre_falling[pre_falling["walks_30d"] >= median_walks_pre]
        lo_walk_pre = pre_falling[pre_falling["walks_30d"] < median_walks_pre]
        print(f"  Median walks_30d: {median_walks_pre:.0f}")
        print(f"  High-walk (n={len(hi_walk_pre)}): TDEE/RMR={hi_walk_pre['tdee_ratio'].mean():.4f}, "
              f"TDEE-RMR={hi_walk_pre['tdee_residual'].mean():.1f}, FM={hi_walk_pre['fat_mass_lbs'].mean():.1f}, "
              f"walks={hi_walk_pre['walks_30d'].mean():.1f}")
        print(f"  Low-walk  (n={len(lo_walk_pre)}): TDEE/RMR={lo_walk_pre['tdee_ratio'].mean():.4f}, "
              f"TDEE-RMR={lo_walk_pre['tdee_residual'].mean():.1f}, FM={lo_walk_pre['fat_mass_lbs'].mean():.1f}, "
              f"walks={lo_walk_pre['walks_30d'].mean():.1f}")
        print(f"  Δ TDEE/RMR: {hi_walk_pre['tdee_ratio'].mean() - lo_walk_pre['tdee_ratio'].mean():.4f}")
        print(f"  Δ TDEE-RMR: {hi_walk_pre['tdee_residual'].mean() - lo_walk_pre['tdee_residual'].mean():.1f} cal/day")

    # --- 8. Quantify rescue: how much of AJ's lost defense do walks recover? ---
    print("\n" + "=" * 70)
    print("RESCUE QUANTIFICATION")
    print("=" * 70)

    # AJ showed: pre-tirz falling TDEE/RMR = +0.075, tirz falling = +0.012
    # Lost defense = 0.063 on ratio
    # Walk coefficient from regression: walks_30d effect on ratio
    # How many walks to recover 0.063?

    if walk_effect_tirz > 0:
        walks_to_recover = 0.063 / walk_effect_tirz
        cal_per_walk = walk_effect_tirz * mean_rmr
        cal_lost = 0.063 * mean_rmr
        print(f"\n  AJ's lost expenditure defense: 0.063 on ratio ≈ {cal_lost:.0f} cal/day")
        print(f"  Walk effect on tirz: {walk_effect_tirz:.6f} per walk/30d ≈ {cal_per_walk:.1f} cal/day")
        print(f"  Walks needed to fully recover: {walks_to_recover:.0f} sessions/30 days ({walks_to_recover/4.3:.0f}/week)")
    else:
        print(f"\n  Walk effect on tirz is zero or negative ({walk_effect_tirz:.6f}). No rescue.")

    # --- 9. Monthly timeline: walks and TDEE through the tirz era ---
    print("\n--- Monthly tirz-era timeline ---")
    tirz_monthly = tirz_days.copy()
    tirz_monthly["month"] = tirz_monthly["date"].dt.to_period("M")
    monthly = tirz_monthly.groupby("month").agg(
        n=("date", "count"),
        walks=("walks_30d", "mean"),
        tdee_ratio=("tdee_ratio", "mean"),
        tdee_residual=("tdee_residual", "mean"),
        fm=("fat_mass_lbs", "mean"),
        intake=("calories", "mean"),
        eff_level=("effective_level", "mean"),
    ).reset_index()

    print(f"{'Month':>8}  {'n':>3}  {'walks':>6}  {'ratio':>7}  {'TDEE-RMR':>9}  {'FM':>6}  {'intake':>7}  {'eff_lvl':>8}")
    for _, row in monthly.iterrows():
        print(f"{str(row['month']):>8}  {row['n']:3.0f}  {row['walks']:6.1f}  {row['tdee_ratio']:7.4f}  "
              f"{row['tdee_residual']:9.1f}  {row['fm']:6.1f}  {row['intake']:7.0f}  {row['eff_level']:8.2f}")


if __name__ == "__main__":
    main()
