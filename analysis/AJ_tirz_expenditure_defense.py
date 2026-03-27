#!/usr/bin/env python3
"""AJ. Tirzepatide and the expenditure defense direction.

F shows 206 cal metabolic clawback on tirzepatide (RMR 1930→1750).
K shows falling phases have elevated TDEE (+202 cal vs stable).
But tirzepatide is a falling phase — shouldn't TDEE be elevated, not suppressed?

Three hypotheses:
  (a) The drug suppresses the expenditure defense too
  (b) The clawback is the net after partial K elevation
  (c) K's elevation requires behavioral restriction, not pharmacological

Test: compare TDEE/RMR during tirz-era falling vs pre-tirz falling at matched FM.
"""

import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

TREND_WINDOW_DAYS = 90
TREND_THRESHOLD_LBS = 3.0


def classify_phase(delta_lbs):
    if pd.isna(delta_lbs):
        return np.nan
    if delta_lbs <= -TREND_THRESHOLD_LBS:
        return "falling"
    if delta_lbs >= TREND_THRESHOLD_LBS:
        return "rising"
    return "stable"


def main():
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    df = kalman.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    df = df.merge(tirz[["date", "effective_level", "blood_level"]], on="date", how="left")
    df = df.merge(intake[["date", "calories"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df["blood_level"] = df["blood_level"].fillna(0)
    df["on_tirz"] = df["effective_level"] > 0

    # Phase classification (same as K)
    df["fat_delta"] = df["fat_mass_lbs"].diff(TREND_WINDOW_DAYS)
    df["phase"] = df["fat_delta"].map(classify_phase)

    df["tdee_ratio"] = df["tdee"] / df["expected_rmr"]
    df["tdee_residual"] = df["tdee"] - df["expected_rmr"]

    df = df.dropna(subset=["fat_mass_lbs", "tdee", "expected_rmr", "phase"]).copy()

    print("=" * 70)
    print("TIRZEPATIDE AND THE EXPENDITURE DEFENSE")
    print("=" * 70)

    # --- Phase distribution by era ---
    print("\n--- Phase distribution ---")
    for era, label in [(~df["on_tirz"], "Pre-tirz"), (df["on_tirz"], "On tirz")]:
        sub = df[era]
        print(f"\n{label} (n={len(sub)}):")
        for phase in ["falling", "stable", "rising"]:
            n = (sub["phase"] == phase).sum()
            print(f"  {phase:>8}: {n:5d} ({100*n/len(sub):.1f}%)")

    # --- TDEE/RMR by phase × era ---
    print("\n--- TDEE/RMR ratio by phase × era ---")
    print(f"{'Era':>12}  {'Phase':>8}  {'n':>5}  {'TDEE/RMR':>9}  {'TDEE-RMR':>9}  {'FM':>6}  {'RMR':>6}")
    for era_mask, era_label in [(~df["on_tirz"], "Pre-tirz"), (df["on_tirz"], "On tirz")]:
        for phase in ["falling", "stable", "rising"]:
            sub = df[era_mask & (df["phase"] == phase)]
            if len(sub) < 10:
                continue
            print(f"{era_label:>12}  {phase:>8}  {len(sub):5d}  {sub['tdee_ratio'].mean():9.4f}  "
                  f"{sub['tdee_residual'].mean():9.1f}  {sub['fat_mass_lbs'].mean():6.1f}  "
                  f"{sub['expected_rmr'].mean():6.0f}")

    # --- FM-matched comparison: tirz falling vs pre-tirz falling ---
    print("\n--- FM-matched falling: tirz vs pre-tirz ---")
    tirz_falling = df[df["on_tirz"] & (df["phase"] == "falling")].copy()
    pre_falling = df[~df["on_tirz"] & (df["phase"] == "falling")].copy()

    print(f"Tirz falling: n={len(tirz_falling)}, FM range {tirz_falling['fat_mass_lbs'].min():.1f}-{tirz_falling['fat_mass_lbs'].max():.1f}")
    print(f"Pre-tirz falling: n={len(pre_falling)}, FM range {pre_falling['fat_mass_lbs'].min():.1f}-{pre_falling['fat_mass_lbs'].max():.1f}")

    # Find overlapping FM range
    fm_lo = max(tirz_falling["fat_mass_lbs"].quantile(0.05), pre_falling["fat_mass_lbs"].quantile(0.05))
    fm_hi = min(tirz_falling["fat_mass_lbs"].quantile(0.95), pre_falling["fat_mass_lbs"].quantile(0.95))
    print(f"Overlap FM band: {fm_lo:.1f}-{fm_hi:.1f}")

    tirz_in = tirz_falling[(tirz_falling["fat_mass_lbs"] >= fm_lo) & (tirz_falling["fat_mass_lbs"] <= fm_hi)]
    pre_in = pre_falling[(pre_falling["fat_mass_lbs"] >= fm_lo) & (pre_falling["fat_mass_lbs"] <= fm_hi)]

    if len(tirz_in) > 5 and len(pre_in) > 5:
        print(f"\nIn overlap band:")
        print(f"  Tirz falling: n={len(tirz_in)}, mean FM={tirz_in['fat_mass_lbs'].mean():.1f}, "
              f"TDEE/RMR={tirz_in['tdee_ratio'].mean():.4f}, TDEE-RMR={tirz_in['tdee_residual'].mean():.1f}")
        print(f"  Pre  falling: n={len(pre_in)}, mean FM={pre_in['fat_mass_lbs'].mean():.1f}, "
              f"TDEE/RMR={pre_in['tdee_ratio'].mean():.4f}, TDEE-RMR={pre_in['tdee_residual'].mean():.1f}")
        print(f"  Δ TDEE/RMR: {tirz_in['tdee_ratio'].mean() - pre_in['tdee_ratio'].mean():.4f}")
        print(f"  Δ TDEE-RMR: {tirz_in['tdee_residual'].mean() - pre_in['tdee_residual'].mean():.1f}")

    # --- FM-binned comparison ---
    print("\n--- FM-binned falling comparison ---")
    fm_bands = [(40, 55), (55, 70), (65, 85)]
    print(f"{'FM band':>10}  {'n_tirz':>7}  {'n_pre':>6}  {'ratio_tirz':>11}  {'ratio_pre':>10}  {'Δ':>8}")
    for lo, hi in fm_bands:
        t = tirz_falling[(tirz_falling["fat_mass_lbs"] >= lo) & (tirz_falling["fat_mass_lbs"] < hi)]
        p = pre_falling[(pre_falling["fat_mass_lbs"] >= lo) & (pre_falling["fat_mass_lbs"] < hi)]
        if len(t) >= 10 and len(p) >= 10:
            delta = t["tdee_ratio"].mean() - p["tdee_ratio"].mean()
            print(f"{lo}-{hi:>3}  {len(t):7d}  {len(p):6d}  {t['tdee_ratio'].mean():11.4f}  "
                  f"{p['tdee_ratio'].mean():10.4f}  {delta:8.4f}")

    # --- Rate of FM change: is tirz falling faster or slower? ---
    print("\n--- Rate of fat loss during falling phases ---")
    for era_mask, label in [(~df["on_tirz"], "Pre-tirz"), (df["on_tirz"], "On tirz")]:
        falling = df[era_mask & (df["phase"] == "falling")].copy()
        if len(falling) < 30:
            continue
        # lbs/month = fat_delta / TREND_WINDOW_DAYS * 30
        rate = falling["fat_delta"].mean() / TREND_WINDOW_DAYS * 30
        cal_deficit = (falling["calories"] - falling["tdee"]).mean()
        print(f"  {label}: {rate:.2f} lbs/month, mean cal deficit {cal_deficit:.0f} cal/day, "
              f"mean intake {falling['calories'].mean():.0f}")

    # --- Regression: TDEE/RMR ~ FM + phase + on_tirz + phase×on_tirz ---
    print("\n--- Regression: TDEE ratio ~ FM + phase + tirz + phase×tirz ---")
    falling_days = df[df["phase"] == "falling"].copy()
    if len(falling_days[falling_days["on_tirz"]]) >= 10:
        X = np.column_stack([
            np.ones(len(falling_days)),
            falling_days["fat_mass_lbs"].values,
            falling_days["on_tirz"].astype(float).values,
        ])
        y = falling_days["tdee_ratio"].values
        coef = np.linalg.lstsq(X, y, rcond=None)[0]
        print(f"  Falling days only (n={len(falling_days)}):")
        print(f"  intercept={coef[0]:.4f}, FM={coef[1]:.6f}, on_tirz={coef[2]:.4f}")
        print(f"  Tirz coefficient: {coef[2]:.4f} ({'suppressed' if coef[2] < 0 else 'elevated'})")

    # --- Full regression with all phases ---
    valid = df.dropna(subset=["tdee_ratio", "fat_mass_lbs"]).copy()
    X = np.column_stack([
        np.ones(len(valid)),
        valid["fat_mass_lbs"].values,
        (valid["phase"] == "falling").astype(float).values,
        (valid["phase"] == "rising").astype(float).values,
        valid["on_tirz"].astype(float).values,
        ((valid["phase"] == "falling") & valid["on_tirz"]).astype(float).values,
    ])
    y = valid["tdee_ratio"].values
    coef = np.linalg.lstsq(X, y, rcond=None)[0]
    labels = ["intercept", "FM", "falling", "rising", "on_tirz", "falling×tirz"]
    print(f"\n  Full model (n={len(valid)}):")
    for name, c in zip(labels, coef):
        print(f"    {name:>16}: {c:.6f}")

    tirz_falling_effect = coef[2] + coef[4] + coef[5]  # falling + tirz + interaction
    pre_falling_effect = coef[2]  # falling only
    print(f"\n  Pre-tirz falling effect on ratio: {pre_falling_effect:+.4f}")
    print(f"  Tirz falling total effect:        {tirz_falling_effect:+.4f}")
    print(f"  Drug suppression of falling bonus: {coef[5]:+.4f}")

    # --- Calorimetry validation ---
    print("\n--- Direct calorimetry check ---")
    rmr = pd.read_csv(ROOT / "RMR" / "rmr.csv", parse_dates=["date"])
    rmr_merged = rmr.merge(df[["date", "fat_mass_lbs", "expected_rmr", "on_tirz", "phase", "tdee", "tdee_ratio"]],
                           on="date", how="inner")
    print(f"{'Date':>12}  {'RMR_meas':>9}  {'RMR_pred':>9}  {'Δ':>6}  {'TDEE':>6}  {'ratio':>6}  {'tirz':>5}  {'phase':>8}")
    for _, row in rmr_merged.iterrows():
        delta = row["rmr_kcal"] - row["expected_rmr"]
        print(f"{row['date'].strftime('%Y-%m-%d'):>12}  {row['rmr_kcal']:9.0f}  {row['expected_rmr']:9.0f}  "
              f"{delta:6.0f}  {row['tdee']:6.0f}  {row['tdee_ratio']:6.3f}  "
              f"{'yes' if row['on_tirz'] else 'no':>5}  {row['phase']:>8}")

    # Summarize pre-tirz falling calorimetry vs tirz falling
    pre_cal = rmr_merged[(~rmr_merged["on_tirz"]) & (rmr_merged["phase"] == "falling")]
    tirz_cal = rmr_merged[rmr_merged["on_tirz"]]
    if len(pre_cal) > 0:
        print(f"\nPre-tirz falling calorimetry (n={len(pre_cal)}): "
              f"mean RMR={pre_cal['rmr_kcal'].mean():.0f}, predicted={pre_cal['expected_rmr'].mean():.0f}, "
              f"Δ={pre_cal['rmr_kcal'].mean()-pre_cal['expected_rmr'].mean():.0f}")
    if len(tirz_cal) > 0:
        print(f"On-tirz calorimetry (n={len(tirz_cal)}): "
              f"mean RMR={tirz_cal['rmr_kcal'].mean():.0f}, predicted={tirz_cal['expected_rmr'].mean():.0f}, "
              f"Δ={tirz_cal['rmr_kcal'].mean()-tirz_cal['expected_rmr'].mean():.0f}")


if __name__ == "__main__":
    main()
