#!/usr/bin/env python3
"""AT. Baseline-temperature-only retest of the old temperature claims.

This replaces the old raw-vs-mean comparison script. Temperature inference in
the repo now uses only `baseline_mean` from `AS_temperature_baseline.csv`.
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

TREND_WINDOW_DAYS = 90
TREND_THRESHOLD_LBS = 3.0


def partial_r(x, y, z):
    if np.ndim(z) == 1:
        v = ~np.isnan(x) & ~np.isnan(y) & ~np.isnan(z)
    else:
        v = ~np.isnan(x) & ~np.isnan(y) & ~np.isnan(z).any(axis=1)
    if v.sum() < 20:
        return np.nan
    Z = z[v] if np.ndim(z) > 1 else z[v, np.newaxis]
    X = np.column_stack([Z, np.ones(v.sum())])
    rx = x[v] - X @ np.linalg.lstsq(X, x[v], rcond=None)[0]
    ry = y[v] - X @ np.linalg.lstsq(X, y[v], rcond=None)[0]
    return np.corrcoef(rx, ry)[0, 1]


def main():
    baseline = pd.read_csv(ROOT / "analysis" / "AS_temperature_baseline.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    intake_daily = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    df = baseline.merge(kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="left")
    df = df.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    df = df.merge(intake_daily[["date", "calories"]], on="date", how="left")
    df = df.merge(tirz[["date", "effective_level", "blood_level", "days_since_injection"]],
                  on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df["blood_level"] = df["blood_level"].fillna(0)
    df["on_tirz"] = (df["effective_level"] > 0).astype(int)
    df["tdee_resid"] = df["tdee"] - df["expected_rmr"]
    df["tdee_ratio"] = df["tdee"] / df["expected_rmr"]

    for win in [7, 14, 30]:
        df[f"intake_{win}d"] = df["calories"].rolling(win, min_periods=max(1, win // 2)).mean()

    df["fat_delta_90"] = df["fat_mass_lbs"] - df["fat_mass_lbs"].shift(TREND_WINDOW_DAYS)
    df["phase"] = np.where(
        df["fat_delta_90"] <= -TREND_THRESHOLD_LBS,
        "falling",
        np.where(df["fat_delta_90"] >= TREND_THRESHOLD_LBS, "rising", "stable"),
    )
    df["phase_falling"] = (df["phase"] == "falling").astype(int)
    df["phase_rising"] = (df["phase"] == "rising").astype(int)
    df = df.dropna(subset=["baseline_mean", "fat_mass_lbs", "tdee", "expected_rmr"])
    fm = df["fat_mass_lbs"].values

    print(f"Days with baseline temperature + Kalman: {len(df)}")
    print(f"  Pre-tirz: {(df['on_tirz'] == 0).sum()}, On-tirz: {(df['on_tirz'] == 1).sum()}")

    print("\nX1. TRAILING INTAKE -> BASELINE TEMPERATURE")
    print(f"  {'Window':>8} {'r':>8} {'r|FM':>8}")
    for win in [7, 14, 30]:
        col = f"intake_{win}d"
        valid = df[col].notna()
        r = np.corrcoef(df.loc[valid, col], df.loc[valid, "baseline_mean"])[0, 1] if valid.sum() > 30 else np.nan
        pr = partial_r(df["baseline_mean"].values, df[col].values, fm)
        print(f"  {win:>6}d {r:+8.3f} {pr:+8.3f}")

    print("\nX2. DRUG BLOOD LEVEL -> BASELINE TEMPERATURE")
    on = df[df["on_tirz"] == 1]
    if len(on) > 30:
        fm_on = on["fat_mass_lbs"].values
        for label, col in [("blood_level", "blood_level"), ("effective_level", "effective_level")]:
            r = np.corrcoef(on[col], on["baseline_mean"])[0, 1]
            pr = partial_r(on["baseline_mean"].values, on[col].values, fm_on)
            print(f"  {label:>20}: r={r:+.3f}  r|FM={pr:+.3f}")

    print("\nX3. PRE-TIRZ vs ON-TIRZ BASELINE TEMPERATURE")
    pre = df[df["on_tirz"] == 0]
    on = df[df["on_tirz"] == 1]
    print(f"  Pre-tirz mean={pre['baseline_mean'].mean():.3f}  n={len(pre)}")
    print(f"  On-tirz mean ={on['baseline_mean'].mean():.3f}  n={len(on)}")
    print(f"  Delta        ={on['baseline_mean'].mean() - pre['baseline_mean'].mean():+.3f}")

    print("\nX5. BASELINE TEMPERATURE -> KALMAN TDEE")
    for label, col in [("TDEE", "tdee"), ("TDEE residual", "tdee_resid"),
                       ("TDEE/RMR ratio", "tdee_ratio"), ("expected_rmr", "expected_rmr")]:
        r = np.corrcoef(df["baseline_mean"], df[col])[0, 1]
        pr = partial_r(df["baseline_mean"].values, df[col].values, fm)
        print(f"  {label:>20}: r={r:+.3f}  r|FM={pr:+.3f}")

    print("\nX6. PHASE-MATCHED BASELINE TEMPERATURE AT SAME FAT MASS")
    for lo, hi in [(60, 70), (65, 75), (70, 80), (75, 85)]:
        band = df[(df["fat_mass_lbs"] >= lo) & (df["fat_mass_lbs"] < hi)]
        for phase in ["falling", "stable", "rising"]:
            for drug_label, drug_mask in [("pre", band["on_tirz"] == 0), ("on", band["on_tirz"] == 1)]:
                sub = band[drug_mask & (band["phase"] == phase)]
                if len(sub) < 5:
                    continue
                print(f"  FM {lo}-{hi} {phase:>7} {drug_label:>3}: base={sub['baseline_mean'].mean():.3f} n={len(sub)}")

    print("\nX7. REGRESSION -- baseline only")
    X = np.column_stack(
        [
            np.ones(len(df)),
            df["fat_mass_lbs"].values,
            df["on_tirz"].values,
            df["phase_falling"].values,
            df["phase_rising"].values,
            df["effective_level"].values,
        ]
    )
    y = df["baseline_mean"].values
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    labels = ["intercept", "fat_mass_lbs", "on_tirz", "phase_falling", "phase_rising", "effective_level"]
    for label, value in zip(labels, beta):
        print(f"  {label:>14}: {value:+.6f}")

    print("\nInjection-day sawtooth:")
    on_tirz = df[df["on_tirz"] == 1].copy()
    valid = on_tirz["days_since_injection"].notna() & on_tirz["days_since_injection"].between(0, 6)
    if valid.sum() >= 30:
        for day in range(7):
            sub = on_tirz[on_tirz["days_since_injection"] == day]
            if len(sub) >= 5:
                print(f"  day {day}: base={sub['baseline_mean'].mean():+.3f} n={len(sub)}")
        r = np.corrcoef(on_tirz.loc[valid, "days_since_injection"], on_tirz.loc[valid, "baseline_mean"])[0, 1]
        print(f"  days-post-injection -> baseline r={r:+.3f}")


if __name__ == "__main__":
    main()
