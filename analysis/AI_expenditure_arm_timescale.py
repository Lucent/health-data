#!/usr/bin/env python3
"""AI. Expenditure arm timescale.

AG shows the appetite arm (binge frequency) has a 50-day half-life EMA set point.
This script asks: does the expenditure arm (TDEE residual) have its own optimal
half-life, and is it different from 50 days?

Method: for each candidate half-life, compute EMA of fat mass, then correlate
(EMA - FM) with TDEE residual (Kalman TDEE minus composition-predicted RMR).
The half-life that maximizes |r| is the expenditure arm's adaptation timescale.
"""

import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def ema(series, half_life):
    alpha = 1 - np.exp(-np.log(2) / half_life)
    return series.ewm(alpha=alpha, min_periods=30).mean()


def main():
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    df = kalman.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)

    # TDEE residual: actual minus composition-predicted
    df["tdee_residual"] = df["tdee"] - df["expected_rmr"]
    df["tdee_ratio"] = df["tdee"] / df["expected_rmr"]

    # Drop rows without needed data
    df = df.dropna(subset=["fat_mass_lbs", "tdee", "expected_rmr"]).copy()

    # --- Pre-tirzepatide only (matches AG methodology) ---
    pre = df[df["effective_level"] == 0].copy()

    print("=" * 70)
    print("EXPENDITURE ARM TIMESCALE")
    print("=" * 70)
    print(f"\nPre-tirzepatide days: {len(pre)}")

    # --- Half-life sweep for TDEE residual ---
    print("\n--- Half-life sweep: SP distance vs TDEE residual ---")
    print(f"{'HL':>6}  {'r(resid)':>10}  {'r(ratio)':>10}  {'partial_r':>10}")

    results = []
    for hl in list(range(5, 25, 1)) + list(range(25, 105, 5)) + list(range(110, 410, 10)) + list(range(450, 1050, 50)):
        sp = ema(pre["fat_mass_lbs"], hl)
        dist = sp - pre["fat_mass_lbs"]  # positive = FM below SP

        valid = dist.notna() & pre["tdee_residual"].notna()
        if valid.sum() < 100:
            continue

        d = dist[valid].values
        resid = pre.loc[valid, "tdee_residual"].values
        ratio = pre.loc[valid, "tdee_ratio"].values
        fm = pre.loc[valid, "fat_mass_lbs"].values

        r_resid = np.corrcoef(d, resid)[0, 1]
        r_ratio = np.corrcoef(d, ratio)[0, 1]

        # Partial correlation controlling for FM
        X = np.column_stack([fm, np.ones(len(fm))])
        res_d = d - X @ np.linalg.lstsq(X, d, rcond=None)[0]
        res_r = resid - X @ np.linalg.lstsq(X, resid, rcond=None)[0]
        r_partial = np.corrcoef(res_d, res_r)[0, 1]

        results.append({
            "half_life": hl,
            "r_residual": r_resid,
            "r_ratio": r_ratio,
            "r_partial": r_partial,
            "n": int(valid.sum()),
        })

    results_df = pd.DataFrame(results)

    # Print sweep
    for _, row in results_df.iterrows():
        print(f"{row['half_life']:6.0f}  {row['r_residual']:10.4f}  {row['r_ratio']:10.4f}  {row['r_partial']:10.4f}")

    # Best half-lives
    best_resid = results_df.loc[results_df["r_residual"].abs().idxmax()]
    best_ratio = results_df.loc[results_df["r_ratio"].abs().idxmax()]
    best_partial = results_df.loc[results_df["r_partial"].abs().idxmax()]

    print(f"\nBest HL by |r_residual|: {best_resid['half_life']:.0f} days (r={best_resid['r_residual']:.4f})")
    print(f"Best HL by |r_ratio|:    {best_ratio['half_life']:.0f} days (r={best_ratio['r_ratio']:.4f})")
    print(f"Best HL by |r_partial|:  {best_partial['half_life']:.0f} days (r={best_partial['r_partial']:.4f})")

    # --- Compare with AG's appetite arm at same half-lives ---
    print("\n--- Comparison with appetite arm (binge rate) ---")
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    pre2 = pre.merge(intake[["date", "calories"]], on="date", how="left")
    pre2["binge"] = (pre2["calories"] > pre2["tdee"] + 1000).astype(float)
    pre2["binge_rate_90d"] = pre2["binge"].rolling(90, min_periods=30).mean()

    print(f"\n{'HL':>6}  {'r_appetite':>12}  {'r_expenditure':>14}  {'r_exp_partial':>14}")
    for hl in [30, 40, 50, 60, 80, 100, 150, 200, 300, 500]:
        sp = ema(pre2["fat_mass_lbs"], hl)
        dist = sp - pre2["fat_mass_lbs"]

        valid_a = dist.notna() & pre2["binge_rate_90d"].notna()
        valid_e = dist.notna() & pre2["tdee_residual"].notna()

        r_appetite = np.corrcoef(dist[valid_a], pre2.loc[valid_a, "binge_rate_90d"])[0, 1] if valid_a.sum() > 100 else np.nan
        r_expenditure = np.corrcoef(dist[valid_e], pre2.loc[valid_e, "tdee_residual"])[0, 1] if valid_e.sum() > 100 else np.nan

        # Partial r for expenditure
        if valid_e.sum() > 100:
            fm = pre2.loc[valid_e, "fat_mass_lbs"].values
            X = np.column_stack([fm, np.ones(len(fm))])
            d = dist[valid_e].values
            resid = pre2.loc[valid_e, "tdee_residual"].values
            res_d = d - X @ np.linalg.lstsq(X, d, rcond=None)[0]
            res_r = resid - X @ np.linalg.lstsq(X, resid, rcond=None)[0]
            r_exp_partial = np.corrcoef(res_d, res_r)[0, 1]
        else:
            r_exp_partial = np.nan

        print(f"{hl:6d}  {r_appetite:12.4f}  {r_expenditure:14.4f}  {r_exp_partial:14.4f}")

    # --- Effective sample size (autocorrelation correction) ---
    print("\n--- Effective sample size for expenditure arm ---")
    best_hl = int(best_partial["half_life"])
    sp = ema(pre["fat_mass_lbs"], best_hl)
    dist = sp - pre["fat_mass_lbs"]
    valid = dist.notna() & pre["tdee_residual"].notna()
    d_vals = dist[valid].values
    r_vals = pre.loc[valid, "tdee_residual"].values
    fm_vals = pre.loc[valid, "fat_mass_lbs"].values

    # Residualize
    X = np.column_stack([fm_vals, np.ones(len(fm_vals))])
    res_d = d_vals - X @ np.linalg.lstsq(X, d_vals, rcond=None)[0]
    res_r = r_vals - X @ np.linalg.lstsq(X, r_vals, rcond=None)[0]

    # Product of residuals for ACF
    product = res_d * res_r
    product_centered = product - product.mean()
    n = len(product_centered)
    acf_1 = np.corrcoef(product_centered[:-1], product_centered[1:])[0, 1]

    # Bartlett effective n
    if acf_1 > 0 and acf_1 < 1:
        n_eff_bartlett = n * (1 - acf_1) / (1 + acf_1)
    else:
        n_eff_bartlett = n

    # Block-based effective n (90-day blocks, matching AG)
    n_blocks_90 = n // 90

    print(f"N raw: {n}")
    print(f"Lag-1 ACF of residual product: {acf_1:.4f}")
    print(f"N effective (Bartlett): {n_eff_bartlett:.0f}")
    print(f"N blocks (90-day): {n_blocks_90}")
    print(f"r_partial at best HL ({best_hl}d): {best_partial['r_partial']:.4f}")

    # p-value at different effective n
    from scipy import stats
    for n_test in [n_eff_bartlett, n_blocks_90, 60]:
        if n_test > 2:
            t = best_partial["r_partial"] * np.sqrt((n_test - 2) / (1 - best_partial["r_partial"]**2))
            p = 2 * stats.t.sf(abs(t), df=n_test - 2)
            print(f"  p-value at n_eff={n_test:.0f}: {p:.4f}")

    # --- TDEE residual by phase at best HL ---
    print(f"\n--- TDEE residual by SP distance bin (HL={best_hl}d) ---")
    sp = ema(pre["fat_mass_lbs"], best_hl)
    pre_binned = pre.copy()
    pre_binned["sp_dist"] = sp - pre_binned["fat_mass_lbs"]
    pre_binned = pre_binned.dropna(subset=["sp_dist", "tdee_residual"])

    bins = [(-np.inf, -5), (-5, -2.5), (-2.5, 0), (0, 2.5), (2.5, 5), (5, 10), (10, np.inf)]
    labels = ["5+ below", "2.5-5 below", "0-2.5 below", "0-2.5 above", "2.5-5 above", "5-10 above", "10+ above"]

    print(f"{'Bin':>15}  {'n':>6}  {'mean TDEE-RMR':>14}  {'mean ratio':>11}")
    for (lo, hi), label in zip(bins, labels):
        mask = (pre_binned["sp_dist"] >= lo) & (pre_binned["sp_dist"] < hi)
        sub = pre_binned[mask]
        if len(sub) > 10:
            print(f"{label:>15}  {len(sub):6d}  {sub['tdee_residual'].mean():14.1f}  {sub['tdee_ratio'].mean():11.4f}")

    # --- Tautology check: does the signal survive a time lag? ---
    # If short-HL correlation is just Kalman mechanics (TDEE responds to recent
    # weight change), adding a lag between SP distance and TDEE residual should
    # kill it. If it's physiological, it should survive a 7-30 day lag.
    print(f"\n--- Tautology check: lagged SP distance vs TDEE residual ---")
    print(f"{'HL':>6}  {'lag0':>8}  {'lag7':>8}  {'lag14':>8}  {'lag30':>8}  {'lag60':>8}")
    for hl in [10, 15, 20, 30, 50, 80, 100, 200]:
        sp = ema(pre["fat_mass_lbs"], hl)
        dist = sp - pre["fat_mass_lbs"]
        row = [f"{hl:6d}"]
        for lag in [0, 7, 14, 30, 60]:
            d_lagged = dist.shift(lag)
            valid = d_lagged.notna() & pre["tdee_residual"].notna()
            if valid.sum() < 100:
                row.append(f"{'':>8}")
                continue
            fm = pre.loc[valid, "fat_mass_lbs"].values
            X = np.column_stack([fm, np.ones(len(fm))])
            dv = d_lagged[valid].values
            rv = pre.loc[valid, "tdee_residual"].values
            res_d = dv - X @ np.linalg.lstsq(X, dv, rcond=None)[0]
            res_r = rv - X @ np.linalg.lstsq(X, rv, rcond=None)[0]
            r = np.corrcoef(res_d, res_r)[0, 1]
            row.append(f"{r:8.4f}")
        print("  ".join(row))

    # --- Save results ---
    out_path = ROOT / "analysis" / "AI_expenditure_arm_sweep.csv"
    results_df.to_csv(out_path, index=False)
    print(f"\nArtifact: {out_path}")


if __name__ == "__main__":
    main()
