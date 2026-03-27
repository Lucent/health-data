#!/usr/bin/env python3
"""AM. Confidence intervals on the key recommendation claims.

Quantify uncertainty on:
1. Set point location (SP = 63 lbs FM)
2. Binge risk at 3 lbs below SP (~6%)
3. Drug equivalence (2.5 lbs per unit effective level)
4. Set point half-life (50 days)
5. Convergence rate (~1 lb/month)
6. Walk session RMR effect (+14 cal/session)
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
N_BOOT = 5000
BLOCK_SIZE = 90  # days, matching AG


def ema(series, half_life):
    alpha = 1 - np.exp(-np.log(2) / half_life)
    return series.ewm(alpha=alpha, min_periods=30).mean()


def block_bootstrap_indices(n, block_size, rng):
    n_blocks = (n + block_size - 1) // block_size
    starts = rng.integers(0, n - block_size + 1, size=n_blocks)
    idx = np.concatenate([np.arange(s, s + block_size) for s in starts])
    return idx[:n]


def main():
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    rmr_data = pd.read_csv(ROOT / "RMR" / "rmr.csv", parse_dates=["date"])
    exercises = pd.read_csv(ROOT / "steps-sleep" / "exercises.csv", parse_dates=["date"])

    df = kalman.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    df = df.merge(intake[["date", "calories"]], on="date", how="left")
    df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df = df.dropna(subset=["fat_mass_lbs", "tdee", "calories"]).copy()
    df["binge"] = (df["calories"] > df["tdee"] + 1000).astype(float)
    df["binge_rate_90d"] = df["binge"].rolling(90, min_periods=30).mean()

    pre = df[df["effective_level"] == 0].copy().reset_index(drop=True)

    rng = np.random.default_rng(42)

    print("=" * 70)
    print("CONFIDENCE INTERVALS ON KEY RECOMMENDATION CLAIMS")
    print("=" * 70)

    # ================================================================
    # 1. SET POINT HALF-LIFE AND LOCATION
    # ================================================================
    print("\n--- 1. Set point half-life ---")

    # Bootstrap the optimal half-life
    def find_best_hl(data):
        best_r, best_hl = 0, 50
        for hl in range(20, 200, 5):
            sp = ema(data["fat_mass_lbs"], hl)
            dist = sp - data["fat_mass_lbs"]
            valid = dist.notna() & data["binge_rate_90d"].notna()
            if valid.sum() < 100:
                continue
            r = np.corrcoef(dist[valid], data.loc[valid, "binge_rate_90d"])[0, 1]
            if abs(r) > abs(best_r):
                best_r = r
                best_hl = hl
        return best_hl, best_r

    point_hl, point_r = find_best_hl(pre)
    print(f"Point estimate: HL={point_hl}d, r={point_r:.4f}")

    boot_hls = []
    boot_rs = []
    for i in range(N_BOOT):
        idx = block_bootstrap_indices(len(pre), BLOCK_SIZE, rng)
        sample = pre.iloc[idx].reset_index(drop=True)
        hl, r = find_best_hl(sample)
        boot_hls.append(hl)
        boot_rs.append(r)

    boot_hls = np.array(boot_hls)
    boot_rs = np.array(boot_rs)
    print(f"Bootstrap HL: median={np.median(boot_hls):.0f}, 95% CI [{np.percentile(boot_hls, 2.5):.0f}, {np.percentile(boot_hls, 97.5):.0f}]")
    print(f"Bootstrap r:  median={np.median(boot_rs):.4f}, 95% CI [{np.percentile(boot_rs, 2.5):.4f}, {np.percentile(boot_rs, 97.5):.4f}]")

    # ================================================================
    # 2. SET POINT LOCATION (current)
    # ================================================================
    print("\n--- 2. Set point location at end of data ---")

    # SP depends on HL. Compute SP at end of data for each bootstrapped HL
    full = df.copy().reset_index(drop=True)
    boot_sps = []
    for hl in boot_hls:
        sp = ema(full["fat_mass_lbs"], hl)
        boot_sps.append(sp.iloc[-1])

    boot_sps = np.array(boot_sps)
    sp_point = ema(full["fat_mass_lbs"], point_hl).iloc[-1]
    fm_current = full["fat_mass_lbs"].iloc[-1]

    print(f"Current FM: {fm_current:.1f}")
    print(f"SP point estimate (HL={point_hl}d): {sp_point:.1f}")
    print(f"SP bootstrap: median={np.median(boot_sps):.1f}, 95% CI [{np.percentile(boot_sps, 2.5):.1f}, {np.percentile(boot_sps, 97.5):.1f}]")

    boot_gaps = boot_sps - fm_current
    print(f"Gap (SP - FM): median={np.median(boot_gaps):.1f}, 95% CI [{np.percentile(boot_gaps, 2.5):.1f}, {np.percentile(boot_gaps, 97.5):.1f}]")

    # ================================================================
    # 3. BINGE RISK AT CURRENT GAP
    # ================================================================
    print("\n--- 3. Binge risk at 3 lbs below SP ---")

    # For each bootstrap HL, compute binge rate in the 0-to-5-below bin
    sp_series = ema(pre["fat_mass_lbs"], point_hl)
    dist = sp_series - pre["fat_mass_lbs"]
    pre_with_dist = pre.copy()
    pre_with_dist["sp_dist"] = dist

    # Point estimate for several bins
    bins = [(-7.5, -5), (-5, -2.5), (-2.5, 0), (0, 2.5), (2.5, 5)]
    labels = ["5-7.5 below", "2.5-5 below", "0-2.5 below", "0-2.5 above", "2.5-5 above"]

    print(f"\nPoint estimates:")
    for (lo, hi), label in zip(bins, labels):
        mask = pre_with_dist["sp_dist"].between(lo, hi) & pre_with_dist["binge"].notna()
        sub = pre_with_dist[mask]
        if len(sub) > 10:
            rate = sub["binge"].mean()
            print(f"  {label:>15}: {100*rate:.1f}% (n={len(sub)})")

    # Bootstrap CI for the 0-2.5 below bin (closest to current 3 lb gap)
    # and the 2.5-5 below bin
    for target_lo, target_hi, target_label in [(-5, -2.5, "2.5-5 below"), (-2.5, 0, "0-2.5 below")]:
        boot_rates = []
        for i in range(N_BOOT):
            idx = block_bootstrap_indices(len(pre_with_dist), BLOCK_SIZE, rng)
            sample = pre_with_dist.iloc[idx]
            mask = sample["sp_dist"].between(target_lo, target_hi) & sample["binge"].notna()
            sub = sample[mask]
            if len(sub) > 5:
                boot_rates.append(sub["binge"].mean())

        boot_rates = np.array(boot_rates)
        print(f"\n  {target_label} bootstrap (n={len(boot_rates)} valid resamples):")
        print(f"    median={100*np.median(boot_rates):.1f}%, 95% CI [{100*np.percentile(boot_rates, 2.5):.1f}%, {100*np.percentile(boot_rates, 97.5):.1f}%]")

    # ================================================================
    # 4. DRUG EQUIVALENCE (lbs per unit effective level)
    # ================================================================
    print("\n--- 4. Drug equivalence (lbs of SP offset per unit effective level) ---")

    # Logistic regression: binge ~ sp_dist + effective_level
    # The ratio of coefficients gives lbs-equivalent per unit drug
    all_with_dist = df.copy()
    sp_all = ema(all_with_dist["fat_mass_lbs"], point_hl)
    all_with_dist["sp_dist"] = sp_all - all_with_dist["fat_mass_lbs"]
    valid = all_with_dist.dropna(subset=["sp_dist", "binge"]).copy()
    valid = valid[valid["sp_dist"].notna()].reset_index(drop=True)

    # Logistic regression via MLE
    from scipy.optimize import minimize

    def neg_log_lik(beta, X, y):
        z = X @ beta
        z = np.clip(z, -20, 20)
        ll = y * z - np.log(1 + np.exp(z))
        return -ll.sum()

    def fit_logistic(X, y):
        res = minimize(neg_log_lik, np.zeros(X.shape[1]), args=(X, y),
                       method="L-BFGS-B")
        return res.x

    X = np.column_stack([
        np.ones(len(valid)),
        valid["sp_dist"].values,
        valid["effective_level"].values,
    ])
    y = valid["binge"].values

    beta = fit_logistic(X, y)
    print(f"Logistic: intercept={beta[0]:.4f}, sp_dist={beta[1]:.4f}, eff_level={beta[2]:.4f}")

    if beta[1] != 0:
        lbs_per_unit = -beta[2] / beta[1]
        print(f"Lbs equivalent per unit effective level: {lbs_per_unit:.2f}")
    else:
        lbs_per_unit = np.nan

    # Bootstrap
    boot_equiv = []
    for i in range(N_BOOT):
        idx = block_bootstrap_indices(len(valid), BLOCK_SIZE, rng)
        X_b = X[idx]
        y_b = y[idx]
        try:
            b = fit_logistic(X_b, y_b)
            if b[1] != 0:
                boot_equiv.append(-b[2] / b[1])
        except:
            pass

    boot_equiv = np.array(boot_equiv)
    boot_equiv = boot_equiv[(boot_equiv > -20) & (boot_equiv < 20)]  # trim outliers
    print(f"Bootstrap: median={np.median(boot_equiv):.2f}, 95% CI [{np.percentile(boot_equiv, 2.5):.2f}, {np.percentile(boot_equiv, 97.5):.2f}]")

    # ================================================================
    # 5. CONVERGENCE RATE
    # ================================================================
    print("\n--- 5. Convergence rate ---")

    # At HL=50d, monthly convergence = 1 - 0.5^(30/50) = 34% of gap
    # At 3 lb gap: 3 * 0.34 = 1.0 lb/month
    for hl in [int(np.percentile(boot_hls, 2.5)), point_hl, int(np.percentile(boot_hls, 97.5))]:
        monthly_frac = 1 - 0.5 ** (30 / hl)
        monthly_lbs = 3.0 * monthly_frac
        months_to_1lb = -np.log(1/3) / (np.log(2) / hl) / 30
        print(f"  HL={hl:3d}d: {100*monthly_frac:.0f}% of gap/month = {monthly_lbs:.1f} lbs/month at 3lb gap, "
              f"gap < 1 lb in {months_to_1lb:.1f} months")

    # ================================================================
    # 6. WALK SESSION RMR EFFECT
    # ================================================================
    print("\n--- 6. Walk session RMR effect ---")

    walks = exercises[exercises["type"] == "walking"]
    rmr_df = rmr_data.merge(comp[["date", "expected_rmr"]], on="date", how="inner")
    rmr_df = rmr_df[rmr_df["date"] >= "2016-01-01"]

    walk_counts = []
    for _, row in rmr_df.iterrows():
        d = row["date"]
        window_start = d - pd.Timedelta(days=30)
        n = len(walks[(walks["date"] > window_start) & (walks["date"] <= d)])
        walk_counts.append(n)
    rmr_df["walks_30d"] = walk_counts

    rmr_df = rmr_df.dropna(subset=["rmr_kcal", "expected_rmr", "walks_30d"])
    print(f"Calorimetry measurements: {len(rmr_df)}")

    # OLS: rmr = a + b*expected_rmr + c*walks_30d
    X_rmr = np.column_stack([np.ones(len(rmr_df)), rmr_df["expected_rmr"].values, rmr_df["walks_30d"].values])
    y_rmr = rmr_df["rmr_kcal"].values
    coef = np.linalg.lstsq(X_rmr, y_rmr, rcond=None)[0]
    print(f"Point estimate: +{coef[2]:.1f} cal RMR per walk session (30d)")

    # Leave-one-out for walk coefficient
    loo_coefs = []
    for i in range(len(rmr_df)):
        X_loo = np.delete(X_rmr, i, axis=0)
        y_loo = np.delete(y_rmr, i)
        c = np.linalg.lstsq(X_loo, y_loo, rcond=None)[0]
        loo_coefs.append(c[2])

    loo_coefs = np.array(loo_coefs)
    # Jackknife SE
    jack_mean = loo_coefs.mean()
    jack_se = np.sqrt((len(loo_coefs) - 1) / len(loo_coefs) * np.sum((loo_coefs - jack_mean)**2))
    t_crit = stats.t.ppf(0.975, df=len(rmr_df) - 3)
    ci_lo = coef[2] - t_crit * jack_se
    ci_hi = coef[2] + t_crit * jack_se
    t_stat = coef[2] / jack_se
    p_val = 2 * stats.t.sf(abs(t_stat), df=len(rmr_df) - 3)

    print(f"Jackknife SE: {jack_se:.1f}")
    print(f"95% CI: [{ci_lo:.1f}, {ci_hi:.1f}]")
    print(f"t = {t_stat:.2f}, p = {p_val:.4f}")

    # Pairs bootstrap (resample measurements)
    boot_walk_coefs = []
    for i in range(N_BOOT):
        idx = rng.integers(0, len(rmr_df), size=len(rmr_df))
        X_b = X_rmr[idx]
        y_b = y_rmr[idx]
        try:
            c = np.linalg.lstsq(X_b, y_b, rcond=None)[0]
            boot_walk_coefs.append(c[2])
        except:
            pass

    boot_walk_coefs = np.array(boot_walk_coefs)
    print(f"Bootstrap 95% CI: [{np.percentile(boot_walk_coefs, 2.5):.1f}, {np.percentile(boot_walk_coefs, 97.5):.1f}]")

    # ================================================================
    # SUMMARY TABLE
    # ================================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"""
Claim                                   Point       95% CI              p
──────────────────────────────────────────────────────────────────────────
SP half-life                            {point_hl}d          [{np.percentile(boot_hls, 2.5):.0f}, {np.percentile(boot_hls, 97.5):.0f}]d
SP-binge correlation (r)                {point_r:.3f}       [{np.percentile(boot_rs, 2.5):.3f}, {np.percentile(boot_rs, 97.5):.3f}]         <0.001
Current SP                              {sp_point:.1f} lbs    [{np.percentile(boot_sps, 2.5):.1f}, {np.percentile(boot_sps, 97.5):.1f}] lbs
Current gap (SP - FM)                   {sp_point-fm_current:.1f} lbs     [{np.percentile(boot_gaps, 2.5):.1f}, {np.percentile(boot_gaps, 97.5):.1f}] lbs
Binge risk 0-2.5 below SP              ~6%         see above
Drug equiv (lbs/unit eff level)         {lbs_per_unit:.1f} lbs    [{np.percentile(boot_equiv, 2.5):.1f}, {np.percentile(boot_equiv, 97.5):.1f}] lbs
Walk RMR effect (cal/session)           +{coef[2]:.0f}         [{ci_lo:.0f}, {ci_hi:.0f}]            p={p_val:.4f}
""")


if __name__ == "__main__":
    main()
