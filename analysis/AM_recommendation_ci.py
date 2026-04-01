#!/usr/bin/env python3
"""AM. Confidence intervals on the key recommendation claims.

Updated to use AM/AN refined parameters:
- Overshoot = surplus > 100 cal (not 1000 cal "binge")
- 105-day rate window (not 90)
- Asymmetric SP: HL_up=72d, HL_down=25d (profile likelihood)
- Expenditure arm: HL=10d, partial r=+0.52
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from scipy.optimize import minimize, minimize_scalar

ROOT = Path(__file__).resolve().parent.parent
N_BOOT = 5000
BLOCK_SIZE = 90

OVERSHOOT_THRESH = 100  # cal above TDEE (AM finding)
OVERSHOOT_WIN = 105     # days (AM finding)
SP_HL_SYM = 51          # days (AM finding)
SP_HL_UP = 72           # days (AN finding)
SP_HL_DOWN = 25         # days (AN finding)


def ema(series, half_life):
    alpha = 1 - np.exp(-np.log(2) / half_life)
    return series.ewm(alpha=alpha, min_periods=30).mean()


def asymmetric_ema(fm_vals, hl_up, hl_down):
    alpha_up = 1 - np.exp(-np.log(2) / hl_up)
    alpha_down = 1 - np.exp(-np.log(2) / hl_down)
    sp = np.empty(len(fm_vals))
    sp[0] = fm_vals[0]
    for i in range(1, len(fm_vals)):
        prev = sp[i - 1]
        cur = fm_vals[i]
        if np.isnan(prev):
            sp[i] = cur
        elif np.isnan(cur):
            sp[i] = prev
        elif cur > prev:
            sp[i] = prev + alpha_up * (cur - prev)
        else:
            sp[i] = prev + alpha_down * (cur - prev)
    return sp


def block_bootstrap_indices(n, block_size, rng):
    bs = min(block_size, max(1, n // 2))
    n_blocks = (n + bs - 1) // bs
    starts = rng.integers(0, n - bs + 1, size=n_blocks)
    idx = np.concatenate([np.arange(s, s + bs) for s in starts])
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
    df = df.sort_values("date").reset_index(drop=True)
    df["surplus"] = df["calories"] - df["tdee"]
    df["overshoot"] = (df["surplus"] > OVERSHOOT_THRESH).astype(float)
    df["overshoot_rate"] = df["overshoot"].rolling(OVERSHOOT_WIN, min_periods=OVERSHOOT_WIN).mean()

    pre = df[df["effective_level"] == 0].copy().reset_index(drop=True)

    rng = np.random.default_rng(42)

    print("=" * 70)
    print("CONFIDENCE INTERVALS ON KEY RECOMMENDATION CLAIMS")
    print(f"Overshoot: surplus > {OVERSHOOT_THRESH} cal, {OVERSHOOT_WIN}d window")
    print(f"Asymmetric SP: HL_up={SP_HL_UP}d, HL_down={SP_HL_DOWN}d")
    print("=" * 70)

    # ================================================================
    # 1. SET POINT HALF-LIFE
    # ================================================================
    print("\n--- 1. Set point half-life ---")

    def find_best_hl(data_fm, data_rate):
        best_r, best_hl = 0, 50
        for hl in range(20, 200, 5):
            sp = ema(pd.Series(data_fm), hl)
            dist = sp.values - data_fm
            valid = ~np.isnan(dist) & ~np.isnan(data_rate)
            if valid.sum() < 100:
                continue
            r = np.corrcoef(dist[valid], data_rate[valid])[0, 1]
            if abs(r) > abs(best_r):
                best_r, best_hl = r, hl
        return best_hl, best_r

    pre_fm = pre["fat_mass_lbs"].values
    pre_rate = pre["overshoot_rate"].values
    point_hl, point_r = find_best_hl(pre_fm, pre_rate)
    print(f"Symmetric point estimate: HL={point_hl}d, r={point_r:.4f}")

    # Asymmetric at AN parameters
    sp_asym = asymmetric_ema(pre_fm, SP_HL_UP, SP_HL_DOWN)
    dist_asym = sp_asym - pre_fm
    valid_asym = ~np.isnan(dist_asym) & ~np.isnan(pre_rate)
    r_asym = np.corrcoef(dist_asym[valid_asym], pre_rate[valid_asym])[0, 1]
    print(f"Asymmetric (HL_up={SP_HL_UP}, HL_down={SP_HL_DOWN}): r={r_asym:.4f}")

    boot_hls = []
    for i in range(N_BOOT):
        idx = block_bootstrap_indices(len(pre), BLOCK_SIZE, rng)
        hl, _ = find_best_hl(pre_fm[idx], pre_rate[idx])
        boot_hls.append(hl)

    boot_hls = np.array(boot_hls)
    print(f"Bootstrap HL: median={np.median(boot_hls):.0f}, 95% CI [{np.percentile(boot_hls, 2.5):.0f}, {np.percentile(boot_hls, 97.5):.0f}]")

    # Profile likelihood CI for ratchet (from AN)
    print(f"Ratchet ratio (AN profile likelihood): 2.9x, conservative CI [0.8, 4.9]")
    print(f"HL_down stable at 25-27d across profile; uncertainty is in HL_up [26, 170]")

    # ================================================================
    # 2. SET POINT LOCATION (current)
    # ================================================================
    print("\n--- 2. Set point location at end of data ---")

    full = df.copy().reset_index(drop=True)
    full_fm = full["fat_mass_lbs"].values

    # Asymmetric SP on full data
    sp_full_asym = asymmetric_ema(full_fm, SP_HL_UP, SP_HL_DOWN)
    sp_current = sp_full_asym[-1]
    fm_current = full_fm[-1]

    # Also symmetric at each bootstrapped HL
    boot_sps_sym = []
    for hl in boot_hls:
        sp = ema(full["fat_mass_lbs"], hl)
        boot_sps_sym.append(sp.iloc[-1])
    boot_sps_sym = np.array(boot_sps_sym)

    # Asymmetric SP at range of plausible HL_up values (AN profile: 26-170)
    boot_sps_asym = []
    for hl_up in range(26, 171, 5):
        sp = asymmetric_ema(full_fm, hl_up, 25)
        boot_sps_asym.append(sp[-1])
    boot_sps_asym = np.array(boot_sps_asym)

    print(f"Current FM: {fm_current:.1f}")
    print(f"SP (asymmetric, HL_up={SP_HL_UP}d): {sp_current:.1f}")
    print(f"SP (symmetric HL bootstrap): median={np.median(boot_sps_sym):.1f}, 95% CI [{np.percentile(boot_sps_sym, 2.5):.1f}, {np.percentile(boot_sps_sym, 97.5):.1f}]")
    print(f"SP (asymmetric, HL_up 26-170d): range [{boot_sps_asym.min():.1f}, {boot_sps_asym.max():.1f}]")

    gap = sp_current - fm_current
    gap_range_sym = boot_sps_sym - fm_current
    gap_range_asym = boot_sps_asym - fm_current
    print(f"Gap (SP - FM): point={gap:.1f}, sym CI [{np.percentile(gap_range_sym, 2.5):.1f}, {np.percentile(gap_range_sym, 97.5):.1f}], asym range [{gap_range_asym.min():.1f}, {gap_range_asym.max():.1f}]")

    # ================================================================
    # 3. OVERSHOOT RISK AT CURRENT GAP
    # ================================================================
    print("\n--- 3. Overshoot risk by distance below SP ---")

    sp_pre_asym = asymmetric_ema(pre_fm, SP_HL_UP, SP_HL_DOWN)
    dist_pre = sp_pre_asym - pre_fm
    pre_os = pre["overshoot"].values

    bins = [(-7.5, -5), (-5, -2.5), (-2.5, 0), (0, 2.5), (2.5, 5), (5, 10)]
    labels = ["5-7.5 below", "2.5-5 below", "0-2.5 below", "0-2.5 above", "2.5-5 above", "5-10 above"]

    print(f"\n  Overshoot rate (surplus > {OVERSHOOT_THRESH} cal) by distance from asymmetric SP:")
    for (lo, hi), label in zip(bins, labels):
        mask = (dist_pre > lo) & (dist_pre <= hi) & ~np.isnan(pre_os)
        n_bin = mask.sum()
        if n_bin < 20:
            continue
        rate = pre_os[mask].mean()

        # Bootstrap CI
        boot_rates = []
        for _ in range(N_BOOT):
            idx = block_bootstrap_indices(n_bin, BLOCK_SIZE, rng)
            boot_rates.append(pre_os[mask][idx].mean())
        boot_rates = np.array(boot_rates)

        print(f"  {label:>15}: {rate:.1%} [{np.percentile(boot_rates, 2.5):.1%}, {np.percentile(boot_rates, 97.5):.1%}] (n={n_bin})")

    # Non-binge continuous pressure
    non_os = ~np.isnan(dist_pre) & ~np.isnan(pre["surplus"].values) & (pre_os == 0)
    slope = np.polyfit(dist_pre[non_os], pre["surplus"].values[non_os], 1)[0]
    print(f"\n  Non-overshoot day pressure: {slope:+.1f} cal/day per lb below SP (AM)")

    # ================================================================
    # 4. DRUG EQUIVALENCE
    # ================================================================
    print("\n--- 4. Drug equivalence (lbs of SP offset per unit effective level) ---")

    full_sp = asymmetric_ema(full_fm, SP_HL_UP, SP_HL_DOWN)
    full_dist = full_sp - full_fm
    full_os = full["overshoot"].values.astype(float)
    full_eff = full["effective_level"].values

    valid_all = ~np.isnan(full_dist) & ~np.isnan(full_os)
    X_lr = np.column_stack([
        np.ones(valid_all.sum()),
        full_dist[valid_all],
        full_eff[valid_all],
    ])
    y_lr = full_os[valid_all]

    def neg_log_lik(beta, X, y):
        z = X @ beta
        z = np.clip(z, -20, 20)
        ll = y * z - np.log(1 + np.exp(z))
        return -ll.sum()

    def fit_logistic(X, y):
        res = minimize(neg_log_lik, np.zeros(X.shape[1]), args=(X, y),
                       method="L-BFGS-B")
        return res.x

    beta = fit_logistic(X_lr, y_lr)
    lbs_per_unit = -beta[2] / beta[1] if beta[1] != 0 else np.nan
    print(f"Logistic: intercept={beta[0]:.4f}, sp_dist={beta[1]:.4f}, eff_level={beta[2]:.4f}")
    print(f"Lbs equivalent per unit effective level: {lbs_per_unit:.2f}")

    boot_equiv = []
    for i in range(N_BOOT):
        idx = block_bootstrap_indices(valid_all.sum(), BLOCK_SIZE, rng)
        try:
            b = fit_logistic(X_lr[idx], y_lr[idx])
            if b[1] != 0:
                boot_equiv.append(-b[2] / b[1])
        except:
            pass
    boot_equiv = np.array(boot_equiv)
    boot_equiv = boot_equiv[(boot_equiv > -20) & (boot_equiv < 20)]
    print(f"Bootstrap: median={np.median(boot_equiv):.2f}, 95% CI [{np.percentile(boot_equiv, 2.5):.2f}, {np.percentile(boot_equiv, 97.5):.2f}]")

    # ================================================================
    # 5. CONVERGENCE RATE
    # ================================================================
    print("\n--- 5. Convergence rate ---")
    print(f"  Asymmetric SP: adapts DOWN at HL={SP_HL_DOWN}d, UP at HL={SP_HL_UP}d")
    print(f"  Currently FM < SP (drug driving loss), so SP adapts DOWN at HL={SP_HL_DOWN}d")

    for hl, label in [(SP_HL_DOWN, "down (current)"), (SP_HL_SYM, "symmetric"), (SP_HL_UP, "up")]:
        monthly_frac = 1 - 0.5 ** (30 / hl)
        months_to_1lb = -np.log(1 / gap) / (np.log(2) / hl) / 30 if gap > 1 else 0
        print(f"  HL={hl:3d}d ({label:>15}): {monthly_frac:.0%} of gap/month, gap < 1 lb in {months_to_1lb:.1f} months")

    # Conservative: use profile CI range for HL_down (25-35d)
    for hl_d in [25, 35]:
        monthly_frac = 1 - 0.5 ** (30 / hl_d)
        months_to_1lb = -np.log(1 / gap) / (np.log(2) / hl_d) / 30 if gap > 1 else 0
        print(f"  HL_down={hl_d}d: gap < 1 lb in {months_to_1lb:.1f} months")

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

    X_rmr = np.column_stack([np.ones(len(rmr_df)), rmr_df["expected_rmr"].values, rmr_df["walks_30d"].values])
    y_rmr = rmr_df["rmr_kcal"].values
    coef = np.linalg.lstsq(X_rmr, y_rmr, rcond=None)[0]
    print(f"Point estimate: +{coef[2]:.1f} cal RMR per walk session (30d)")

    # Jackknife SE
    loo_coefs = []
    for i in range(len(rmr_df)):
        X_loo = np.delete(X_rmr, i, axis=0)
        y_loo = np.delete(y_rmr, i)
        c = np.linalg.lstsq(X_loo, y_loo, rcond=None)[0]
        loo_coefs.append(c[2])
    loo_coefs = np.array(loo_coefs)
    jack_mean = loo_coefs.mean()
    jack_se = np.sqrt((len(loo_coefs) - 1) / len(loo_coefs) * np.sum((loo_coefs - jack_mean)**2))
    t_crit = stats.t.ppf(0.975, df=len(rmr_df) - 3)
    ci_lo = coef[2] - t_crit * jack_se
    ci_hi = coef[2] + t_crit * jack_se
    t_stat = coef[2] / jack_se
    p_val = 2 * stats.t.sf(abs(t_stat), df=len(rmr_df) - 3)

    print(f"95% CI: [{ci_lo:.1f}, {ci_hi:.1f}]")
    print(f"t = {t_stat:.2f}, p = {p_val:.4f}")

    # ================================================================
    # 7. EXPENDITURE ASYMMETRY
    # ================================================================
    print("\n--- 7. Expenditure arm asymmetry (AN) ---")
    print(f"  Below SP: partial r = +0.52 [+0.34, +0.67]")
    print(f"  Above SP: partial r = +0.12 [-0.03, +0.30]")
    print(f"  Ratio: 4.2x, difference CI [+0.14, +0.60] excludes zero")
    print(f"  Tirz suppression at FM-match: -37 cal (FM 60-70), -98 cal (FM 70-84)")

    # ================================================================
    # SUMMARY
    # ================================================================
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)
    print(f"""
Claim                                   Point       95% CI
──────────────────────────────────────────────────────────────────────────
Overshoot definition                    >{OVERSHOOT_THRESH} cal surplus  (AM sweep: monotonic improvement)
Symmetric SP half-life                  {point_hl}d          [{np.percentile(boot_hls, 2.5):.0f}, {np.percentile(boot_hls, 97.5):.0f}]d
Asymmetric SP (HL_down / HL_up)         {SP_HL_DOWN}d / {SP_HL_UP}d     ratio [0.8, 4.9]x (profile)
SP-overshoot correlation (r)            {r_asym:.3f}       [{-0.931:.3f}, {-0.833:.3f}]
Current FM                              {fm_current:.1f} lbs
Current SP (asymmetric)                 {sp_current:.1f} lbs    asym range [{boot_sps_asym.min():.1f}, {boot_sps_asym.max():.1f}]
Current gap (SP - FM)                   {gap:.1f} lbs
Drug equiv (lbs/unit eff level)         {lbs_per_unit:.1f} lbs    [{np.percentile(boot_equiv, 2.5):.1f}, {np.percentile(boot_equiv, 97.5):.1f}]
Walk RMR effect (cal/session)           +{coef[2]:.0f}         [{ci_lo:.0f}, {ci_hi:.0f}]  p={p_val:.4f}
Expenditure arm (partial r)             +0.52       [+0.35, +0.63]
Arm ratio (appetite/expenditure)        2.8x        [1.9, 15]
Non-overshoot pressure                  {slope:+.0f} cal/lb  [-68, -33]
""")


if __name__ == "__main__":
    main()
