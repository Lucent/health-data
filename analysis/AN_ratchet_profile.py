#!/usr/bin/env python3
"""AN. Profile likelihood for the asymmetric ratchet ratio.

AM found the ratchet ratio CI was [0.5, 46] because bootstrap grid search
finds different optima on flat surfaces. This script uses profile likelihood
and continuous optimization to get a proper CI.

Method:
  1. For each fixed HL_up, optimize HL_down via scipy minimize_scalar
  2. Record the best r at each HL_up → likelihood profile
  3. CI from Fisher z-transform: pairs where Δ(Fisher z) < 1.96/√n_eff
  4. Also: expenditure arm asymmetry (new finding)
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parent.parent
BLOCK_SIZE = 90
N_BOOT = 1000


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


def block_bootstrap_indices(n, n_boot=N_BOOT, block_size=BLOCK_SIZE):
    bs = min(block_size, max(1, n // 2))
    n_blocks = max(1, n // bs)
    all_idx = []
    for _ in range(n_boot):
        starts = np.random.randint(0, n - bs + 1, size=n_blocks)
        idx = np.concatenate([np.arange(s, s + bs) for s in starts])[:n]
        all_idx.append(idx)
    return all_idx


def partial_corr(dist, outcome, fm):
    X = np.column_stack([fm, np.ones(len(fm))])
    res_d = dist - X @ np.linalg.lstsq(X, dist, rcond=None)[0]
    res_o = outcome - X @ np.linalg.lstsq(X, outcome, rcond=None)[0]
    return np.corrcoef(res_d, res_o)[0, 1]


def load_data():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    df = intake[["date", "calories"]].merge(
        kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="inner")
    df = df.merge(comp[["date", "expected_rmr", "ffm_lbs"]], on="date", how="left")
    df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df = df.sort_values("date").reset_index(drop=True)
    df["surplus"] = df["calories"] - df["tdee"]
    df["tdee_resid"] = df["tdee"] - df["expected_rmr"]
    return df


def main():
    np.random.seed(42)
    df = load_data()
    pre = df[df["effective_level"] == 0].copy().reset_index(drop=True)

    # Use optimized overshoot threshold from AM
    THRESH = 100
    WIN = 105
    pre["overshoot"] = (pre["calories"] > pre["tdee"] + THRESH).astype(float)
    pre["overshoot_rate"] = pre["overshoot"].rolling(WIN, min_periods=WIN).mean()

    fm = pre["fat_mass_lbs"].values
    br = pre["overshoot_rate"].values
    resid = pre["tdee_resid"].values

    # ═══════════════════════════════════════════════════════════════════
    print("=" * 70)
    print("1. PROFILE LIKELIHOOD for asymmetric half-lives")
    print("=" * 70)

    def neg_abs_r(hl_down, hl_up):
        """For a given HL_up, find the r at this HL_down."""
        if hl_down < 3 or hl_down > 200:
            return 0
        sp = asymmetric_ema(fm, hl_up, hl_down)
        dist = sp - fm
        valid = ~np.isnan(dist) & ~np.isnan(br)
        if valid.sum() < 300:
            return 0
        r = np.corrcoef(dist[valid], br[valid])[0, 1]
        return r  # more negative is better, so minimize r (= maximize |r|)

    # Profile: for each HL_up, find optimal HL_down and record r
    print(f"\n  {'HL_up':>6} {'HL_down_opt':>12} {'Ratio':>7} {'r':>8}")
    profile = []
    for hl_up in list(range(20, 60, 2)) + list(range(60, 150, 5)) + list(range(150, 301, 10)):
        res = minimize_scalar(lambda hd: neg_abs_r(hd, hl_up),
                              bounds=(3, 150), method='bounded',
                              options={'xatol': 0.5})
        best_hd = res.x
        best_r = res.fun  # this is r (negative)
        ratio = hl_up / best_hd if best_hd > 0 else np.inf
        profile.append({"hl_up": hl_up, "hl_down": best_hd, "ratio": ratio, "r": best_r})
        if hl_up % 10 == 0 or hl_up < 40:
            print(f"  {hl_up:>4}d {best_hd:>10.1f}d {ratio:>7.1f}x {best_r:+8.4f}")

    profile_df = pd.DataFrame(profile)

    # Find optimum
    best_idx = profile_df["r"].idxmin()
    best = profile_df.iloc[best_idx]
    print(f"\n  Optimum: HL_up={best['hl_up']:.0f}d, HL_down={best['hl_down']:.1f}d, ratio={best['ratio']:.1f}x, r={best['r']:+.4f}")

    # Fisher z-transform CI
    # The CI comes from: which (hl_up, hl_down) pairs have r not significantly
    # worse than the best? Use Fisher z: z = arctanh(r), var(z) ≈ 1/n_eff
    # Compute effective n from autocorrelation
    sp_opt = asymmetric_ema(fm, best["hl_up"], best["hl_down"])
    dist_opt = sp_opt - fm
    valid_opt = ~np.isnan(dist_opt) & ~np.isnan(br)
    product = dist_opt[valid_opt] * br[valid_opt]
    product_c = product - product.mean()
    acf1 = np.corrcoef(product_c[:-1], product_c[1:])[0, 1]
    n_raw = valid_opt.sum()
    n_eff = n_raw * (1 - acf1) / (1 + acf1) if 0 < acf1 < 1 else n_raw
    n_blocks = n_raw / BLOCK_SIZE

    print(f"\n  N raw: {n_raw}, ACF1: {acf1:.3f}, N_eff (Bartlett): {n_eff:.0f}, N_blocks: {n_blocks:.0f}")

    # For profile likelihood CI: Δz = 1.96 / √n_eff
    z_best = np.arctanh(abs(best["r"]))
    z_threshold = z_best - 1.96 / np.sqrt(n_eff)  # 95% CI in Fisher z
    r_threshold = np.tanh(z_threshold)  # back to r

    in_ci = profile_df[profile_df["r"].abs() >= r_threshold]
    if len(in_ci) > 0:
        ci_ratio_lo = in_ci["ratio"].min()
        ci_ratio_hi = in_ci["ratio"].max()
        ci_up_lo = in_ci["hl_up"].min()
        ci_up_hi = in_ci["hl_up"].max()
        ci_down_lo = in_ci["hl_down"].min()
        ci_down_hi = in_ci["hl_down"].max()
    else:
        ci_ratio_lo = ci_ratio_hi = best["ratio"]
        ci_up_lo = ci_up_hi = best["hl_up"]
        ci_down_lo = ci_down_hi = best["hl_down"]

    print(f"\n  Fisher z threshold: r >= {r_threshold:.4f} (Δz = 1.96/√{n_eff:.0f} = {1.96 / np.sqrt(n_eff):.4f})")
    print(f"  HL_up CI:   [{ci_up_lo:.0f}, {ci_up_hi:.0f}]")
    print(f"  HL_down CI: [{ci_down_lo:.1f}, {ci_down_hi:.1f}]")
    print(f"  Ratio CI:   [{ci_ratio_lo:.1f}, {ci_ratio_hi:.1f}]")

    # More conservative: use n_blocks instead of n_eff
    z_threshold_b = z_best - 1.96 / np.sqrt(n_blocks)
    r_threshold_b = np.tanh(z_threshold_b)
    in_ci_b = profile_df[profile_df["r"].abs() >= r_threshold_b]
    if len(in_ci_b) > 0:
        print(f"\n  Conservative (N_blocks={n_blocks:.0f}):")
        print(f"  HL_up CI:   [{in_ci_b['hl_up'].min():.0f}, {in_ci_b['hl_up'].max():.0f}]")
        print(f"  HL_down CI: [{in_ci_b['hl_down'].min():.1f}, {in_ci_b['hl_down'].max():.1f}]")
        print(f"  Ratio CI:   [{in_ci_b['ratio'].min():.1f}, {in_ci_b['ratio'].max():.1f}]")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("2. BOOTSTRAP with continuous optimization (tighter than grid)")
    print("=" * 70)

    boot_idx = block_bootstrap_indices(len(pre), 500, BLOCK_SIZE)
    boot_ratios = []
    boot_ups = []
    boot_downs = []
    for idx in boot_idx:
        boot_fm = fm[idx]
        boot_br = br[idx]
        valid_b = ~np.isnan(boot_br)
        if valid_b.sum() < 200:
            continue

        def neg_r_boot(params):
            hu, hd = params
            if hu < 10 or hd < 3 or hu > 300 or hd > 150:
                return 0
            sp = asymmetric_ema(boot_fm, hu, hd)
            dist = sp - boot_fm
            vb = ~np.isnan(dist) & valid_b
            if vb.sum() < 100:
                return 0
            return np.corrcoef(dist[vb], boot_br[vb])[0, 1]

        # Coarse grid to find starting point, then refine
        b_best_r, b_best_u, b_best_d = 0, 50, 50
        for hu in range(30, 200, 30):
            for hd in range(5, 60, 15):
                r = neg_r_boot([hu, hd])
                if r < b_best_r:
                    b_best_r, b_best_u, b_best_d = r, hu, hd

        # Refine with continuous optimization around the best
        for hu in range(max(10, b_best_u - 20), b_best_u + 21, 5):
            res = minimize_scalar(lambda hd: neg_r_boot([hu, hd]),
                                  bounds=(3, 100), method='bounded',
                                  options={'xatol': 1.0})
            if res.fun < b_best_r:
                b_best_r, b_best_u, b_best_d = res.fun, hu, res.x

        boot_ups.append(b_best_u)
        boot_downs.append(b_best_d)
        boot_ratios.append(b_best_u / b_best_d if b_best_d > 0 else np.nan)

    boot_ratios = np.array([x for x in boot_ratios if not np.isnan(x)])
    boot_ups = np.array(boot_ups)
    boot_downs = np.array(boot_downs)

    print(f"\n  {len(boot_ratios)} valid bootstrap resamples")
    print(f"  HL_up:   {best['hl_up']:.0f}d [{np.percentile(boot_ups, 2.5):.0f}, {np.percentile(boot_ups, 97.5):.0f}]")
    print(f"  HL_down: {best['hl_down']:.1f}d [{np.percentile(boot_downs, 2.5):.0f}, {np.percentile(boot_downs, 97.5):.0f}]")
    print(f"  Ratio:   {best['ratio']:.1f}x [{np.percentile(boot_ratios, 2.5):.1f}, {np.percentile(boot_ratios, 97.5):.1f}]")
    print(f"  Median ratio: {np.median(boot_ratios):.1f}x")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("3. EXPENDITURE ARM ASYMMETRY — does TDEE defense differ above vs below SP?")
    print("=" * 70)

    # For each HL, split: FM > SP (above) vs FM < SP (below)
    # Compare partial r of TDEE residual vs distance in each half
    print(f"\n  {'HL':>4} {'r_below':>9} {'r_above':>9} {'r_all':>7} {'ratio':>7}")
    for hl in [5, 8, 10, 15, 20, 30, 50, 80, 100]:
        sp = ema(pre["fat_mass_lbs"], hl).values
        dist = sp - fm
        valid = ~np.isnan(dist) & ~np.isnan(resid) & ~np.isnan(fm)

        below = valid & (dist > 0)  # FM below SP → dist positive
        above = valid & (dist <= 0)  # FM above SP → dist negative/zero

        r_all = partial_corr(dist[valid], resid[valid], fm[valid]) if valid.sum() > 200 else np.nan
        r_below = partial_corr(dist[below], resid[below], fm[below]) if below.sum() > 100 else np.nan
        r_above = partial_corr(dist[above], resid[above], fm[above]) if above.sum() > 100 else np.nan

        ratio_str = f"{abs(r_below) / abs(r_above):.1f}x" if (not np.isnan(r_below) and not np.isnan(r_above) and abs(r_above) > 0.01) else "n/a"
        print(f"  {hl:>4} {r_below:+9.3f} {r_above:+9.3f} {r_all:+7.3f} {ratio_str:>7}" if not np.isnan(r_all) else
              f"  {hl:>4}  {'n/a':>8}  {'n/a':>8}  {'n/a':>6}  {'n/a':>6}")

    # Bootstrap CI for the above/below asymmetry at HL=10
    print(f"\n  --- Bootstrap CI for expenditure asymmetry at HL=10d ---")
    sp10 = ema(pre["fat_mass_lbs"], 10).values
    dist10 = sp10 - fm
    valid10 = ~np.isnan(dist10) & ~np.isnan(resid) & ~np.isnan(fm)

    below10 = valid10 & (dist10 > 0)
    above10 = valid10 & (dist10 <= 0)

    r_below_pt = partial_corr(dist10[below10], resid[below10], fm[below10])
    r_above_pt = partial_corr(dist10[above10], resid[above10], fm[above10])

    n_below = below10.sum()
    n_above = above10.sum()
    print(f"  Below SP (n={n_below}): partial r = {r_below_pt:+.3f}")
    print(f"  Above SP (n={n_above}): partial r = {r_above_pt:+.3f}")

    # Bootstrap
    below_boots, above_boots, diff_boots = [], [], []
    boot_idx_b = block_bootstrap_indices(n_below, 500, BLOCK_SIZE)
    boot_idx_a = block_bootstrap_indices(n_above, 500, BLOCK_SIZE)

    d_below = dist10[below10]
    r_below_arr = resid[below10]
    f_below = fm[below10]
    d_above = dist10[above10]
    r_above_arr = resid[above10]
    f_above = fm[above10]

    for i in range(500):
        rb = partial_corr(d_below[boot_idx_b[i]], r_below_arr[boot_idx_b[i]], f_below[boot_idx_b[i]])
        ra = partial_corr(d_above[boot_idx_a[i]], r_above_arr[boot_idx_a[i]], f_above[boot_idx_a[i]])
        below_boots.append(rb)
        above_boots.append(ra)
        diff_boots.append(abs(rb) - abs(ra))

    below_boots = np.array(below_boots)
    above_boots = np.array(above_boots)
    diff_boots = np.array(diff_boots)

    print(f"  Below SP: {r_below_pt:+.3f} [{np.percentile(below_boots, 2.5):+.3f}, {np.percentile(below_boots, 97.5):+.3f}]")
    print(f"  Above SP: {r_above_pt:+.3f} [{np.percentile(above_boots, 2.5):+.3f}, {np.percentile(above_boots, 97.5):+.3f}]")
    print(f"  |below| - |above|: {abs(r_below_pt) - abs(r_above_pt):+.3f} [{np.percentile(diff_boots, 2.5):+.3f}, {np.percentile(diff_boots, 97.5):+.3f}]")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("4. APPETITE ARM ASYMMETRY — does overshoot response differ above vs below SP?")
    print("=" * 70)

    for hl in [30, 43, 51, 70, 100]:
        sp = ema(pre["fat_mass_lbs"], hl).values
        dist = sp - fm
        valid = ~np.isnan(dist) & ~np.isnan(br) & ~np.isnan(fm)

        below = valid & (dist > 0)
        above = valid & (dist <= 0)

        r_all = partial_corr(dist[valid], br[valid], fm[valid]) if valid.sum() > 200 else np.nan
        r_below = partial_corr(dist[below], br[below], fm[below]) if below.sum() > 100 else np.nan
        r_above = partial_corr(dist[above], br[above], fm[above]) if above.sum() > 100 else np.nan

        print(f"  HL={hl:>3}d: below r={r_below:+.3f} (n={below.sum()}), above r={r_above:+.3f} (n={above.sum()})")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("5. CROSS-VALIDATION: does the model predict held-out years?")
    print("=" * 70)

    years = sorted(pre["date"].dt.year.unique())
    print(f"\n  {'Held-out':>10} {'r_train':>8} {'r_test':>7} {'HL_up':>6} {'HL_down':>8}")
    for hold_yr in years:
        if hold_yr < 2013 or hold_yr > 2024:
            continue
        train = pre[pre["date"].dt.year != hold_yr].reset_index(drop=True)
        test = pre[pre["date"].dt.year == hold_yr].reset_index(drop=True)

        train_fm = train["fat_mass_lbs"].values
        train_br = train["overshoot"].rolling(WIN, min_periods=WIN).mean().values

        # Find best HL on train
        b_r, b_u, b_d = 0, 50, 50
        for hu in range(30, 200, 15):
            for hd in range(5, 60, 10):
                sp = asymmetric_ema(train_fm, hu, hd)
                dist = sp - train_fm
                v = ~np.isnan(dist) & ~np.isnan(train_br)
                if v.sum() < 200:
                    continue
                r = np.corrcoef(dist[v], train_br[v])[0, 1]
                if r < b_r:
                    b_r, b_u, b_d = r, hu, hd

        # Apply to full data and extract test year
        full_fm = pre["fat_mass_lbs"].values
        full_br = pre["overshoot_rate"].values
        sp_full = asymmetric_ema(full_fm, b_u, b_d)
        dist_full = sp_full - full_fm

        test_mask = pre["date"].dt.year == hold_yr
        test_dist = dist_full[test_mask.values]
        test_br = full_br[test_mask.values]
        v_test = ~np.isnan(test_dist) & ~np.isnan(test_br)
        r_test = np.corrcoef(test_dist[v_test], test_br[v_test])[0, 1] if v_test.sum() > 30 else np.nan

        print(f"  {hold_yr:>10} {b_r:+8.3f} {r_test:+7.3f} {b_u:>4}d {b_d:>6}d")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)
    print(f"\n  Profile likelihood ratio CI: [{ci_ratio_lo:.1f}, {ci_ratio_hi:.1f}]")
    if len(in_ci_b) > 0:
        print(f"  Conservative ratio CI:       [{in_ci_b['ratio'].min():.1f}, {in_ci_b['ratio'].max():.1f}]")
    print(f"  Bootstrap (continuous opt):   [{np.percentile(boot_ratios, 2.5):.1f}, {np.percentile(boot_ratios, 97.5):.1f}]")
    print(f"  Point estimate:              {best['ratio']:.1f}x (HL_up={best['hl_up']:.0f}d, HL_down={best['hl_down']:.1f}d)")
    print(f"  Expenditure asymmetry:       below r={r_below_pt:+.3f}, above r={r_above_pt:+.3f}")


if __name__ == "__main__":
    main()
