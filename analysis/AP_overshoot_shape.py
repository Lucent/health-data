#!/usr/bin/env python3
"""AP. What functional form best captures the set point → intake relationship?

AM found binary overshoot at 100 cal beats 1000 cal. But this treats all
overshoots equally. Test whether magnitude-weighted measures do better:

1. Binary threshold (AM baseline)
2. Continuous mean surplus (no threshold)
3. Censored surplus: max(surplus - threshold, 0) — magnitude above threshold
4. Mean surplus clipped to positive: mean(max(surplus, 0))
5. Surplus percentiles (75th, 90th) — captures right tail
6. Log-surplus: mean(log(1 + max(surplus, 0)))
7. Two-component: separate rate and magnitude channels
8. Quadratic: mean surplus + mean surplus²
"""

from pathlib import Path
import numpy as np
import pandas as pd

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


def partial_corr(dist, outcome, fm):
    X = np.column_stack([fm, np.ones(len(fm))])
    res_d = dist - X @ np.linalg.lstsq(X, dist, rcond=None)[0]
    res_o = outcome - X @ np.linalg.lstsq(X, outcome, rcond=None)[0]
    return np.corrcoef(res_d, res_o)[0, 1]


def block_bootstrap_indices(n, n_boot=N_BOOT, block_size=BLOCK_SIZE):
    bs = min(block_size, max(1, n // 2))
    n_blocks = max(1, n // bs)
    all_idx = []
    for _ in range(n_boot):
        starts = np.random.randint(0, n - bs + 1, size=n_blocks)
        idx = np.concatenate([np.arange(s, s + bs) for s in starts])[:n]
        all_idx.append(idx)
    return all_idx


def block_bootstrap_corr(x, y, n_boot=N_BOOT, block_size=BLOCK_SIZE):
    indices = block_bootstrap_indices(len(x), n_boot, block_size)
    rs = np.array([np.corrcoef(x[idx], y[idx])[0, 1] for idx in indices])
    return np.percentile(rs, [2.5, 97.5])


def main():
    np.random.seed(42)
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    df = intake[["date", "calories"]].merge(
        kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="inner")
    df = df.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df = df.sort_values("date").reset_index(drop=True)
    df["surplus"] = df["calories"] - df["tdee"]

    pre = df[df["effective_level"] == 0].copy().reset_index(drop=True)
    n = len(pre)

    fm = pre["fat_mass_lbs"].values
    surplus = pre["surplus"].values

    # Set point (asymmetric from AN)
    sp_asym = asymmetric_ema(fm, 72, 25)
    dist_asym = sp_asym - fm  # positive = below SP
    # Also symmetric for comparison
    sp_sym = ema(pre["fat_mass_lbs"], 51).values

    # ═══════════════════════════════════════════════════════════════════
    print("=" * 70)
    print("1. FUNCTIONAL FORMS — which intake summary best correlates with SP distance?")
    print("=" * 70)

    def test_metric(name, daily_values, windows=[60, 90, 105, 120]):
        """Compute rolling mean of daily_values at each window, correlate with SP distance."""
        best_r, best_win = 0, 90
        for win in windows:
            rolled = pd.Series(daily_values).rolling(win, min_periods=win).mean().values
            valid = ~np.isnan(dist_asym) & ~np.isnan(rolled)
            if valid.sum() < 300:
                continue
            r = np.corrcoef(dist_asym[valid], rolled[valid])[0, 1]
            if abs(r) > abs(best_r):
                best_r, best_win = r, win
        return name, best_r, best_win

    # Define all metrics
    metrics = []

    # Binary thresholds
    for thresh in [0, 100, 200, 500, 1000]:
        daily = (surplus > thresh).astype(float)
        metrics.append(test_metric(f"binary >{thresh}", daily))

    # Continuous mean surplus (raw)
    metrics.append(test_metric("mean surplus", surplus))

    # Positive surplus only: max(surplus, 0)
    pos_surplus = np.maximum(surplus, 0)
    metrics.append(test_metric("mean pos surplus", pos_surplus))

    # Negative surplus only: min(surplus, 0)
    neg_surplus = np.minimum(surplus, 0)
    metrics.append(test_metric("mean neg surplus", neg_surplus))

    # Censored surplus above thresholds: max(surplus - thresh, 0)
    for thresh in [0, 100, 500, 1000]:
        censored = np.maximum(surplus - thresh, 0)
        metrics.append(test_metric(f"cens surplus >{thresh}", censored))

    # Log positive surplus
    log_pos = np.log1p(np.maximum(surplus, 0))
    metrics.append(test_metric("log(1+pos_surplus)", log_pos))

    # Sqrt positive surplus
    sqrt_pos = np.sqrt(np.maximum(surplus, 0))
    metrics.append(test_metric("sqrt(pos_surplus)", sqrt_pos))

    # Absolute surplus (both directions)
    abs_surplus = np.abs(surplus)
    metrics.append(test_metric("mean |surplus|", abs_surplus))

    # Surplus percentiles (rolling)
    for pct in [75, 90, 95]:
        for win in [90, 105]:
            rolled = pd.Series(surplus).rolling(win, min_periods=win).quantile(pct / 100).values
            valid = ~np.isnan(dist_asym) & ~np.isnan(rolled)
            if valid.sum() > 300:
                r = np.corrcoef(dist_asym[valid], rolled[valid])[0, 1]
                metrics.append((f"P{pct} surplus ({win}d)", r, win))

    # Surplus variance (rolling std)
    for win in [60, 90, 120]:
        rolled = pd.Series(surplus).rolling(win, min_periods=win).std().values
        valid = ~np.isnan(dist_asym) & ~np.isnan(rolled)
        if valid.sum() > 300:
            r = np.corrcoef(dist_asym[valid], rolled[valid])[0, 1]
            metrics.append((f"surplus std ({win}d)", r, win))

    # Sort by |r|
    metrics.sort(key=lambda x: -abs(x[1]))
    print(f"\n  {'Metric':>30} {'r':>8} {'Window':>8}")
    for name, r, win in metrics[:25]:
        print(f"  {name:>30} {r:+8.4f} {win:>6}d")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("2. HEAD-TO-HEAD: binary 100 vs continuous vs magnitude-weighted")
    print("=" * 70)

    WIN = 105  # fixed for comparison

    # Binary 100
    binary_100 = pd.Series((surplus > 100).astype(float)).rolling(WIN, min_periods=WIN).mean().values
    # Mean surplus
    mean_surp = pd.Series(surplus).rolling(WIN, min_periods=WIN).mean().values
    # Positive surplus
    mean_pos = pd.Series(pos_surplus).rolling(WIN, min_periods=WIN).mean().values
    # Log positive
    mean_log = pd.Series(log_pos).rolling(WIN, min_periods=WIN).mean().values
    # Sqrt positive
    mean_sqrt = pd.Series(sqrt_pos).rolling(WIN, min_periods=WIN).mean().values

    valid = ~np.isnan(dist_asym) & ~np.isnan(binary_100) & ~np.isnan(mean_surp) & ~np.isnan(fm)
    print(f"\n  All at {WIN}d window, asymmetric SP:")
    print(f"  {'Metric':>25} {'r':>8} {'partial r':>10} {'CI':>18}")

    for label, arr in [("binary >100", binary_100),
                       ("mean surplus", mean_surp),
                       ("mean pos surplus", mean_pos),
                       ("log(1+pos)", mean_log),
                       ("sqrt(pos)", mean_sqrt)]:
        valid_i = valid & ~np.isnan(arr)
        if valid_i.sum() < 300:
            continue
        r = np.corrcoef(dist_asym[valid_i], arr[valid_i])[0, 1]
        pr = partial_corr(dist_asym[valid_i], arr[valid_i], fm[valid_i])
        ci = block_bootstrap_corr(dist_asym[valid_i], arr[valid_i])
        print(f"  {label:>25} {r:+8.4f} {pr:+10.4f} [{ci[0]:+.3f}, {ci[1]:+.3f}]")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("3. TWO-COMPONENT MODEL: rate + magnitude as separate channels")
    print("=" * 70)

    # Can we decompose the signal into frequency and magnitude?
    # On overshoot days: does the magnitude of overshoot also correlate with SP distance?
    print(f"\n  --- Magnitude of overshoot days ---")
    for thresh in [0, 100, 500]:
        overshoot_mask = surplus > thresh
        overshoot_mag = np.where(overshoot_mask, surplus - thresh, np.nan)
        # Rolling mean of magnitude (only on overshoot days)
        # Use: total surplus above thresh / total days = rate × mean_magnitude
        total_excess = pd.Series(np.maximum(surplus - thresh, 0)).rolling(WIN, min_periods=WIN).mean().values
        # Decompose: total = rate × conditional_mean
        rate = pd.Series(overshoot_mask.astype(float)).rolling(WIN, min_periods=WIN).mean().values

        valid_both = ~np.isnan(dist_asym) & ~np.isnan(total_excess) & ~np.isnan(rate) & (rate > 0)
        if valid_both.sum() < 300:
            continue
        conditional_mean = total_excess[valid_both] / rate[valid_both]

        r_total = np.corrcoef(dist_asym[valid_both], total_excess[valid_both])[0, 1]
        r_rate = np.corrcoef(dist_asym[valid_both], rate[valid_both])[0, 1]
        r_cond = np.corrcoef(dist_asym[valid_both], conditional_mean)[0, 1]

        print(f"  Threshold {thresh}: total r={r_total:+.4f}, rate r={r_rate:+.4f}, cond.mean r={r_cond:+.4f}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("4. OPTIMAL HALF-LIFE for each metric")
    print("   Does the best metric change which HL wins?")
    print("=" * 70)

    print(f"\n  {'Metric':>25} {'Best HL':>8} {'r':>8}")
    for label, daily in [("binary >100", (surplus > 100).astype(float)),
                         ("binary >1000", (surplus > 1000).astype(float)),
                         ("mean surplus", surplus),
                         ("mean pos surplus", pos_surplus),
                         ("sqrt(pos surplus)", sqrt_pos),
                         ("log(1+pos surplus)", log_pos)]:
        best_r, best_hl = 0, 50
        for hl in range(15, 200, 2):
            sp = ema(pre["fat_mass_lbs"], hl)
            dist = (sp - pre["fat_mass_lbs"]).values
            rolled = pd.Series(daily).rolling(WIN, min_periods=WIN).mean().values
            valid_i = ~np.isnan(dist) & ~np.isnan(rolled)
            if valid_i.sum() < 300:
                continue
            r = np.corrcoef(dist[valid_i], rolled[valid_i])[0, 1]
            if abs(r) > abs(best_r):
                best_r, best_hl = r, hl
        print(f"  {label:>25} {best_hl:>6}d {best_r:+8.4f}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("5. PARTIAL R² — does magnitude add signal beyond rate?")
    print("=" * 70)

    # Regression: SP_distance ~ rate + conditional_magnitude
    # If magnitude adds signal, the R² increases
    rate_105 = pd.Series((surplus > 100).astype(float)).rolling(WIN, min_periods=WIN).mean().values
    total_105 = pd.Series(pos_surplus).rolling(WIN, min_periods=WIN).mean().values

    valid_reg = ~np.isnan(dist_asym) & ~np.isnan(rate_105) & ~np.isnan(total_105) & ~np.isnan(fm)
    d = dist_asym[valid_reg]
    r_v = rate_105[valid_reg]
    t_v = total_105[valid_reg]
    f_v = fm[valid_reg]

    # Model 1: dist ~ rate only
    X1 = np.column_stack([r_v, f_v, np.ones(valid_reg.sum())])
    beta1 = np.linalg.lstsq(X1, d, rcond=None)[0]
    resid1 = d - X1 @ beta1
    r2_1 = 1 - np.sum(resid1**2) / np.sum((d - d.mean())**2)

    # Model 2: dist ~ total (= rate × magnitude)
    X2 = np.column_stack([t_v, f_v, np.ones(valid_reg.sum())])
    beta2 = np.linalg.lstsq(X2, d, rcond=None)[0]
    resid2 = d - X2 @ beta2
    r2_2 = 1 - np.sum(resid2**2) / np.sum((d - d.mean())**2)

    # Model 3: dist ~ rate + total (both)
    X3 = np.column_stack([r_v, t_v, f_v, np.ones(valid_reg.sum())])
    beta3 = np.linalg.lstsq(X3, d, rcond=None)[0]
    resid3 = d - X3 @ beta3
    r2_3 = 1 - np.sum(resid3**2) / np.sum((d - d.mean())**2)

    print(f"\n  Predicting SP distance (controlling for FM):")
    print(f"  Rate only (binary >100):       R² = {r2_1:.4f}")
    print(f"  Total pos surplus:             R² = {r2_2:.4f}")
    print(f"  Rate + total (both):           R² = {r2_3:.4f}")
    print(f"  Marginal R² from adding total: {r2_3 - r2_1:.4f}")

    # Also: does the REVERSE direction work? (SP distance predicts intake measures)
    print(f"\n  Predicting intake measures from SP distance (controlling for FM):")
    for label, y_var in [("rate >100", r_v), ("mean pos surplus", t_v)]:
        X = np.column_stack([d, f_v, np.ones(len(d))])
        beta = np.linalg.lstsq(X, y_var, rcond=None)[0]
        pred = X @ beta
        r2 = 1 - np.sum((y_var - pred)**2) / np.sum((y_var - y_var.mean())**2)
        print(f"  dist → {label:>20}: R² = {r2:.4f}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)


if __name__ == "__main__":
    main()
