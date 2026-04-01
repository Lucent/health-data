#!/usr/bin/env python3
"""AO. Tirzepatide specifically suppresses the below-SP expenditure defense.

AN found the expenditure arm is 4.2x stronger when below SP than above.
AJ found tirzepatide suppresses the overall expenditure defense.
This script tests: does the drug selectively suppress the below-SP defense,
leaving the already-weak above-SP response unchanged?

If true, this explains the mechanism: the drug doesn't just globally reduce
metabolic rate — it specifically disables the body's strongest defense
against fat loss, the below-SP TDEE elevation.
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
    df["tdee_resid"] = df["tdee"] - df["expected_rmr"]
    df["tdee_ratio"] = df["tdee"] / df["expected_rmr"]
    df["on_tirz"] = (df["effective_level"] > 0).astype(int)

    fm = df["fat_mass_lbs"].values
    resid = df["tdee_resid"].values
    ratio = df["tdee_ratio"].values

    # Compute SP distance at HL=10d (expenditure arm timescale)
    sp10 = ema(df["fat_mass_lbs"], 10).values
    dist10 = sp10 - fm  # positive = FM below SP (losing weight)

    valid = ~np.isnan(dist10) & ~np.isnan(resid) & ~np.isnan(fm)
    pre = df["on_tirz"].values == 0
    on = df["on_tirz"].values == 1

    # ═══════════════════════════════════════════════════════════════════
    print("=" * 70)
    print("1. EXPENDITURE DEFENSE: pre-tirz vs on-tirz, by SP position")
    print("=" * 70)

    # Split into 4 groups: (pre/on tirz) × (below/above SP)
    groups = {
        "Pre-tirz, below SP": valid & pre & (dist10 > 0),
        "Pre-tirz, above SP": valid & pre & (dist10 <= 0),
        "On-tirz, below SP":  valid & on & (dist10 > 0),
        "On-tirz, above SP":  valid & on & (dist10 <= 0),
    }

    print(f"\n  {'Group':>25} {'n':>6} {'TDEE-RMR':>9} {'TDEE/RMR':>9} {'partial r':>10}")
    for label, mask in groups.items():
        n_g = mask.sum()
        if n_g < 30:
            print(f"  {label:>25} {n_g:>6}  too few")
            continue
        mean_resid = resid[mask].mean()
        mean_ratio = ratio[mask].mean()
        r_partial = partial_corr(dist10[mask], resid[mask], fm[mask]) if n_g > 50 else np.nan
        print(f"  {label:>25} {n_g:>6} {mean_resid:+9.0f} {mean_ratio:9.4f} {r_partial:+10.3f}" if not np.isnan(r_partial) else
              f"  {label:>25} {n_g:>6} {mean_resid:+9.0f} {mean_ratio:9.4f}  {'n/a':>9}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("2. BINNED TDEE/RMR by SP distance, pre-tirz vs on-tirz")
    print("=" * 70)

    bins = [(-15, -5), (-5, -2.5), (-2.5, 0), (0, 2.5), (2.5, 5), (5, 10), (10, 20)]
    print(f"\n  {'Distance':>15} {'Pre-tirz':>12} {'On-tirz':>12} {'Δ':>8}")
    for lo, hi in bins:
        mask_pre = valid & pre & (dist10 > lo) & (dist10 <= hi)
        mask_on = valid & on & (dist10 > lo) & (dist10 <= hi)
        n_pre = mask_pre.sum()
        n_on = mask_on.sum()
        if n_pre < 20:
            continue
        r_pre = ratio[mask_pre].mean()
        r_on = ratio[mask_on].mean() if n_on >= 10 else np.nan
        delta = r_on - r_pre if not np.isnan(r_on) else np.nan
        on_str = f"{r_on:12.4f}" if not np.isnan(r_on) else f"{'n/a (n=' + str(n_on) + ')':>12}"
        d_str = f"{delta:+8.4f}" if not np.isnan(delta) else f"{'':>8}"
        print(f"  {lo:+.0f} to {hi:+.0f} {r_pre:12.4f} {on_str} {d_str}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("3. REGRESSION: does tirz × below_SP interaction explain defense suppression?")
    print("=" * 70)

    # Model: TDEE_resid ~ FM + below_SP + on_tirz + below_SP × on_tirz
    mask = valid
    below_sp = (dist10[mask] > 0).astype(float)
    on_tirz_v = df["on_tirz"].values[mask].astype(float)
    fm_v = fm[mask]
    resid_v = resid[mask]
    dist_v = dist10[mask]

    # Full regression
    X = np.column_stack([
        fm_v,
        dist_v,
        on_tirz_v,
        dist_v * on_tirz_v,
        np.ones(mask.sum()),
    ])
    labels = ["FM", "SP_distance", "on_tirz", "dist×tirz", "intercept"]
    beta = np.linalg.lstsq(X, resid_v, rcond=None)[0]
    pred = X @ beta
    ss_res = np.sum((resid_v - pred) ** 2)
    ss_tot = np.sum((resid_v - resid_v.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot

    print(f"\n  TDEE_resid ~ FM + SP_distance + on_tirz + distance×tirz")
    print(f"  R² = {r2:.4f}, n = {mask.sum()}")
    print(f"\n  {'Feature':>15} {'Coefficient':>12}")
    for l, b in zip(labels, beta):
        print(f"  {l:>15} {b:+12.3f}")

    print(f"\n  Interpretation:")
    print(f"  SP_distance: each lb below SP adds {beta[1]:+.1f} cal to TDEE (pre-tirz)")
    print(f"  dist×tirz:   on tirz, this effect changes by {beta[3]:+.1f} cal/lb")
    print(f"  Net on tirz: {beta[1] + beta[3]:+.1f} cal per lb below SP")
    if abs(beta[1]) > 0:
        print(f"  Drug suppresses {abs(beta[3]) / abs(beta[1]) * 100:.0f}% of the per-lb TDEE defense")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("4. FM-MATCHED COMPARISON — same fat mass, different drug status")
    print("=" * 70)

    # Find FM overlap between pre and on-tirz
    pre_fm = fm[valid & pre]
    on_fm = fm[valid & on]
    overlap_lo = max(pre_fm.min(), on_fm.min())
    overlap_hi = min(pre_fm.max(), on_fm.max())
    print(f"\n  FM overlap: {overlap_lo:.0f} to {overlap_hi:.0f} lbs")

    # Within overlap, compare below-SP TDEE residual
    for band_lo, band_hi in [(overlap_lo, overlap_lo + 10), (overlap_lo + 10, overlap_hi)]:
        mask_pre_band = valid & pre & (fm >= band_lo) & (fm < band_hi) & (dist10 > 0)
        mask_on_band = valid & on & (fm >= band_lo) & (fm < band_hi) & (dist10 > 0)
        n_p = mask_pre_band.sum()
        n_o = mask_on_band.sum()
        if n_p < 20 or n_o < 10:
            print(f"  FM {band_lo:.0f}-{band_hi:.0f}: pre n={n_p}, on n={n_o} — too few")
            continue
        r_pre = resid[mask_pre_band].mean()
        r_on = resid[mask_on_band].mean()
        fm_pre = fm[mask_pre_band].mean()
        fm_on = fm[mask_on_band].mean()
        d_pre = dist10[mask_pre_band].mean()
        d_on = dist10[mask_on_band].mean()
        print(f"  FM {band_lo:.0f}-{band_hi:.0f}:")
        print(f"    Pre-tirz (n={n_p}): mean FM={fm_pre:.1f}, dist={d_pre:+.1f}, TDEE-RMR={r_pre:+.0f}")
        print(f"    On-tirz  (n={n_o}): mean FM={fm_on:.1f}, dist={d_on:+.1f}, TDEE-RMR={r_on:+.0f}")
        print(f"    Δ TDEE-RMR: {r_on - r_pre:+.0f} cal")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("5. APPETITE ARM: does tirz suppress below-SP overshoot differently?")
    print("=" * 70)

    sp50 = ema(df["fat_mass_lbs"], 51).values
    dist50 = sp50 - fm
    THRESH = 100
    WIN = 105
    df["overshoot"] = (df["calories"] > df["tdee"] + THRESH).astype(float)
    df["overshoot_rate"] = df["overshoot"].rolling(WIN, min_periods=WIN).mean()
    or_v = df["overshoot_rate"].values

    valid_a = ~np.isnan(dist50) & ~np.isnan(or_v)
    for label, mask_era in [("Pre-tirz", valid_a & pre), ("On-tirz", valid_a & on)]:
        below = mask_era & (dist50 > 0)
        above = mask_era & (dist50 <= 0)
        r_below = np.corrcoef(dist50[below], or_v[below])[0, 1] if below.sum() > 50 else np.nan
        r_above = np.corrcoef(dist50[above], or_v[above])[0, 1] if above.sum() > 50 else np.nan
        or_below = or_v[below].mean() if below.sum() > 50 else np.nan
        or_above = or_v[above].mean() if above.sum() > 50 else np.nan
        print(f"  {label}:")
        print(f"    Below SP (n={below.sum()}): overshoot rate = {or_below:.1%}, r = {r_below:+.3f}" if not np.isnan(r_below) else
              f"    Below SP (n={below.sum()}): too few")
        print(f"    Above SP (n={above.sum()}): overshoot rate = {or_above:.1%}, r = {r_above:+.3f}" if not np.isnan(r_above) else
              f"    Above SP (n={above.sum()}): too few")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)
    print(f"\n  The expenditure defense (TDEE-RMR) is strongly asymmetric:")
    print(f"  Pre-tirz below SP: TDEE-RMR = {resid[valid & pre & (dist10 > 0)].mean():+.0f} cal (defending)")
    print(f"  Pre-tirz above SP: TDEE-RMR = {resid[valid & pre & (dist10 <= 0)].mean():+.0f} cal (passive)")
    print(f"  On-tirz below SP:  TDEE-RMR = {resid[valid & on & (dist10 > 0)].mean():+.0f} cal (suppressed)")
    print(f"  On-tirz above SP:  TDEE-RMR = {resid[valid & on & (dist10 <= 0)].mean():+.0f} cal (unchanged)")
    print(f"\n  Tirz suppresses {abs(beta[3]) / abs(beta[1]) * 100:.0f}% of the per-lb TDEE defense")
    print(f"  Net TDEE defense per lb below SP: pre-tirz {beta[1]:+.1f}, on-tirz {beta[1] + beta[3]:+.1f}")


if __name__ == "__main__":
    main()
