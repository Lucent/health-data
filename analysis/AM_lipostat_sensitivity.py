#!/usr/bin/env python3
"""AM. Lipostat sensitivity analysis — tighten confidence intervals.

Sweeps every free parameter in the set point model (AG, AH, AI) to find
the definitions and windows that maximize signal and report bootstrap CIs.

Parameters swept:
  1. Binge threshold: relative (TDEE + X) and absolute, continuous sweep
  2. Binge rate window: 30-180 days
  3. Set point half-life: fine 1-day grid
  4. Asymmetric half-lives: fine grid with bootstrap CIs
  5. Restriction definition: threshold and minimum duration
  6. Expenditure arm half-life

All key results report 95% block bootstrap confidence intervals (90-day
blocks) to respect autocorrelation.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent.parent
N_BOOT = 1000
BLOCK_SIZE = 90  # days per bootstrap block


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


def block_bootstrap_corr(x, y, n_boot=N_BOOT, block_size=BLOCK_SIZE):
    indices = block_bootstrap_indices(len(x), n_boot, block_size)
    rs = np.array([np.corrcoef(x[idx], y[idx])[0, 1] for idx in indices])
    return np.percentile(rs, [2.5, 97.5])


def partial_corr(dist, outcome, fm):
    X = np.column_stack([fm, np.ones(len(fm))])
    res_d = dist - X @ np.linalg.lstsq(X, dist, rcond=None)[0]
    res_o = outcome - X @ np.linalg.lstsq(X, outcome, rcond=None)[0]
    return np.corrcoef(res_d, res_o)[0, 1]


def block_bootstrap_partial_corr(dist, outcome, fm, n_boot=N_BOOT, block_size=BLOCK_SIZE):
    indices = block_bootstrap_indices(len(dist), n_boot, block_size)
    rs = np.array([partial_corr(dist[idx], outcome[idx], fm[idx]) for idx in indices])
    return np.percentile(rs, [2.5, 97.5])


def fit_sigmoid(centers, rates):
    from scipy.optimize import curve_fit
    def sigmoid(x, a, k, d0, base):
        return a / (1 + np.exp(-k * (x - d0))) + base

    best_result = None
    best_val = np.inf
    for a0 in [0.05, 0.10, 0.15]:
        for k0 in [0.2, 0.5, 1.0]:
            for d00 in [-5, 0, 5]:
                for b0 in [0.01, 0.03]:
                    try:
                        popt, _ = curve_fit(sigmoid, centers, rates,
                                            p0=[a0, k0, d00, b0],
                                            bounds=([0.005, 0.01, -20, 0],
                                                    [0.5, 5.0, 20, 0.15]),
                                            maxfev=5000)
                        pred = sigmoid(centers, *popt)
                        err = np.sum((pred - rates) ** 2)
                        if err < best_val:
                            best_val = err
                            best_result = popt
                    except (RuntimeError, ValueError):
                        continue
    if best_result is None:
        return (0.1, 0.3, 0, 0.02)  # fallback
    return tuple(best_result)


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
    n = len(df)
    pre = df[df["effective_level"] == 0].copy().reset_index(drop=True)
    print(f"Total days: {n}, pre-tirzepatide: {len(pre)}")

    # ═══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("1. OVERSHOOT THRESHOLD SWEEP")
    print("   Reframing 'binge' as any day above TDEE + threshold")
    print("=" * 70)

    print(f"\n  --- Relative threshold (cal above TDEE) ---")
    print(f"  {'Threshold':>10} {'N days':>7} {'Rate':>6} {'r(sym50)':>9}")

    best_rel_r, best_rel_thresh = 0, 0
    for thresh in range(100, 2200, 50):
        overshoot = (pre["calories"] > pre["tdee"] + thresh).astype(float)
        overshoot_rate = overshoot.rolling(90, min_periods=90).mean()
        sp = ema(pre["fat_mass_lbs"], 50)
        dist = sp - pre["fat_mass_lbs"]
        valid = dist.notna() & overshoot_rate.notna()
        if valid.sum() < 300:
            continue
        r = np.corrcoef(dist[valid], overshoot_rate[valid])[0, 1]
        nb = int(overshoot.sum())
        rate = overshoot.mean()
        if abs(r) > abs(best_rel_r):
            best_rel_r, best_rel_thresh = r, thresh
        if thresh % 200 == 0 or thresh == 1000:
            print(f"  {thresh:>8} cal {nb:>7} {rate:5.1%} {r:+9.3f}")

    print(f"\n  Best relative threshold: {best_rel_thresh} cal (r = {best_rel_r:+.3f})")

    # Absolute threshold
    print(f"\n  --- Absolute threshold (total cal/day) ---")
    best_abs_r, best_abs_thresh = 0, 0
    for thresh in range(1500, 4500, 100):
        overshoot = (pre["calories"] > thresh).astype(float)
        overshoot_rate = overshoot.rolling(90, min_periods=90).mean()
        sp = ema(pre["fat_mass_lbs"], 50)
        dist = sp - pre["fat_mass_lbs"]
        valid = dist.notna() & overshoot_rate.notna()
        if valid.sum() < 300:
            continue
        r = np.corrcoef(dist[valid], overshoot_rate[valid])[0, 1]
        if abs(r) > abs(best_abs_r):
            best_abs_r, best_abs_thresh = r, thresh
    print(f"  Best absolute threshold: {best_abs_thresh} cal (r = {best_abs_r:+.3f})")
    print(f"  Relative wins: {abs(best_rel_r) > abs(best_abs_r)}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("2. RATE WINDOW SWEEP")
    print("=" * 70)

    # Sweep at multiple thresholds to check interaction
    print(f"\n  {'Window':>8} {'r@'+str(best_rel_thresh):>8} {'r@500':>8} {'r@1000':>8}")
    best_win_r, best_win = 0, 90  # default to 90
    for win in [30, 45, 60, 75, 90, 105, 120, 150, 180]:
        row = f"  {win:>6}d"
        for thresh in [best_rel_thresh, 500, 1000]:
            overshoot = (pre["calories"] > pre["tdee"] + thresh).astype(float)
            overshoot_rate = overshoot.rolling(win, min_periods=win).mean()
            sp = ema(pre["fat_mass_lbs"], 50)
            dist = sp - pre["fat_mass_lbs"]
            valid = dist.notna() & overshoot_rate.notna()
            if valid.sum() < 300:
                row += f" {'n/a':>8}"
                continue
            r = np.corrcoef(dist[valid], overshoot_rate[valid])[0, 1]
            row += f" {r:+8.3f}"
            if thresh == best_rel_thresh and abs(r) > abs(best_win_r):
                best_win_r, best_win = r, win
        print(row)

    print(f"\n  Best window at {best_rel_thresh}-cal threshold: {best_win} days (r = {best_win_r:+.3f})")

    # Also test: does the CONTINUOUS surplus correlate better than any threshold?
    print(f"\n  --- Continuous surplus (no threshold) ---")
    for win in [30, 60, 90, 120]:
        mean_surplus = pre["surplus"].rolling(win, min_periods=win).mean()
        sp = ema(pre["fat_mass_lbs"], 50)
        dist = sp - pre["fat_mass_lbs"]
        valid = dist.notna() & mean_surplus.notna()
        if valid.sum() < 300:
            continue
        r = np.corrcoef(dist[valid], mean_surplus[valid])[0, 1]
        print(f"  {win:>4}d mean surplus: r = {r:+.3f}")

    # Use best settings going forward
    THRESH = best_rel_thresh
    WIN = best_win
    print(f"\n  Using: surplus > {THRESH} cal, {WIN}-day window")

    pre["overshoot"] = (pre["calories"] > pre["tdee"] + THRESH).astype(float)
    pre["overshoot_rate"] = pre["overshoot"].rolling(WIN, min_periods=WIN).mean()

    # Also compute at original 1000-cal threshold for sigmoid (where shape is meaningful)
    pre["binge"] = (pre["calories"] > pre["tdee"] + 1000).astype(float)
    pre["binge_rate_90d"] = pre["binge"].rolling(90, min_periods=90).mean()

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("3. SYMMETRIC HALF-LIFE — fine sweep with bootstrap CIs")
    print("=" * 70)

    fine_hl = []
    for hl in range(15, 400):
        sp = ema(pre["fat_mass_lbs"], hl)
        dist = sp - pre["fat_mass_lbs"]
        valid = dist.notna() & pre["overshoot_rate"].notna()
        if valid.sum() < 300:
            continue
        r = np.corrcoef(dist[valid], pre.loc[valid, "overshoot_rate"])[0, 1]
        fine_hl.append((hl, r))

    fine_hl.sort(key=lambda x: -abs(x[1]))
    best_sym_hl = fine_hl[0][0]
    best_sym_r = fine_hl[0][1]
    print(f"\n  Top 5 half-lives (1-day resolution):")
    for hl, r in fine_hl[:5]:
        print(f"    {hl:>4}d: r = {r:+.4f}")

    # Also check at original 1000-cal threshold
    fine_hl_1k = []
    for hl in range(15, 400):
        sp = ema(pre["fat_mass_lbs"], hl)
        dist = sp - pre["fat_mass_lbs"]
        valid = dist.notna() & pre["binge_rate_90d"].notna()
        if valid.sum() < 300:
            continue
        r = np.corrcoef(dist[valid], pre.loc[valid, "binge_rate_90d"])[0, 1]
        fine_hl_1k.append((hl, r))
    fine_hl_1k.sort(key=lambda x: -abs(x[1]))
    print(f"\n  At 1000-cal threshold (original AG):")
    for hl, r in fine_hl_1k[:3]:
        print(f"    {hl:>4}d: r = {r:+.4f}")

    sp = ema(pre["fat_mass_lbs"], best_sym_hl)
    dist = (sp - pre["fat_mass_lbs"]).values
    br = pre["overshoot_rate"].values
    fm = pre["fat_mass_lbs"].values
    valid = ~np.isnan(dist) & ~np.isnan(br)
    ci = block_bootstrap_corr(dist[valid], br[valid])
    print(f"\n  Best symmetric HL: {best_sym_hl}d, r = {best_sym_r:+.4f} [{ci[0]:+.4f}, {ci[1]:+.4f}]")

    valid_p = valid & ~np.isnan(fm)
    r_partial = partial_corr(dist[valid_p], br[valid_p], fm[valid_p])
    ci_p = block_bootstrap_partial_corr(dist[valid_p], br[valid_p], fm[valid_p])
    print(f"  Partial r (|FM): {r_partial:+.4f} [{ci_p[0]:+.4f}, {ci_p[1]:+.4f}]")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("4. ASYMMETRIC HALF-LIVES — fine sweep with bootstrap CIs")
    print("=" * 70)

    best_ar, best_au, best_ad = 1, 50, 50  # init to worst possible
    for hl_up in range(30, 300, 10):
        for hl_down in range(5, 100, 10):
            sp_a = asymmetric_ema(pre["fat_mass_lbs"].values, hl_up, hl_down)
            dist_a = sp_a - pre["fat_mass_lbs"].values
            valid_a = ~np.isnan(dist_a) & ~np.isnan(br)
            if valid_a.sum() < 300:
                continue
            r = np.corrcoef(dist_a[valid_a], br[valid_a])[0, 1]
            if r < best_ar:
                best_ar, best_au, best_ad = r, hl_up, hl_down

    # Fine-tune around optimum (1-day)
    for hl_up in range(max(10, best_au - 20), best_au + 21):
        for hl_down in range(max(3, best_ad - 15), best_ad + 16):
            sp_a = asymmetric_ema(pre["fat_mass_lbs"].values, hl_up, hl_down)
            dist_a = sp_a - pre["fat_mass_lbs"].values
            valid_a = ~np.isnan(dist_a) & ~np.isnan(br)
            if valid_a.sum() < 300:
                continue
            r = np.corrcoef(dist_a[valid_a], br[valid_a])[0, 1]
            if r < best_ar:
                best_ar, best_au, best_ad = r, hl_up, hl_down

    print(f"\n  Best asymmetric: HL_up={best_au}d, HL_down={best_ad}d, r = {best_ar:+.4f}")
    print(f"  Ratio (up/down): {best_au / best_ad:.1f}x")
    print(f"  87% adaptation: down={best_ad * 3}d ({best_ad * 3 / 30:.1f}mo), up={best_au * 3}d ({best_au * 3 / 30:.1f}mo)")

    sp_a = asymmetric_ema(pre["fat_mass_lbs"].values, best_au, best_ad)
    dist_a = sp_a - pre["fat_mass_lbs"].values
    valid_a = ~np.isnan(dist_a) & ~np.isnan(br)
    ci_a = block_bootstrap_corr(dist_a[valid_a], br[valid_a])
    print(f"  r = {best_ar:+.4f} [{ci_a[0]:+.4f}, {ci_a[1]:+.4f}]")

    valid_ap = valid_a & ~np.isnan(fm)
    r_partial_a = partial_corr(dist_a[valid_ap], br[valid_ap], fm[valid_ap])
    ci_ap = block_bootstrap_partial_corr(dist_a[valid_ap], br[valid_ap], fm[valid_ap])
    print(f"  Partial r (|FM): {r_partial_a:+.4f} [{ci_ap[0]:+.4f}, {ci_ap[1]:+.4f}]")

    # Improvement over symmetric
    print(f"  Δr vs symmetric: {abs(best_ar) - abs(best_sym_r):.4f}")

    # Bootstrap CI for the half-lives themselves (coarse grid, 500 boots)
    print(f"\n  --- Bootstrap CI for optimal half-lives ---")
    hl_up_boots, hl_down_boots = [], []
    n_pre = len(pre)
    boot_idx = block_bootstrap_indices(n_pre, 500, BLOCK_SIZE)
    for idx in boot_idx:
        boot_fm = pre["fat_mass_lbs"].values[idx]
        boot_br = br[idx]
        valid_b = ~np.isnan(boot_br)
        if valid_b.sum() < 200:
            continue
        b_best_r, b_best_u, b_best_d = 1, 50, 50
        for hl_up in range(30, 250, 20):
            for hl_down in range(5, 80, 15):
                sp_b = asymmetric_ema(boot_fm, hl_up, hl_down)
                dist_b = sp_b - boot_fm
                vb = ~np.isnan(dist_b) & valid_b
                if vb.sum() < 100:
                    continue
                r = np.corrcoef(dist_b[vb], boot_br[vb])[0, 1]
                if r < b_best_r:
                    b_best_r, b_best_u, b_best_d = r, hl_up, hl_down
        hl_up_boots.append(b_best_u)
        hl_down_boots.append(b_best_d)

    hl_up_boots = np.array(hl_up_boots)
    hl_down_boots = np.array(hl_down_boots)
    ratio_boots = hl_up_boots / hl_down_boots

    print(f"  HL_up:   {best_au}d [{np.percentile(hl_up_boots, 2.5):.0f}, {np.percentile(hl_up_boots, 97.5):.0f}]")
    print(f"  HL_down: {best_ad}d [{np.percentile(hl_down_boots, 2.5):.0f}, {np.percentile(hl_down_boots, 97.5):.0f}]")
    print(f"  Ratio:   {best_au / best_ad:.1f}x [{np.percentile(ratio_boots, 2.5):.1f}, {np.percentile(ratio_boots, 97.5):.1f}]")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("5. EXPENDITURE ARM — half-life sweep with bootstrap CIs")
    print("=" * 70)

    resid = pre["tdee_resid"].values
    exp_results = []
    for hl in list(range(5, 50, 1)) + list(range(50, 400, 5)):
        sp = ema(pre["fat_mass_lbs"], hl)
        dist_e = (sp - pre["fat_mass_lbs"]).values
        valid_e = ~np.isnan(dist_e) & ~np.isnan(resid) & ~np.isnan(fm)
        if valid_e.sum() < 200:
            continue
        r_pe = partial_corr(dist_e[valid_e], resid[valid_e], fm[valid_e])
        exp_results.append((hl, r_pe))

    exp_results.sort(key=lambda x: -abs(x[1]))
    best_exp_hl = exp_results[0][0]
    best_exp_r = exp_results[0][1]
    print(f"\n  Top 5 expenditure half-lives:")
    for hl, r in exp_results[:5]:
        print(f"    {hl:>4}d: partial r = {r:+.4f}")

    sp = ema(pre["fat_mass_lbs"], best_exp_hl)
    dist_e = (sp - pre["fat_mass_lbs"]).values
    valid_e = ~np.isnan(dist_e) & ~np.isnan(resid) & ~np.isnan(fm)
    ci_e = block_bootstrap_partial_corr(dist_e[valid_e], resid[valid_e], fm[valid_e])
    print(f"\n  Expenditure arm: HL={best_exp_hl}d, partial r = {best_exp_r:+.4f} [{ci_e[0]:+.4f}, {ci_e[1]:+.4f}]")
    print(f"  Appetite arm:    HL={best_sym_hl}d, partial r = {r_partial:+.4f} [{ci_p[0]:+.4f}, {ci_p[1]:+.4f}]")
    ratio_arms = abs(r_partial) / abs(best_exp_r)
    print(f"  Arm strength ratio: {ratio_arms:.1f}x (appetite / expenditure)")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("6. RESTRICTION RUNS — sweep threshold and minimum duration")
    print("=" * 70)

    sp_best = asymmetric_ema(pre["fat_mass_lbs"].values, best_au, best_ad)

    def run_analysis(deficit_thresh, min_days, post_days=30):
        restricting = pre["calories"] < (pre["tdee"] - deficit_thresh)
        runs = []
        in_run, start = False, 0
        for i in range(len(pre)):
            if restricting.iloc[i]:
                if not in_run:
                    in_run, start = True, i
            else:
                if in_run:
                    if i - start >= min_days:
                        runs.append((start, i - 1))
                    in_run = False
        if in_run and len(pre) - start >= min_days:
            runs.append((start, len(pre) - 1))

        records = []
        for s, e in runs:
            if np.isnan(sp_best[e]):
                continue
            dist_end = sp_best[e] - pre["fat_mass_lbs"].iloc[e]
            post = pre[(pre.index > e) & (pre.index <= e + post_days + 2)]
            fm_post = post["fat_mass_lbs"].iloc[-1] if len(post) > post_days - 5 else np.nan
            records.append({
                "sp_dist_end": dist_end,
                "rebound": fm_post - pre["fat_mass_lbs"].iloc[e] if not np.isnan(fm_post) else np.nan,
            })
        if len(records) < 10:
            return len(runs), np.nan
        rdf = pd.DataFrame(records).dropna(subset=["rebound"])
        if len(rdf) < 10:
            return len(runs), np.nan
        return len(runs), np.corrcoef(rdf["sp_dist_end"], rdf["rebound"])[0, 1]

    print(f"\n  {'Deficit':>8} {'MinDays':>8} {'N runs':>7} {'r(dist→rebound)':>16}")
    best_run_r, best_def, best_min = 0, 200, 3
    for deficit in [100, 150, 200, 300, 400, 500]:
        for min_d in [2, 3, 4, 5, 7]:
            n_runs, r = run_analysis(deficit, min_d)
            if not np.isnan(r) and abs(r) > abs(best_run_r):
                best_run_r, best_def, best_min = r, deficit, min_d
            print(f"  {deficit:>6} cal {min_d:>6}d {n_runs:>7} {r:+16.3f}" if not np.isnan(r) else
                  f"  {deficit:>6} cal {min_d:>6}d {n_runs:>7} {'n/a':>16}")

    print(f"\n  Best: deficit > {best_def} cal, min {best_min} days → r = {best_run_r:+.3f}")

    # Follow-up window
    print(f"\n  --- Follow-up window sweep (deficit={best_def}, min={best_min}d) ---")
    print(f"  {'Post days':>10} {'N runs':>7} {'r':>8}")
    for post_d in [7, 14, 21, 30, 45, 60, 90]:
        nr, r = run_analysis(best_def, best_min, post_d)
        print(f"  {post_d:>8}d {nr:>7} {r:+8.3f}" if not np.isnan(r) else
              f"  {post_d:>8}d {nr:>7} {'n/a':>8}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("7. DUAL DEFENSE STRENGTH — bootstrap CIs for arm ratio")
    print("=" * 70)

    sp_opt = asymmetric_ema(pre["fat_mass_lbs"].values, best_au, best_ad)
    dist_opt = sp_opt - pre["fat_mass_lbs"].values
    valid_both = (~np.isnan(dist_opt) & ~np.isnan(br) &
                  ~np.isnan(resid) & ~np.isnan(fm))
    d_v = dist_opt[valid_both]
    br_v = br[valid_both]
    resid_v = resid[valid_both]
    fm_v = fm[valid_both]

    r_appetite = partial_corr(d_v, br_v, fm_v)
    r_expenditure = partial_corr(d_v, resid_v, fm_v)

    n_v = len(d_v)
    boot_idx_dual = block_bootstrap_indices(n_v, N_BOOT, BLOCK_SIZE)
    app_boots, exp_boots, ratio_arm_boots = [], [], []
    for idx in boot_idx_dual:
        ra = partial_corr(d_v[idx], br_v[idx], fm_v[idx])
        re = partial_corr(d_v[idx], resid_v[idx], fm_v[idx])
        app_boots.append(ra)
        exp_boots.append(re)
        if abs(re) > 0.01:
            ratio_arm_boots.append(abs(ra) / abs(re))

    app_boots = np.array(app_boots)
    exp_boots = np.array(exp_boots)
    ratio_arm_boots = np.array(ratio_arm_boots)

    print(f"\n  Appetite arm partial r:    {r_appetite:+.4f} [{np.percentile(app_boots, 2.5):+.4f}, {np.percentile(app_boots, 97.5):+.4f}]")
    print(f"  Expenditure arm partial r: {r_expenditure:+.4f} [{np.percentile(exp_boots, 2.5):+.4f}, {np.percentile(exp_boots, 97.5):+.4f}]")
    print(f"  Arm ratio (|app|/|exp|):   {abs(r_appetite) / abs(r_expenditure):.1f}x [{np.percentile(ratio_arm_boots, 2.5):.1f}, {np.percentile(ratio_arm_boots, 97.5):.1f}]")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("8. SIGMOID RESPONSE CURVE at 1000-cal threshold")
    print("   (shape is meaningful only at thresholds where baseline < 10%)")
    print("=" * 70)

    # Use asymmetric SP for sigmoid
    binge_1k = pre["binge"].values
    valid_s = ~np.isnan(dist_opt) & ~np.isnan(binge_1k)
    sub_dist = dist_opt[valid_s]
    sub_binge = binge_1k[valid_s]

    bins = np.arange(-25, 25, 2.5)
    bin_centers, bin_rates, bin_ns = [], [], []
    print(f"\n  {'Distance':>20} {'Overshoot %':>12} {'n':>5}")
    for i in range(len(bins) - 1):
        mask = (sub_dist > bins[i]) & (sub_dist <= bins[i + 1])
        if mask.sum() < 20:
            continue
        center = (bins[i] + bins[i + 1]) / 2
        rate = sub_binge[mask].mean()
        bin_centers.append(center)
        bin_rates.append(rate)
        bin_ns.append(int(mask.sum()))
        print(f"  {bins[i]:+5.1f} to {bins[i+1]:+5.1f} {rate * 100:11.1f}% {mask.sum():5d}")

    bin_centers = np.array(bin_centers)
    bin_rates = np.array(bin_rates)

    a, k, d0, base = fit_sigmoid(bin_centers, bin_rates)
    print(f"\n  Sigmoid: binge% = {a * 100:.1f}% / (1+exp(-{k:.3f}*(dist-{d0:.1f}))) + {base * 100:.2f}%")
    print(f"  Baseline (at/above SP):   {base * 100:.2f}%")
    print(f"  Maximum (far below SP):   {(a + base) * 100:.1f}%")
    print(f"  Inflection:               {d0:.1f} lbs below set point")
    print(f"  Per-lb gradient at steep:  {a * k / 4 * 100:.2f}%/lb")

    # Bootstrap CIs for sigmoid (resample bins)
    print(f"\n  --- Bootstrap CIs (200 resamples) ---")
    a_b, k_b, d0_b, base_b = [], [], [], []
    for _ in range(200):
        idx = np.random.randint(0, len(bin_centers), size=len(bin_centers))
        p = fit_sigmoid(bin_centers[idx], bin_rates[idx])
        if p is not None and p[0] > 0.001 and p[1] > 0.001:
            a_b.append(p[0])
            k_b.append(p[1])
            d0_b.append(p[2])
            base_b.append(p[3])
    a_b, k_b, d0_b, base_b = [np.array(x) for x in [a_b, k_b, d0_b, base_b]]
    grad_b = a_b * k_b / 4

    print(f"  Baseline:    {base * 100:.2f}% [{np.percentile(base_b, 2.5) * 100:.2f}, {np.percentile(base_b, 97.5) * 100:.2f}]%")
    print(f"  Maximum:     {(a + base) * 100:.1f}% [{np.percentile(a_b + base_b, 2.5) * 100:.1f}, {np.percentile(a_b + base_b, 97.5) * 100:.1f}]%")
    print(f"  Inflection:  {d0:.1f} lbs [{np.percentile(d0_b, 2.5):.1f}, {np.percentile(d0_b, 97.5):.1f}]")
    print(f"  Gradient:    {a * k / 4 * 100:.2f}%/lb [{np.percentile(grad_b, 2.5) * 100:.2f}, {np.percentile(grad_b, 97.5) * 100:.2f}]%/lb")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("9. NON-BINGE CONTINUOUS PRESSURE — per-lb intake shift")
    print("=" * 70)

    valid_nb = (~np.isnan(dist_opt) & pre["surplus"].notna().values &
                (pre["binge"].values == 0))
    d_nb = dist_opt[valid_nb]
    s_nb = pre["surplus"].values[valid_nb]

    slope = np.polyfit(d_nb, s_nb, 1)[0]
    n_nb = len(d_nb)
    boot_idx_nb = block_bootstrap_indices(n_nb, N_BOOT, BLOCK_SIZE)
    slope_boots = np.array([np.polyfit(d_nb[idx], s_nb[idx], 1)[0] for idx in boot_idx_nb])
    print(f"\n  Non-binge surplus vs SP distance (asymmetric SP):")
    print(f"  Slope: {slope:+.1f} cal/day per lb [{np.percentile(slope_boots, 2.5):+.1f}, {np.percentile(slope_boots, 97.5):+.1f}]")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("SUMMARY — optimized parameters with 95% CIs")
    print("=" * 70)
    print(f"\n  Overshoot threshold: surplus > {THRESH} cal (optimal from sweep)")
    print(f"  Overshoot rate window: {WIN} days")
    print(f"  Original 1000-cal threshold retains sigmoid shape")
    print(f"")
    print(f"  Symmetric half-life:  {best_sym_hl}d, r = {best_sym_r:+.4f} [{ci[0]:+.4f}, {ci[1]:+.4f}]")
    print(f"  Asymmetric HL_up:     {best_au}d [{np.percentile(hl_up_boots, 2.5):.0f}, {np.percentile(hl_up_boots, 97.5):.0f}]")
    print(f"  Asymmetric HL_down:   {best_ad}d [{np.percentile(hl_down_boots, 2.5):.0f}, {np.percentile(hl_down_boots, 97.5):.0f}]")
    print(f"  Ratchet ratio:        {best_au / best_ad:.1f}x [{np.percentile(ratio_boots, 2.5):.1f}, {np.percentile(ratio_boots, 97.5):.1f}]")
    print(f"  Asymmetric r:         {best_ar:+.4f} [{ci_a[0]:+.4f}, {ci_a[1]:+.4f}]")
    print(f"  Partial r (|FM):      {r_partial_a:+.4f} [{ci_ap[0]:+.4f}, {ci_ap[1]:+.4f}]")
    print(f"  Expenditure arm:      HL={best_exp_hl}d, partial r = {best_exp_r:+.4f} [{ci_e[0]:+.4f}, {ci_e[1]:+.4f}]")
    print(f"  Appetite arm:         partial r = {r_appetite:+.4f} [{np.percentile(app_boots, 2.5):+.4f}, {np.percentile(app_boots, 97.5):+.4f}]")
    print(f"  Arm ratio:            {abs(r_appetite) / abs(r_expenditure):.1f}x [{np.percentile(ratio_arm_boots, 2.5):.1f}, {np.percentile(ratio_arm_boots, 97.5):.1f}]")
    print(f"  Sigmoid baseline:     {base * 100:.2f}%")
    print(f"  Sigmoid gradient:     {a * k / 4 * 100:.2f}%/lb")
    print(f"  Non-binge pressure:   {slope:+.1f} cal/lb [{np.percentile(slope_boots, 2.5):+.1f}, {np.percentile(slope_boots, 97.5):+.1f}]")
    print(f"  Restriction rebound:  r = {best_run_r:+.3f} (deficit > {best_def}, min {best_min}d)")


if __name__ == "__main__":
    main()
