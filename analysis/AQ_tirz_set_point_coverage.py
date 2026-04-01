#!/usr/bin/env python3
"""AQ. Tirzepatide set point coverage — how much gap does the drug offset?

Quantify the drug's interaction with the set point model:
1. Does tirz change the SP half-life? (AG said 165d on-drug vs 50d off)
2. How does the drug shift the surplus-vs-distance curve?
   - Intercept shift (constant cal/day suppression)?
   - Slope change (attenuates the per-lb pressure)?
   - Both?
3. At current drug levels, how many lbs of SP distance does it offset?
4. Separate appetite (mean surplus) and expenditure (TDEE residual) arms.
5. How does tachyphylaxis erode coverage over time?
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
    df = df.merge(tirz[["date", "effective_level", "blood_level", "dose_mg",
                         "days_since_injection"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df["blood_level"] = df["blood_level"].fillna(0)
    df = df.sort_values("date").reset_index(drop=True)
    df["surplus"] = df["calories"] - df["tdee"]
    df["tdee_resid"] = df["tdee"] - df["expected_rmr"]
    df["on_tirz"] = (df["effective_level"] > 0).astype(int)

    fm = df["fat_mass_lbs"].values
    surplus = df["surplus"].values
    resid = df["tdee_resid"].values

    # SP using asymmetric model (AN parameters)
    sp = asymmetric_ema(fm, 72, 25)
    dist = sp - fm  # positive = below SP

    # 90-day mean surplus (AP best metric)
    mean_surplus_90 = pd.Series(surplus).rolling(90, min_periods=90).mean().values

    pre_mask = df["on_tirz"].values == 0
    on_mask = df["on_tirz"].values == 1

    # ═══════════════════════════════════════════════════════════════════
    print("=" * 70)
    print("1. SP HALF-LIFE ON vs OFF DRUG — using mean surplus metric")
    print("=" * 70)

    # Sweep symmetric HL for pre-tirz and on-tirz separately
    print(f"\n  {'HL':>6} {'r (pre-tirz)':>13} {'r (on-tirz)':>13}")
    for hl in [20, 30, 40, 50, 60, 80, 100, 120, 150, 200, 300]:
        sp_test = ema(df["fat_mass_lbs"], hl).values
        dist_test = sp_test - fm

        for label, mask in [("pre", pre_mask), ("on", on_mask)]:
            rolled = mean_surplus_90.copy()
            valid = ~np.isnan(dist_test) & ~np.isnan(rolled) & mask
            if valid.sum() < 100:
                continue

        valid_pre = ~np.isnan(dist_test) & ~np.isnan(mean_surplus_90) & pre_mask
        valid_on = ~np.isnan(dist_test) & ~np.isnan(mean_surplus_90) & on_mask

        r_pre = np.corrcoef(dist_test[valid_pre], mean_surplus_90[valid_pre])[0, 1] if valid_pre.sum() > 100 else np.nan
        r_on = np.corrcoef(dist_test[valid_on], mean_surplus_90[valid_on])[0, 1] if valid_on.sum() > 50 else np.nan

        print(f"  {hl:>4}d {r_pre:+13.4f} {r_on:+13.4f}" if not np.isnan(r_on) else
              f"  {hl:>4}d {r_pre:+13.4f} {'n/a':>13}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("2. SURPLUS vs SP DISTANCE — separate curves for pre/on-tirz")
    print("=" * 70)

    valid = ~np.isnan(dist) & ~np.isnan(surplus)

    # Bin by distance, show mean daily surplus
    bins = [(-10, -5), (-5, -2.5), (-2.5, 0), (0, 2.5), (2.5, 5), (5, 10), (10, 20)]
    print(f"\n  {'Distance':>15} {'Pre-tirz':>12} {'On-tirz':>12} {'Drug Δ':>10}")
    for lo, hi in bins:
        mask_pre = valid & pre_mask & (dist > lo) & (dist <= hi)
        mask_on = valid & on_mask & (dist > lo) & (dist <= hi)
        n_pre = mask_pre.sum()
        n_on = mask_on.sum()
        if n_pre < 30:
            continue
        s_pre = surplus[mask_pre].mean()
        s_on = surplus[mask_on].mean() if n_on >= 10 else np.nan
        delta = s_on - s_pre if not np.isnan(s_on) else np.nan
        on_str = f"{s_on:+12.0f}" if not np.isnan(s_on) else f"{'n/a (n=' + str(n_on) + ')':>12}"
        d_str = f"{delta:+10.0f}" if not np.isnan(delta) else ""
        print(f"  {lo:+.0f} to {hi:+.0f} lbs {s_pre:+12.0f} {on_str} {d_str}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("3. REGRESSION — intercept shift vs slope change")
    print("=" * 70)

    # Use 90-day smoothed surplus for regressions (matches AP finding)
    smoothed_surplus = mean_surplus_90
    valid_reg = ~np.isnan(dist) & ~np.isnan(smoothed_surplus) & ~np.isnan(fm)
    d_v = dist[valid_reg]
    s_v = smoothed_surplus[valid_reg]
    f_v = fm[valid_reg]
    on_v = df["on_tirz"].values[valid_reg].astype(float)
    eff_v = df["effective_level"].values[valid_reg]

    # Model A: surplus ~ dist + FM + on_tirz (intercept shift only)
    XA = np.column_stack([d_v, f_v, on_v, np.ones(valid_reg.sum())])
    betaA = np.linalg.lstsq(XA, s_v, rcond=None)[0]
    residA = s_v - XA @ betaA
    r2A = 1 - np.sum(residA**2) / np.sum((s_v - s_v.mean())**2)

    # Model B: surplus ~ dist + FM + on_tirz + dist×on_tirz (slope change)
    XB = np.column_stack([d_v, f_v, on_v, d_v * on_v, np.ones(valid_reg.sum())])
    betaB = np.linalg.lstsq(XB, s_v, rcond=None)[0]
    residB = s_v - XB @ betaB
    r2B = 1 - np.sum(residB**2) / np.sum((s_v - s_v.mean())**2)

    # Model C: surplus ~ dist + FM + effective_level (continuous drug level)
    XC = np.column_stack([d_v, f_v, eff_v, np.ones(valid_reg.sum())])
    betaC = np.linalg.lstsq(XC, s_v, rcond=None)[0]
    residC = s_v - XC @ betaC
    r2C = 1 - np.sum(residC**2) / np.sum((s_v - s_v.mean())**2)

    # Model D: surplus ~ dist + FM + eff_level + dist×eff_level
    XD = np.column_stack([d_v, f_v, eff_v, d_v * eff_v, np.ones(valid_reg.sum())])
    betaD = np.linalg.lstsq(XD, s_v, rcond=None)[0]
    residD = s_v - XD @ betaD
    r2D = 1 - np.sum(residD**2) / np.sum((s_v - s_v.mean())**2)

    print(f"\n  Model A (intercept shift): R² = {r2A:.4f}")
    print(f"    dist: {betaA[0]:+.1f} cal/lb, FM: {betaA[1]:+.1f}, on_tirz: {betaA[2]:+.0f} cal")
    print(f"\n  Model B (slope change): R² = {r2B:.4f}")
    print(f"    dist: {betaB[0]:+.1f} cal/lb, on_tirz: {betaB[2]:+.0f}, dist×tirz: {betaB[3]:+.1f}")
    print(f"\n  Model C (continuous drug level): R² = {r2C:.4f}")
    print(f"    dist: {betaC[0]:+.1f} cal/lb, FM: {betaC[1]:+.1f}, eff_level: {betaC[2]:+.1f} cal/unit")
    print(f"\n  Model D (continuous + interaction): R² = {r2D:.4f}")
    print(f"    dist: {betaD[0]:+.1f} cal/lb, eff_level: {betaD[2]:+.1f}, dist×eff: {betaD[3]:+.1f}")

    # Lbs-equivalent: how many lbs of SP distance does 1 unit of eff_level offset?
    if betaC[0] != 0:
        lbs_per_unit = -betaC[2] / betaC[0]
        print(f"\n  Drug equivalence (Model C): 1 unit effective level = {lbs_per_unit:.1f} lbs of SP offset")

    # Bootstrap CI for drug equivalence
    boot_equiv = []
    boot_idx = block_bootstrap_indices(valid_reg.sum(), N_BOOT, BLOCK_SIZE)
    for idx in boot_idx:
        try:
            b = np.linalg.lstsq(XC[idx], s_v[idx], rcond=None)[0]
            if abs(b[0]) > 0.1:
                boot_equiv.append(-b[2] / b[0])
        except:
            pass
    boot_equiv = np.array(boot_equiv)
    boot_equiv = boot_equiv[(boot_equiv > -50) & (boot_equiv < 50)]
    print(f"  Bootstrap 95% CI: [{np.percentile(boot_equiv, 2.5):.1f}, {np.percentile(boot_equiv, 97.5):.1f}] lbs/unit")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("4. EXPENDITURE ARM — drug effect on TDEE defense")
    print("=" * 70)

    # TDEE residual — use daily (already smooth from Kalman)
    valid_exp_full = ~np.isnan(dist) & ~np.isnan(resid) & ~np.isnan(fm)

    d_e = dist[valid_exp_full]
    r_e = resid[valid_exp_full]
    f_e = fm[valid_exp_full]
    on_e = df["on_tirz"].values[valid_exp_full].astype(float)
    eff_e = df["effective_level"].values[valid_exp_full]

    # Model: TDEE_resid ~ dist + FM + eff_level
    XE = np.column_stack([d_e, f_e, eff_e, np.ones(valid_exp_full.sum())])
    betaE = np.linalg.lstsq(XE, r_e, rcond=None)[0]
    residE = r_e - XE @ betaE
    r2E = 1 - np.sum(residE**2) / np.sum((r_e - r_e.mean())**2)

    # Model with interaction
    XF = np.column_stack([d_e, f_e, eff_e, d_e * eff_e, np.ones(valid_exp_full.sum())])
    betaF = np.linalg.lstsq(XF, r_e, rcond=None)[0]
    residF = r_e - XF @ betaF
    r2F = 1 - np.sum(residF**2) / np.sum((r_e - r_e.mean())**2)

    print(f"\n  TDEE_resid ~ dist + FM + eff_level: R² = {r2E:.4f}")
    print(f"    dist: {betaE[0]:+.1f} cal/lb, FM: {betaE[1]:+.1f}, eff_level: {betaE[2]:+.1f} cal/unit")

    print(f"\n  With interaction: R² = {r2F:.4f}")
    print(f"    dist: {betaF[0]:+.1f}, eff_level: {betaF[2]:+.1f}, dist×eff: {betaF[3]:+.1f}")

    if betaE[0] != 0:
        lbs_exp = -betaE[2] / betaE[0]
        print(f"\n  Expenditure arm: 1 unit eff_level = {lbs_exp:.1f} lbs of SP offset")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("5. COVERAGE TABLE — at current and future drug levels")
    print("=" * 70)

    # Current state
    current_fm = fm[-1]
    current_sp = sp[-1]
    current_gap = current_sp - current_fm
    current_eff = df["effective_level"].values[-1]
    current_blood = df["blood_level"].values[-1]

    print(f"\n  Current state:")
    print(f"    FM = {current_fm:.1f}, SP = {current_sp:.1f}, gap = {current_gap:+.1f} lbs")
    print(f"    Effective level = {current_eff:.1f}, blood level = {current_blood:.1f}")

    # Without drug: expected 90d mean surplus from the gap
    surplus_no_drug = betaC[0] * current_gap
    surplus_with_drug = betaC[0] * current_gap + betaC[2] * current_eff

    print(f"\n  Expected 90d mean surplus at current gap:")
    print(f"    Without drug: {surplus_no_drug:+.0f} cal/day (from {current_gap:+.1f} lb gap × {betaC[0]:+.1f} cal/lb)")
    print(f"    With drug:    {surplus_with_drug:+.0f} cal/day (drug adds {betaC[2] * current_eff:+.0f} cal)")

    # Coverage at various drug levels
    print(f"\n  Appetite arm coverage:")
    print(f"  {'Eff level':>10} {'Offset (lbs)':>13} {'Max gap covered':>16} {'Cal/day suppression':>20}")
    for eff in [2, 4, 6, 8, 10, 12]:
        offset = lbs_per_unit * eff
        suppression = betaC[2] * eff
        print(f"  {eff:>10.0f} {offset:>+13.1f} {offset:>14.1f} lbs {suppression:>+18.0f}")

    # Expenditure arm coverage
    print(f"\n  Expenditure arm coverage:")
    print(f"  {'Eff level':>10} {'TDEE suppression':>17}")
    for eff in [2, 4, 6, 8, 10, 12]:
        suppression = betaE[2] * eff
        print(f"  {eff:>10.0f} {suppression:>+15.0f} cal")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("6. TACHYPHYLAXIS EROSION — how coverage degrades over time")
    print("=" * 70)

    # Tachyphylaxis model: effective = blood × exp(-0.0217 × weeks_on_dose)
    # At 12.5mg: blood level at trough ≈ 8-10, peak ≈ 12-15
    # Tachyphylaxis HL = 32 weeks
    tachy_hl_weeks = 32
    print(f"\n  Tachyphylaxis half-life: {tachy_hl_weeks} weeks")
    print(f"\n  {'Weeks on dose':>14} {'Effectiveness':>14} {'Eff level (mean)':>17} {'Appetite offset':>16} {'TDEE suppression':>17}")
    for weeks in [0, 8, 16, 24, 32, 48, 64, 80]:
        effectiveness = np.exp(-0.0217 * weeks)
        mean_blood = 10  # approximate mean across sawtooth
        eff_level = mean_blood * effectiveness
        app_offset = lbs_per_unit * eff_level
        tdee_supp = betaE[2] * eff_level
        print(f"  {weeks:>12} wk {effectiveness:>13.0%} {eff_level:>15.1f} {app_offset:>+14.1f} lbs {tdee_supp:>+15.0f} cal")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("7. INJECTION CYCLE — within-week sawtooth")
    print("=" * 70)

    on_tirz_df = df[df["on_tirz"] == 1].copy()
    if "days_since_injection" in on_tirz_df.columns:
        print(f"\n  {'Day post-inj':>13} {'Mean cal':>9} {'Mean surplus':>13} {'Eff level':>10} {'n':>5}")
        for day in range(8):
            mask = on_tirz_df["days_since_injection"] == day
            if mask.sum() < 10:
                continue
            sub = on_tirz_df[mask]
            print(f"  {day:>13} {sub['calories'].mean():>9.0f} {sub['surplus'].mean():>+13.0f} {sub['effective_level'].mean():>10.1f} {len(sub):>5}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)
    print(f"\n  Appetite arm:")
    print(f"    Per lb below SP: {betaC[0]:+.1f} cal/day surplus")
    print(f"    Per unit eff level: {betaC[2]:+.1f} cal/day suppression")
    print(f"    Drug equivalence: 1 unit = {lbs_per_unit:.1f} lbs [{np.percentile(boot_equiv, 2.5):.1f}, {np.percentile(boot_equiv, 97.5):.1f}]")
    print(f"  Expenditure arm:")
    print(f"    Per lb below SP: {betaE[0]:+.1f} cal/day TDEE elevation")
    print(f"    Per unit eff level: {betaE[2]:+.1f} cal/day TDEE suppression")
    print(f"  Current coverage:")
    print(f"    Gap: {current_gap:+.1f} lbs, eff level: {current_eff:.1f}")
    print(f"    Appetite: gap pressure {surplus_no_drug:+.0f}, drug offset {betaC[2] * current_eff:+.0f}, net {surplus_with_drug:+.0f}")


if __name__ == "__main__":
    main()
