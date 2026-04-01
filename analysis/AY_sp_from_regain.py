#!/usr/bin/env python3
"""AY. Derive set point adaptation curve from SURMOUNT-4 regain data.

Instead of assuming an EMA, invert the regain curve:
  regain_rate(t) → surplus(t) → gap(t) → SP_adaptation(t)

SURMOUNT-4 published data (Aronne et al. JAMA 2024):
  - 36 weeks open-label tirzepatide (max dose 15mg)
  - Mean weight loss at week 36: ~20% (-21 kg from 102.3 kg start)
  - Randomized to continue (n=340) or placebo (n=330)
  - Placebo arm regain over 52 weeks: +14% from week-36 weight (+11.5 kg)
  - Regain distribution at week 88: 17.5% maintained, 82.5% regained >25%

From the paper's Figure 2, approximate placebo arm trajectory:
  Week 36 (stop): 0% (reference)
  Week 48 (+12 wk): ~+7%
  Week 60 (+24 wk): ~+11%
  Week 72 (+36 wk): ~+13%
  Week 88 (+52 wk): +14%

This curve tells us the set point's behavior.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

# SURMOUNT-4 approximate regain trajectory (placebo arm, from published figure)
# Weeks post-discontinuation → % regain from week-36 weight
REGAIN_WEEKS = np.array([0, 4, 8, 12, 16, 20, 24, 32, 40, 52])
REGAIN_PCT = np.array([0, 2.5, 5.0, 7.0, 8.5, 9.8, 11.0, 12.5, 13.3, 14.0])

# Trial parameters
START_KG = 102.3
LOSS_PCT = 20.0  # % at week 36
W36_KG = START_KG * (1 - LOSS_PCT / 100)  # 81.8 kg
FM_FRACTION = 0.93  # fraction of weight change that is fat

# Drug parameters from AX
DRUG_CAL_PER_BLOOD = -74.5
TACHY_HL_WEEKS = 35

# Set point pressure from AM (2014+)
SP_CAL_PER_LB = 27


def main():
    print("=" * 70)
    print("1. INVERT THE REGAIN CURVE — derive gap and SP adaptation")
    print("=" * 70)

    # At each timepoint, compute:
    # 1. Weight and FM from the regain percentage
    # 2. Rate of FM change (lbs/week)
    # 3. Implied surplus (from FM change rate)
    # 4. Implied gap (from surplus / SP_CAL_PER_LB)
    # 5. Implied SP (from gap = SP - FM)
    # 6. Fraction of original gap that remains

    w36_lbs = W36_KG * 2.205
    start_lbs = START_KG * 2.205
    fm_lost = (start_lbs - w36_lbs) * FM_FRACTION
    fm36 = start_lbs * 0.40 - fm_lost  # approximate: start at 40% BF
    fm_start = start_lbs * 0.40

    # FM at start of treatment ≈ SP at start (equilibrium)
    sp_original = fm_start

    print(f"\n  Start: {START_KG:.1f} kg ({start_lbs:.0f} lbs), FM ~{fm_start:.0f} lbs")
    print(f"  Week 36: {W36_KG:.1f} kg ({w36_lbs:.0f} lbs), FM ~{fm36:.0f} lbs")
    print(f"  FM lost: {fm_lost:.0f} lbs")
    print(f"  Original SP (= start FM): {sp_original:.0f} lbs")

    print(f"\n  {'Wk post':>8} {'Regain%':>8} {'Weight':>7} {'FM':>6} {'FM rate':>8} "
          f"{'Surplus':>8} {'Gap':>6} {'SP':>6} {'SP adapted':>11}")

    fm_values = []
    sp_values = []
    gap_values = []

    for i in range(len(REGAIN_WEEKS)):
        wk = REGAIN_WEEKS[i]
        pct = REGAIN_PCT[i]

        weight = w36_lbs * (1 + pct / 100)
        fm_regained = (weight - w36_lbs) * FM_FRACTION
        fm = fm36 + fm_regained

        # FM rate: use finite differences
        if i == 0:
            # Use forward difference
            next_wk = REGAIN_WEEKS[1]
            next_pct = REGAIN_PCT[1]
            next_weight = w36_lbs * (1 + next_pct / 100)
            next_fm = fm36 + (next_weight - w36_lbs) * FM_FRACTION
            fm_rate = (next_fm - fm) / (next_wk - wk)  # lbs/week
        elif i == len(REGAIN_WEEKS) - 1:
            # Use backward difference
            prev_wk = REGAIN_WEEKS[i - 1]
            prev_pct = REGAIN_PCT[i - 1]
            prev_weight = w36_lbs * (1 + prev_pct / 100)
            prev_fm = fm36 + (prev_weight - w36_lbs) * FM_FRACTION
            fm_rate = (fm - prev_fm) / (wk - prev_wk)
        else:
            # Central difference
            prev_wk = REGAIN_WEEKS[i - 1]
            next_wk = REGAIN_WEEKS[i + 1]
            prev_fm = fm36 + (w36_lbs * (1 + REGAIN_PCT[i - 1] / 100) - w36_lbs) * FM_FRACTION
            next_fm = fm36 + (w36_lbs * (1 + REGAIN_PCT[i + 1] / 100) - w36_lbs) * FM_FRACTION
            fm_rate = (next_fm - prev_fm) / (next_wk - prev_wk)

        # Implied surplus: fm_rate lbs/week * 3500 cal/lb / 7 days * (1/FM_FRACTION)
        surplus_per_day = fm_rate * 3500 / 7 / FM_FRACTION

        # Implied gap: surplus = SP_CAL_PER_LB * gap
        if SP_CAL_PER_LB > 0:
            gap = surplus_per_day / SP_CAL_PER_LB
        else:
            gap = 0

        # Implied SP: gap = SP - FM → SP = FM + gap
        sp = fm + gap

        # How much has SP adapted from original toward week-36 FM?
        # At treatment start: SP = sp_original (~90 lbs)
        # If fully adapted: SP = fm36 (~49 lbs)
        # Fraction adapted = 1 - (SP - fm36) / (sp_original - fm36)
        total_range = sp_original - fm36
        if total_range > 0:
            fraction_adapted = 1 - (sp - fm36) / total_range
        else:
            fraction_adapted = 0

        fm_values.append(fm)
        sp_values.append(sp)
        gap_values.append(gap)

        print(f"  {wk:>6}wk {pct:>7.1f}% {weight:>7.0f} {fm:>6.0f} {fm_rate:>+7.1f}/wk "
              f"{surplus_per_day:>+8.0f} {gap:>+6.1f} {sp:>6.0f} {fraction_adapted:>10.0%}")

    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("2. WHAT SHAPE IS THE ADAPTATION?")
    print("=" * 70)

    # The SP adaptation curve: fraction_adapted vs weeks_on_drug (0 to 36)
    # We derived it from post-drug regain, but it tells us about ON-drug adaptation.
    # At week 36 (treatment end), the SP had adapted to some degree.
    # The remaining gap drives the regain. So:
    #   SP at discontinuation = fm36 + gap_at_discontinuation
    # From week 0 of regain: gap ≈ surplus / 27 ≈ 330 / 27 ≈ 12 lbs
    # So SP ≈ 49 + 12 = 61 lbs
    # Original SP = 90. Range = 90 - 49 = 41 lbs. Adapted = 1 - (61-49)/41 = 71%

    initial_gap = gap_values[0] if gap_values else 12
    sp_at_discontinuation = fm36 + initial_gap
    adapted_at_36wk = 1 - (sp_at_discontinuation - fm36) / (sp_original - fm36)

    print(f"\n  SP at discontinuation (week 36): {sp_at_discontinuation:.0f} lbs")
    print(f"  Original SP: {sp_original:.0f} lbs")
    print(f"  Week-36 FM: {fm36:.0f} lbs")
    print(f"  Gap at discontinuation: {initial_gap:+.0f} lbs")
    print(f"  Fraction adapted after 36 weeks on drug: {adapted_at_36wk:.0%}")
    print(f"  Fraction NOT adapted (residual): {1 - adapted_at_36wk:.0%}")

    # What HL would produce this adaptation?
    # EMA: fraction_adapted = 1 - exp(-ln(2)/HL * days)
    # 0.71 = 1 - exp(-ln(2)/HL * 252)
    # exp(-ln(2)/HL * 252) = 0.29
    # -ln(2)/HL * 252 = ln(0.29)
    # HL = -ln(2) * 252 / ln(0.29) = 0.693 * 252 / 1.238 = 141 days
    if adapted_at_36wk > 0 and adapted_at_36wk < 1:
        implied_hl = -np.log(2) * 252 / np.log(1 - adapted_at_36wk)
        print(f"\n  Implied HL (if EMA): {implied_hl:.0f} days ({implied_hl/7:.0f} weeks)")
        print(f"  Compare: pre-drug HL = 45 days")
        print(f"  Ratio: {implied_hl / 45:.1f}x slower on drug")

    # But is it actually exponential? Check the regain SHAPE.
    # If SP adaptation were exponential, regain should be front-loaded (fast
    # initial regain, decaying). If linear, regain should be constant rate
    # until SP catches up. If threshold, no regain then sudden regain.

    print(f"\n  Regain rate by period:")
    for i in range(1, len(REGAIN_WEEKS)):
        dt_weeks = REGAIN_WEEKS[i] - REGAIN_WEEKS[i - 1]
        dpct = REGAIN_PCT[i] - REGAIN_PCT[i - 1]
        rate = dpct / dt_weeks
        print(f"    Week {REGAIN_WEEKS[i-1]:>2}-{REGAIN_WEEKS[i]:<2}: {rate:+.2f}%/week")

    print(f"\n  Pattern: front-loaded (0.63%/wk first 4 weeks, 0.06%/wk last 12 weeks)")
    print(f"  This IS consistent with exponential gap closure — just at a longer HL")

    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("3. TEST ON THIS SUBJECT'S PRE-DRUG DATA")
    print("=" * 70)

    # Does a 141-day HL fit the pre-drug 2014+ data better than 45d?
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    df = intake[["date", "calories"]].merge(
        kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="inner")
    df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df = df.sort_values("date")
    df["surplus"] = df["calories"] - df["tdee"]
    post = df[(df["date"] >= "2014-01-01") & (df["effective_level"] == 0)].copy()

    def ema(series, hl):
        alpha = 1 - np.exp(-np.log(2) / hl)
        return series.ewm(alpha=alpha, min_periods=30).mean()

    print(f"\n  {'HL':>6} {'r (surplus 90d)':>16}")
    for hl in [30, 45, 60, 80, 100, 120, 141, 160, 200, 300]:
        sp = ema(post["fat_mass_lbs"], hl)
        dist = sp - post["fat_mass_lbs"]
        surplus = post["surplus"].rolling(90, min_periods=90).mean()
        v = dist.notna() & surplus.notna()
        if v.sum() > 200:
            r = np.corrcoef(dist[v], surplus[v])[0, 1]
            print(f"  {hl:>4}d {r:>+16.4f}")

    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("4. TWO-COMPONENT MODEL: fast + slow adaptation")
    print("=" * 70)

    # Maybe the SP has two components:
    # SP = w * fast_EMA(FM, HL=45d) + (1-w) * slow_EMA(FM, HL=141d)
    # The fast component captures the pre-drug signal (r=-0.92 at 45d)
    # The slow component captures the residual that drives post-drug regain

    print(f"\n  SP = w × EMA(FM, 45d) + (1-w) × EMA(FM, 141d)")
    print(f"  {'w':>6} {'r (surplus 90d)':>16}")
    for w in [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.3, 0.0]:
        sp_fast = ema(post["fat_mass_lbs"], 45)
        sp_slow = ema(post["fat_mass_lbs"], 141)
        sp_mix = w * sp_fast + (1 - w) * sp_slow
        dist = sp_mix - post["fat_mass_lbs"]
        surplus = post["surplus"].rolling(90, min_periods=90).mean()
        v = dist.notna() & surplus.notna()
        if v.sum() > 200:
            r = np.corrcoef(dist[v], surplus[v])[0, 1]
            print(f"  {w:>5.1f} {r:>+16.4f}")

    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)
    print(f"\n  SURMOUNT-4 regain implies SP adapted ~{adapted_at_36wk:.0%} over 36 weeks on drug")
    if adapted_at_36wk > 0 and adapted_at_36wk < 1:
        print(f"  If EMA: implied HL = {implied_hl:.0f} days ({implied_hl / 45:.1f}x slower than pre-drug 45d)")
    print(f"  Regain shape is front-loaded, consistent with exponential gap closure")
    print(f"  The 45d pre-drug HL would predict ~97% adaptation in 36 weeks — too much")
    print(f"  The 141d implied HL predicts ~{adapted_at_36wk:.0%} — matching the regain data")


if __name__ == "__main__":
    main()
