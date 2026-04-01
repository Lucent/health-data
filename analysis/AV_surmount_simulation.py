#!/usr/bin/env python3
"""AV. Forward simulation of SURMOUNT-1 weight loss from set point model.

Uses the set point parameters derived entirely from this dataset (AG, AM, AQ)
to predict weight loss trajectories for SURMOUNT-1 trial participants at
5mg, 10mg, and 15mg tirzepatide, then compares against published results.

Model parameters (all from this dataset, none fitted to trial data):
  - Set point: EMA of FM, HL = 50d (AG/AM)
  - Eating pressure: -49 cal per unit effective level (F)
  - Metabolic suppression: -9 cal per unit effective level (AQ)
  - Tachyphylaxis: 32-week HL (F)
  - PK: FDA one-compartment (t½=5d, Tmax=24h, ka=3.31/day)
  - Set point eating pressure: -45 cal/lb (AM) — but this is from one subject.
    For the trial we use -49 cal/unit directly since we don't know trial
    subjects' set point dynamics.
  - Ratchet: HL_down=25d (AM/AN)

Published SURMOUNT-1 results (Jastreboff et al. NEJM 2022):
  - Start: 104.8 kg mean, BMI 38
  - 72 weeks: placebo -2.4%, 5mg -16.0%, 10mg -21.4%, 15mg -22.5%
  - Dose escalation: all start 2.5mg, increase q4w to target

SURMOUNT-4 discontinuation (Aronne et al. JAMA 2024):
  - 36 weeks on-drug: ~20% loss
  - 52 weeks post-discontinuation: +14% regain (from week-36 weight)
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

# ═══════════════════════════════════════════════════════════════════
# MODEL PARAMETERS — all from this dataset
# ═══════════════════════════════════════════════════════════════════

# PK model (FDA label)
PK_HALF_LIFE_DAYS = 5.0
PK_KA = 3.31  # absorption rate constant (1/day)
PK_TMAX_HOURS = 24

# Drug effects (AX clean within-week identification)
APPETITE_CAL_PER_BLOOD = -74.5  # cal/day per unit BLOOD level (AX within-week)
METABOLIC_CAL_PER_UNIT = 0     # within-week: zero; metabolic arm is slower
TACHYPHYLAXIS_HL_WEEKS = 35    # weeks, cumulative exposure (AX)

# Set point (AM, 2014+ natural dynamics data)
SP_HL_DOWN = 45   # days — symmetric (ratchet retracted)
SP_HL_UP = 45     # days — symmetric
SP_EATING_PRESSURE = 27  # cal/day per lb of gap (AM 2014+: -27 cal/lb)

# Body composition
FORBES_FAT_FRACTION = 0.85  # fraction of weight change that is fat at high FM
CAL_PER_LB_FAT = 3500

# Baseline metabolic rate (approximate for trial population)
# SURMOUNT-1 mean: 104.8 kg, ~38 BMI, ~67% female
# Rough TDEE for sedentary 104.8 kg person: ~2200 cal/day
BASELINE_TDEE = 2200
PLACEBO_DRIFT = -0.024 / (72 * 7) * 104.8 * 2.205  # -2.4% over 72 weeks in lbs/day


def pk_blood_level(dose_schedule, day):
    """Sum contributions from all prior injections using one-compartment SC model."""
    ke = np.log(2) / PK_HALF_LIFE_DAYS
    level = 0.0
    for inj_day, dose_mg in dose_schedule:
        if inj_day > day:
            break
        dt = day - inj_day
        if dt < 0:
            continue
        # One-compartment: C(t) = dose * ka / (ka - ke) * (exp(-ke*t) - exp(-ka*t))
        if dt == 0:
            level += 0
        else:
            level += dose_mg * PK_KA / (PK_KA - ke) * (np.exp(-ke * dt) - np.exp(-PK_KA * dt))
    return level


def tachyphylaxis_factor(weeks_on_current_dose):
    """Effectiveness decay: exp(-ln(2)/HL * weeks)."""
    return np.exp(-np.log(2) / TACHYPHYLAXIS_HL_WEEKS * weeks_on_current_dose)


def build_dose_schedule(target_mg, n_weeks):
    """SURMOUNT-1 dose escalation: start 2.5mg, increase q4w."""
    schedule = []
    escalation = {
        5: [(0, 2.5), (4, 5.0)],
        10: [(0, 2.5), (4, 5.0), (8, 7.5), (12, 10.0)],
        15: [(0, 2.5), (4, 5.0), (8, 7.5), (12, 10.0), (16, 12.5), (20, 15.0)],
    }
    dose_changes = escalation[target_mg]
    current_dose = 2.5
    for week_start, dose in dose_changes:
        current_dose = dose

    # Build weekly injection schedule
    dose_idx = 0
    current_dose = dose_changes[0][1]
    for week in range(n_weeks):
        # Check for dose increase
        if dose_idx + 1 < len(dose_changes) and week >= dose_changes[dose_idx + 1][0]:
            dose_idx += 1
            current_dose = dose_changes[dose_idx][1]
        schedule.append((week * 7, current_dose))

    return schedule


def simulate(start_weight_lbs, target_mg, n_weeks=72, start_fm_fraction=0.40,
             sp_gap_lbs=0.0, on_drug_sp_hl_down=None):
    """Simulate weight loss trajectory.

    sp_gap_lbs: starting gap between FM and SP. Negative means FM is above SP
    (subject was gaining, SP hasn't caught up — this HELPS weight loss).
    Positive means FM is below SP (subject was recently dieting — headwind).
    on_drug_sp_hl_down: if set, overrides SP_HL_DOWN while drug is active.
    """
    sp_hl_down_active = on_drug_sp_hl_down if on_drug_sp_hl_down else SP_HL_DOWN
    n_days = n_weeks * 7
    start_fm = start_weight_lbs * start_fm_fraction
    start_lean = start_weight_lbs - start_fm

    dose_schedule = build_dose_schedule(target_mg, n_weeks)

    # Track state
    fm = np.zeros(n_days)
    sp = np.zeros(n_days)
    weight = np.zeros(n_days)
    intake = np.zeros(n_days)
    tdee = np.zeros(n_days)
    eff_level = np.zeros(n_days)

    fm[0] = start_fm
    sp[0] = start_fm - sp_gap_lbs  # negative gap → SP below FM → helps loss
    weight[0] = start_weight_lbs

    # Track weeks on current dose for tachyphylaxis
    current_dose = 0
    weeks_on_dose = 0
    last_dose_change_day = 0

    for d in range(1, n_days):
        # PK: blood level from all prior injections
        blood = pk_blood_level(dose_schedule, d)

        # Tachyphylaxis: cumulative weeks on drug (not per-dose reset)
        cumulative_weeks = d / 7
        tachy = tachyphylaxis_factor(cumulative_weeks)
        eff = blood * tachy
        eff_level[d] = eff

        # Set point update (asymmetric EMA)
        # Ratchet tracks DIRECTION OF FM CHANGE, not gap direction.
        # FM falling → body accepts lower weight quickly (HL_down)
        # FM rising → body resists raising defended weight (HL_up)
        fm_direction = fm[d - 1] - fm[max(0, d - 7)]  # 7-day FM trend
        if fm_direction <= 0:
            # FM falling or stable → SP adapts down quickly
            alpha = 1 - np.exp(-np.log(2) / sp_hl_down_active)
        else:
            # FM rising → SP adapts up slowly
            alpha = 1 - np.exp(-np.log(2) / SP_HL_UP)
        sp[d] = sp[d - 1] + alpha * (fm[d - 1] - sp[d - 1])

        # Gap and eating pressure
        gap = sp[d] - fm[d - 1]  # positive = FM below SP, eating pressure up
        eating_pressure = SP_EATING_PRESSURE * gap  # cal/day

        # TDEE: baseline adjusted for body composition
        weight_ratio = weight[d - 1] / start_weight_lbs
        base_tdee = BASELINE_TDEE * (weight_ratio ** 0.75)

        # Drug effects
        appetite_effect = APPETITE_CAL_PER_BLOOD * eff  # negative = less eating
        # Metabolic arm: drug suppresses the body's calorie-burning boost
        # (zero within-week from AX, but structurally included for the simulator)
        metabolic_effect = METABOLIC_CAL_PER_UNIT * eff  # modifies TDEE independently

        # Actual TDEE = base + metabolic drug effect
        tdee[d] = base_tdee + metabolic_effect

        # Intake = base_tdee + eating pressure + drug appetite effect
        # (Intake is driven off base_tdee, NOT the drug-modified actual tdee)
        intake[d] = base_tdee + eating_pressure + appetite_effect

        # Net surplus: intake - actual TDEE (metabolic arm affects this independently)
        surplus = intake[d] - tdee[d]
        # = (base + pressure + appetite) - (base + metabolic)
        # = pressure + appetite - metabolic
        # All three terms contribute independently.

        # FM change
        fm_change = surplus * FORBES_FAT_FRACTION / CAL_PER_LB_FAT
        fm[d] = fm[d - 1] + fm_change

        # Lean mass: assume constant (SURMOUNT-1: ~93% fat loss from this dataset)
        lean = start_lean - (start_fm - fm[d]) * (1 - FORBES_FAT_FRACTION) / FORBES_FAT_FRACTION
        weight[d] = fm[d] + lean

    return {
        "day": np.arange(n_days),
        "week": np.arange(n_days) / 7,
        "weight": weight,
        "fm": fm,
        "sp": sp,
        "intake": intake,
        "tdee": tdee,
        "eff_level": eff_level,
        "pct_change": (weight - weight[0]) / weight[0] * 100,
    }


def simulate_discontinuation(start_weight_lbs, target_mg, on_weeks=36, off_weeks=52,
                              start_fm_fraction=0.40, sp_gap_lbs=0.0,
                              on_drug_sp_hl_down=None):
    """Simulate SURMOUNT-4: on-drug phase then discontinuation."""
    # Phase 1: on-drug (may have slower SP adaptation)
    on_result = simulate(start_weight_lbs, target_mg, on_weeks, start_fm_fraction,
                         sp_gap_lbs=sp_gap_lbs,
                         on_drug_sp_hl_down=on_drug_sp_hl_down)

    # Phase 2: off-drug (continue from on-drug endpoint)
    n_off_days = off_weeks * 7
    fm_end = on_result["fm"][-1]
    sp_end = on_result["sp"][-1]
    weight_end = on_result["weight"][-1]
    lean_end = weight_end - fm_end

    fm = np.zeros(n_off_days)
    sp = np.zeros(n_off_days)
    weight = np.zeros(n_off_days)
    fm[0] = fm_end
    sp[0] = sp_end
    weight[0] = weight_end

    for d in range(1, n_off_days):
        # No drug — off-drug SP uses standard HL_DOWN=25d
        fm_direction = fm[d - 1] - fm[max(0, d - 7)]
        if fm_direction <= 0:
            alpha = 1 - np.exp(-np.log(2) / SP_HL_DOWN)
        else:
            alpha = 1 - np.exp(-np.log(2) / SP_HL_UP)
        sp[d] = sp[d - 1] + alpha * (fm[d - 1] - sp[d - 1])

        gap = sp[d] - fm[d - 1]
        eating_pressure = SP_EATING_PRESSURE * gap

        weight_ratio = weight[d - 1] / start_weight_lbs
        tdee_d = BASELINE_TDEE * (weight_ratio ** 0.75)

        intake_d = tdee_d + eating_pressure
        surplus = intake_d - tdee_d

        fm_change = surplus * FORBES_FAT_FRACTION / CAL_PER_LB_FAT
        fm[d] = fm[d - 1] + fm_change
        lean = lean_end  # constant lean
        weight[d] = fm[d] + lean

    return on_result, {
        "day": np.arange(n_off_days),
        "week": np.arange(n_off_days) / 7,
        "weight": weight,
        "fm": fm,
        "sp": sp,
        "pct_from_baseline": (weight - start_weight_lbs) / start_weight_lbs * 100,
        "pct_from_week36": (weight - weight_end) / weight_end * 100,
    }


def main():
    # SURMOUNT-1 parameters
    start_kg = 104.8
    start_lbs = start_kg * 2.205
    # BMI 38 at ~104.8 kg → approximate 40% body fat
    start_fm_frac = 0.40

    print("=" * 70)
    print("SURMOUNT-1 FORWARD SIMULATION")
    print(f"Start: {start_kg:.1f} kg ({start_lbs:.0f} lbs), ~{start_fm_frac*100:.0f}% body fat")
    print("Model parameters: all from this dataset (AG, AM, AQ, F)")
    print("Zero parameters fitted to trial data")
    print("=" * 70)

    # Simulate each dose arm
    results = {}
    for dose in [5, 10, 15]:
        r = simulate(start_lbs, dose, n_weeks=72, start_fm_fraction=start_fm_frac)
        results[dose] = r

    # Print trajectories at key weeks
    print(f"\n{'Week':>6} {'Placebo':>9} {'5mg':>9} {'10mg':>9} {'15mg':>9}  (% change from baseline)")
    for week in [0, 12, 24, 36, 48, 60, 72]:
        day = week * 7
        placebo_pct = -2.4 * week / 72  # linear approximation
        row = f"{week:>6} {placebo_pct:>+8.1f}%"
        for dose in [5, 10, 15]:
            if day < len(results[dose]["pct_change"]):
                pct = results[dose]["pct_change"][day]
                row += f" {pct:>+8.1f}%"
            else:
                row += f" {'n/a':>9}"
        print(row)

    # Published SURMOUNT-1 results
    print(f"\n{'Published SURMOUNT-1 (72 weeks)':}")
    print(f"  Placebo: -2.4%")
    print(f"  5mg:    -16.0%")
    print(f"  10mg:   -21.4%")
    print(f"  15mg:   -22.5%")

    # Compare
    print(f"\n{'Dose':>6} {'Simulated':>10} {'Published':>10} {'Δ':>8}")
    published = {5: -16.0, 10: -21.4, 15: -22.5}
    for dose in [5, 10, 15]:
        sim = results[dose]["pct_change"][-1]
        pub = published[dose]
        print(f"{dose:>4}mg {sim:>+9.1f}% {pub:>+9.1f}% {sim - pub:>+7.1f}%")

    # Detail: what's happening inside the model at week 72?
    print(f"\n--- Model internals at week 72 ---")
    print(f"{'':>6} {'FM (lbs)':>9} {'SP (lbs)':>9} {'Gap':>6} {'Eff lvl':>8} {'Intake':>7} {'TDEE':>6}")
    for dose in [5, 10, 15]:
        r = results[dose]
        d = -1
        gap = r["sp"][d] - r["fm"][d]
        print(f"{dose:>4}mg {r['fm'][d]:>9.1f} {r['sp'][d]:>9.1f} {gap:>+6.1f} {r['eff_level'][d]:>8.1f} "
              f"{r['intake'][d]:>7.0f} {r['tdee'][d]:>6.0f}")

    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("SURMOUNT-4 DISCONTINUATION SIMULATION")
    print("36 weeks on 15mg → 52 weeks off")
    print("=" * 70)

    on_result, off_result = simulate_discontinuation(
        start_lbs, target_mg=15, on_weeks=36, off_weeks=52,
        start_fm_fraction=start_fm_frac)
    # Also with best-fit gap
    on_g, off_g = simulate_discontinuation(
        start_lbs, target_mg=15, on_weeks=36, off_weeks=52,
        start_fm_fraction=start_fm_frac)

    print(f"\n  On-drug phase (0-36 weeks):")
    print(f"    Weight change: {on_result['pct_change'][-1]:+.1f}%")
    print(f"    FM: {on_result['fm'][0]:.0f} → {on_result['fm'][-1]:.0f} lbs")
    print(f"    SP at week 36: {on_result['sp'][-1]:.0f} lbs")
    print(f"    Gap at discontinuation: {on_result['sp'][-1] - on_result['fm'][-1]:+.1f} lbs")

    print(f"\n  Off-drug phase (36-88 weeks):")
    print(f"    Weight regain from week-36: {off_result['pct_from_week36'][-1]:+.1f}%")
    print(f"    Total from baseline: {off_result['pct_from_baseline'][-1]:+.1f}%")
    print(f"    FM: {off_result['fm'][0]:.0f} → {off_result['fm'][-1]:.0f} lbs")
    print(f"    SP at week 88: {off_result['sp'][-1]:.0f} lbs")

    print(f"\n  Published SURMOUNT-4:")
    print(f"    On-drug (36 wk): ~-20%")
    print(f"    Post-discontinuation regain (52 wk): +14%")

    # Sweep on-drug SP half-life for discontinuation
    print(f"\n  --- Discontinuation: on-drug SP HL determines regain ---")
    print(f"  {'On-drug HL':>11} {'On-drug %':>10} {'Regain %':>9} {'Gap at stop':>12}")
    for on_hl in [25, 50, 80, 120, 180]:
        on_r, off_r = simulate_discontinuation(
            start_lbs, target_mg=15, on_weeks=36, off_weeks=52,
            start_fm_fraction=start_fm_frac, on_drug_sp_hl_down=on_hl)
        gap_at_stop = on_r['sp'][-1] - on_r['fm'][-1]
        print(f"  {on_hl:>9}d {on_r['pct_change'][-1]:>+9.1f}% {off_r['pct_from_week36'][-1]:>+8.1f}% "
              f"{gap_at_stop:>+11.1f} lbs")
    print(f"  Published:     ~-20.0%    +14.0%")

    # Regain trajectory at on-drug HL=120d (best match for +14% regain)
    print(f"\n  --- Post-discontinuation trajectory (on-drug HL=120d) ---")
    on_r, off_r = simulate_discontinuation(
        start_lbs, target_mg=15, on_weeks=36, off_weeks=52,
        start_fm_fraction=start_fm_frac, on_drug_sp_hl_down=120)
    print(f"  {'Week':>6} {'Regain from wk36':>17} {'FM':>8} {'SP':>8} {'Gap':>6}")
    for wk in [0, 4, 8, 12, 24, 36, 52]:
        d = wk * 7
        if d >= len(off_r["weight"]):
            break
        regain_pct = off_r["pct_from_week36"][d]
        print(f"  {36 + wk:>6} {regain_pct:>+16.1f}% {off_r['fm'][d]:>8.1f} {off_r['sp'][d]:>8.1f} "
              f"{off_r['sp'][d] - off_r['fm'][d]:>+6.1f}")

    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("STARTING GAP SWEEP — reverse-engineering the hidden set point")
    print("=" * 70)

    print(f"\n  If trial subjects were slowly gaining (FM above SP), the set point's")
    print(f"  eating pressure provides a tailwind. How much gap reproduces the data?")
    print(f"\n  --- SURMOUNT-1 (non-diabetic, published: 5mg -16.0%, 10mg -21.4%, 15mg -22.5%) ---")
    print(f"  {'Gap (lbs)':>10} {'5mg':>8} {'10mg':>8} {'15mg':>8}")
    for gap in [0, -2, -3, -4, -5, -7, -10]:
        row = f"  {gap:>+8} lb"
        for dose in [5, 10, 15]:
            r = simulate(start_lbs, dose, n_weeks=72, start_fm_fraction=start_fm_frac,
                         sp_gap_lbs=gap)
            row += f" {r['pct_change'][-1]:>+7.1f}%"
        print(row)

    # Best-fit gap for SURMOUNT-1
    print(f"\n  Published:          -16.0%   -21.4%   -22.5%")

    # SURMOUNT-2 comparison
    s2_kg = 100.7
    s2_lbs = s2_kg * 2.205
    print(f"\n  --- SURMOUNT-2 (diabetic, start {s2_kg} kg, published: 10mg -12.8%, 15mg -14.7%) ---")
    print(f"  {'Gap (lbs)':>10} {'10mg':>8} {'15mg':>8}")
    for gap in [+3, +2, +1, 0, -1, -2, -3]:
        row = f"  {gap:>+8} lb"
        for dose in [10, 15]:
            r = simulate(s2_lbs, dose, n_weeks=72, start_fm_fraction=0.38,
                         sp_gap_lbs=gap)
            row += f" {r['pct_change'][-1]:>+7.1f}%"
        print(row)
    print(f"\n  Published:          -12.8%   -14.7%")

    # BMI subgroup comparison (pooled tirzepatide, from Horn et al.)
    print(f"\n  --- SURMOUNT-1 by BMI category (pooled tirz, Horn et al. 2025) ---")
    print(f"  {'BMI':>12} {'Start kg':>9} {'Published':>10} {'Sim gap=0':>10} {'Sim gap=-3':>11}")
    for bmi_label, start_kg_est, fm_frac, published in [
        ("27-30", 81, 0.30, -18.0),
        ("30-35", 92, 0.34, -21.5),
        ("35-40", 106, 0.38, -23.1),
        ("40+", 130, 0.44, -23.0),
    ]:
        s_lbs = start_kg_est * 2.205
        r0 = simulate(s_lbs, 10, n_weeks=72, start_fm_fraction=fm_frac, sp_gap_lbs=0)
        r3 = simulate(s_lbs, 10, n_weeks=72, start_fm_fraction=fm_frac, sp_gap_lbs=-3)
        print(f"  {bmi_label:>12} {start_kg_est:>7} kg {published:>+9.1f}% {r0['pct_change'][-1]:>+9.1f}% {r3['pct_change'][-1]:>+10.1f}%")

    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("WHY 15mg BARELY BEATS 10mg")
    print("=" * 70)

    gap_best = -3  # approximate best-fit gap for SURMOUNT-1
    for dose in [5, 10, 15]:
        r = simulate(start_lbs, dose, n_weeks=72, start_fm_fraction=start_fm_frac,
                     sp_gap_lbs=gap_best)
        weekly_pct = np.diff(r["pct_change"][::7])
        plateau_week = None
        for w in range(12, 72):
            if w < len(weekly_pct) and abs(weekly_pct[w]) < 0.1:
                plateau_week = w
                break
        eff_at_plateau = r["eff_level"][plateau_week * 7] if plateau_week else r["eff_level"][-1]
        gap_at_plateau = r["sp"][plateau_week * 7] - r["fm"][plateau_week * 7] if plateau_week else r["sp"][-1] - r["fm"][-1]
        print(f"  {dose:>2}mg: plateau ~week {plateau_week or '>72'}, eff_level={eff_at_plateau:.1f}, gap={gap_at_plateau:+.1f} lbs, final={r['pct_change'][-1]:+.1f}%")

    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)


if __name__ == "__main__":
    main()
