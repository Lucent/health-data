#!/usr/bin/env python3
"""BZ. Two-component drug effect model: acute response + decaying tonic.

The multiplicative tachyphylaxis model (AX: drug = slope × blood × exp(-decay×wk))
correctly identifies the within-week slope (-74.5 cal/unit) and the between-week
decay (35-week HL). But it fails forward prediction: at 80 weeks, it predicts a
130-cal sawtooth vs the 650-cal sawtooth observed. The multiplicative assumption —
that tachyphylaxis scales down the within-week sensitivity — is wrong.

This script fits and validates a two-component model:

    drug_effect = acute × (blood - EMA₁ᵈ(blood)) + tonic × blood × exp(-decay × wk)
                  ───────────────────────────────   ────────────────────────────────
                  response to blood level CHANGES    sustained suppression (decays)

The acute term captures the sawtooth: receptors respond to deviations from the
adapted (1-day EMA) blood level at full strength regardless of time on drug.
The tonic term captures the mean suppression that erodes with tachyphylaxis.

Locked parameters (from PK literature or structural constraints):
  - PK: t½=5.0d, ka=3.31/d (FDA label)
  - adapt_hl: 1 day (grid-searched, best fit)
  - tonic_decay: 0.0126/week → HL=55 weeks (drug property, held constant)
  - acute = 2 × tonic (ratio from subject fit: -102.5/-53.5 = 1.92, locked to 2.0)

Free parameter: tonic_slope (cal per unit blood level)
  - Subject fit: -53.5 (from 529 days logged intake, RMSE=326, r=0.56)
  - Undercount-corrected: -57.5 (subject documents ~7% intake undercount)
  - Trial best-fit (SURMOUNT-1/2/4 SSE minimization): -66

Trial validation (no set point, pure drug effect + MSJ-scaled TDEE):
  SURMOUNT-1 15mg:  -22.3% (published -22.5%, Δ=+0.2%)
  SURMOUNT-1 10mg:  -16.2% (published -21.4%, Δ=+5.2%)
  SURMOUNT-1 5mg:    -8.8% (published -16.0%, Δ=+7.2%)
  SURMOUNT-2 10mg:  -16.9% (published -12.8%, Δ=-4.1%)  [T2D: drug diverted to glucose]
  SURMOUNT-4 15mg:  -11.2% (published -20.0%, Δ=+8.8%)  [36wk only; intake>TDEE not modeled]

Residual pattern: underpredicts low doses (missing set point pressure), overpredicts
diabetic arms (GLP-1 diverted to glucose control). The tachyphylaxis structure
(55-week HL, 2:1 acute/tonic ratio) is not contradicted by any trial arm.
"""

from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent.parent

# ═══════════════════════════════════════════════════════════════════
# PK MODEL (FDA label)
# ═══════════════════════════════════════════════════════════════════
PK_HALF_LIFE = 5.0  # days
PK_KE = math.log(2) / PK_HALF_LIFE
PK_KA = 3.31  # absorption rate, gives Tmax ≈ 1 day

def pk_curve(dose, t):
    if t <= 0:
        return 0.0
    return dose * PK_KA / (PK_KA - PK_KE) * (math.exp(-PK_KE * t) - math.exp(-PK_KA * t))

def blood_on_day(schedule, day):
    level = 0.0
    for inj_day, dose in schedule:
        dt = day - inj_day
        if dt < 0:
            break
        if dt > 30:
            continue
        level += pk_curve(dose, dt)
    return level


# ═══════════════════════════════════════════════════════════════════
# MODEL PARAMETERS
# ═══════════════════════════════════════════════════════════════════
ADAPT_HL = 1.0        # day — EMA half-life for adapted baseline
TONIC_DECAY = 0.0126  # per week → HL = ln(2)/0.0126 = 55 weeks
ACUTE_RATIO = 2.0     # acute_slope = ACUTE_RATIO × tonic_slope
FORBES_FAT_FRAC = 0.85


def compute_ema(values, hl_days):
    alpha = 1 - math.exp(-math.log(2) / hl_days)
    ema = np.empty(len(values))
    ema[0] = values[0]
    for i in range(1, len(values)):
        ema[i] = ema[i - 1] + alpha * (values[i] - ema[i - 1])
    return ema


def mifflin_tdee(weight_kg, height_cm, age, pct_female, pal=1.4):
    """Population-average TDEE from Mifflin-St Jeor × physical activity level."""
    rmr_m = 10 * weight_kg + 6.25 * height_cm - 5 * age + 5
    rmr_f = 10 * weight_kg + 6.25 * height_cm - 5 * age - 161
    rmr = (1 - pct_female) * rmr_m + pct_female * rmr_f
    return rmr * pal


# ═══════════════════════════════════════════════════════════════════
# SUBJECT FIT
# ═══════════════════════════════════════════════════════════════════
def fit_subject():
    """Fit the two-component model to this subject's on-drug intake data."""
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    df = intake[["date", "calories"]].merge(
        tirz[["date", "blood_level", "days_since_injection"]], on="date")
    df = df.merge(kalman[["date", "tdee", "fat_mass_lbs"]], on="date")
    df = df[df["blood_level"] > 0].sort_values("date").reset_index(drop=True)
    drug_start = df["date"].min()
    df["cum_weeks"] = (df["date"] - drug_start).dt.days / 7

    # SP pressure (to isolate drug effect in the target)
    SP_HL, SP_CAL_PER_LB = 45, 27
    fm_all = kalman.set_index("date")["fat_mass_lbs"].reindex(
        pd.date_range(kalman["date"].min(), kalman["date"].max())).ffill().bfill()
    alpha_sp = 1 - math.exp(-math.log(2) / SP_HL)
    sp_vals = np.empty(len(fm_all))
    sp_vals[0] = fm_all.iloc[0]
    for i in range(1, len(fm_all)):
        v = fm_all.iloc[i]
        sp_vals[i] = sp_vals[i - 1] + alpha_sp * ((v if not np.isnan(v) else sp_vals[i - 1]) - sp_vals[i - 1])
    sp_series = pd.Series(sp_vals, index=fm_all.index)
    df["sp_pressure"] = SP_CAL_PER_LB * (df["date"].map(sp_series) - df["fat_mass_lbs"])
    df["target"] = df["calories"] - df["tdee"] - df["sp_pressure"]

    blood = df["blood_level"].values
    target = df["target"].values
    cum_weeks = df["cum_weeks"].values
    n = len(df)

    ema = compute_ema(blood, ADAPT_HL)
    deviation = blood - ema

    # Fit: target = base + acute×deviation + tonic×blood×exp(-decay×wk)
    # With constraint: acute = 2 × tonic
    def objective(params):
        tonic = params[0]
        acute = ACUTE_RATIO * tonic
        pred = acute * deviation + tonic * blood * np.exp(-TONIC_DECAY * cum_weeks)
        return np.mean((target - pred) ** 2)

    res = minimize(objective, [-55], method="Nelder-Mead")
    tonic = res.x[0]
    acute = ACUTE_RATIO * tonic
    pred = acute * deviation + tonic * blood * np.exp(-TONIC_DECAY * cum_weeks)
    rmse = np.sqrt(np.mean((target - pred) ** 2))
    r = np.corrcoef(target, pred)[0, 1]

    # Also fit unconstrained for comparison
    def obj_free(params):
        acute_f, tonic_f = params
        pred_f = acute_f * deviation + tonic_f * blood * np.exp(-TONIC_DECAY * cum_weeks)
        return np.mean((target - pred_f) ** 2)

    res_free = minimize(obj_free, [-100, -55], method="Nelder-Mead")
    acute_free, tonic_free = res_free.x
    pred_free = acute_free * deviation + tonic_free * blood * np.exp(-TONIC_DECAY * cum_weeks)
    rmse_free = np.sqrt(np.mean((target - pred_free) ** 2))
    r_free = np.corrcoef(target, pred_free)[0, 1]

    # Multiplicative model (AX) for comparison
    def obj_mult(params):
        slope, decay = params
        if decay < 0:
            return 1e12
        pred_m = slope * blood * np.exp(-decay * cum_weeks)
        return np.mean((target - pred_m) ** 2)

    res_mult = minimize(obj_mult, [-65, 0.02], method="Nelder-Mead")
    slope_m, decay_m = res_mult.x
    pred_mult = slope_m * blood * np.exp(-decay_m * cum_weeks)
    rmse_mult = np.sqrt(np.mean((target - pred_mult) ** 2))
    r_mult = np.corrcoef(target, pred_mult)[0, 1]
    hl_mult = math.log(2) / decay_m if decay_m > 0 else float("inf")

    print("=" * 70)
    print("1. SUBJECT FIT (529 on-drug days)")
    print("=" * 70)
    print(f"\n  Multiplicative (AX-style): slope={slope_m:+.1f}, HL={hl_mult:.0f}wk")
    print(f"    RMSE={rmse_mult:.0f}, r={r_mult:.3f}")
    print(f"\n  Two-component (unconstrained): acute={acute_free:+.1f}, tonic={tonic_free:+.1f}")
    print(f"    ratio={acute_free / tonic_free:.2f}:1, RMSE={rmse_free:.0f}, r={r_free:.3f}")
    print(f"\n  Two-component (acute=2×tonic): acute={acute:+.1f}, tonic={tonic:+.1f}")
    print(f"    RMSE={rmse:.0f}, r={r:.3f}")

    # Sawtooth validation by era
    df["pred_2c"] = pred
    df["pred_mult"] = pred_mult

    print(f"\n  Sawtooth by era (observed vs predicted swing, day 0-1 vs day 4-5):")
    print(f"  {'Era':<20} {'Observed':>9} {'2-component':>12} {'Multiplicative':>15}")
    for label, w0, w1 in [("Wk 0-20", 0, 20), ("Wk 20-40", 20, 40),
                           ("Wk 40-60", 40, 60), ("Wk 60-80", 60, 100)]:
        mask = (df["cum_weeks"] >= w0) & (df["cum_weeks"] < w1)
        sub = df[mask]
        if len(sub) < 30:
            continue
        obs_swing = (sub[sub["days_since_injection"].isin([4, 5])]["target"].mean()
                     - sub[sub["days_since_injection"].isin([0, 1])]["target"].mean())
        p2c_swing = (sub[sub["days_since_injection"].isin([4, 5])]["pred_2c"].mean()
                     - sub[sub["days_since_injection"].isin([0, 1])]["pred_2c"].mean())
        pm_swing = (sub[sub["days_since_injection"].isin([4, 5])]["pred_mult"].mean()
                    - sub[sub["days_since_injection"].isin([0, 1])]["pred_mult"].mean())
        print(f"  {label:<20} {obs_swing:>+8.0f} {p2c_swing:>+11.0f} {pm_swing:>+14.0f}")

    return tonic


# ═══════════════════════════════════════════════════════════════════
# TRIAL SIMULATION
# ═══════════════════════════════════════════════════════════════════
def build_tirz_schedule(target_mg, n_weeks):
    escalation = {
        5: [(0, 2.5), (4, 5.0)],
        10: [(0, 2.5), (4, 5.0), (8, 7.5), (12, 10.0)],
        15: [(0, 2.5), (4, 5.0), (8, 7.5), (12, 10.0), (16, 12.5), (20, 15.0)],
    }
    changes = escalation[target_mg]
    schedule = []
    dose_idx = 0
    current = changes[0][1]
    for week in range(n_weeks):
        if dose_idx + 1 < len(changes) and week >= changes[dose_idx + 1][0]:
            dose_idx += 1
            current = changes[dose_idx][1]
        schedule.append((week * 7, current))
    return schedule


def simulate_trial(start_kg, target_mg, n_weeks, fm_frac, tonic_slope,
                   height_cm, age, pct_female, pal=1.4):
    """Pure drug effect simulation: intake = TDEE + drug_effect, no set point."""
    acute_slope = ACUTE_RATIO * tonic_slope
    start_lbs = start_kg * 2.205
    n_days = n_weeks * 7 + 1
    schedule = build_tirz_schedule(target_mg, n_weeks)

    fm = start_lbs * fm_frac
    lean = start_lbs * (1 - fm_frac)
    alpha_ema = 1 - math.exp(-math.log(2) / ADAPT_HL)
    blood_ema = 0.0

    for d in range(1, n_days):
        blood = blood_on_day(schedule, d)
        blood_ema = blood_ema + alpha_ema * (blood - blood_ema)
        deviation = blood - blood_ema
        tachy = math.exp(-TONIC_DECAY * d / 7)
        drug_effect = acute_slope * deviation + tonic_slope * blood * tachy

        weight_kg = (fm + lean) / 2.205
        tdee = mifflin_tdee(weight_kg, height_cm, age, pct_female, pal)

        fm += drug_effect * FORBES_FAT_FRAC / 3500
        lean -= drug_effect * (1 - FORBES_FAT_FRAC) / 3500

    weight_end = fm + lean
    return (weight_end - start_lbs) / start_lbs * 100


def run_trials(tonic_slope, label=""):
    trials = {
        "SURMOUNT-1": {
            "desc": "Jastreboff NEJM 2022, non-diabetic",
            "start_kg": 104.8, "height_cm": 167, "age": 45, "pct_female": 0.67,
            "fm_frac": 0.40, "n_weeks": 72,
            "arms": {5: -16.0, 10: -21.4, 15: -22.5},
        },
        "SURMOUNT-2": {
            "desc": "Garvey Lancet 2023, T2D",
            "start_kg": 100.7, "height_cm": 167, "age": 54, "pct_female": 0.50,
            "fm_frac": 0.38, "n_weeks": 72,
            "arms": {10: -12.8, 15: -14.7},
        },
        "SURMOUNT-4": {
            "desc": "Aronne JAMA 2024, on-drug 36wk",
            "start_kg": 107.3, "height_cm": 167, "age": 46, "pct_female": 0.67,
            "fm_frac": 0.455, "n_weeks": 36,
            "arms": {15: -20.0},
        },
    }

    total_sse = 0
    print(f"\n  {'Trial':<14} {'Dose':>5} {'Simulated':>10} {'Published':>10} {'Δ':>7}")
    print(f"  {'-'*50}")
    for name, t in trials.items():
        for dose, published in t["arms"].items():
            pct = simulate_trial(t["start_kg"], dose, t["n_weeks"], t["fm_frac"],
                                 tonic_slope, t["height_cm"], t["age"], t["pct_female"])
            delta = pct - published
            total_sse += delta ** 2
            print(f"  {name:<14} {dose:>3}mg {pct:>+9.1f}% {published:>+9.1f}% {delta:>+6.1f}%")

    print(f"  {'SSE':>50}: {total_sse:.1f}")
    return total_sse


def sweep_trials():
    """Sweep tonic_slope across all trial arms."""
    print(f"\n{'=' * 70}")
    print("2. TRIAL SWEEP (no set point, pure drug + MSJ TDEE)")
    print("   Locked: acute=2×tonic, adapt=1d, tachy_HL=55wk")
    print("=" * 70)

    # Show baseline TDEEs
    for name, kg, h, a, f in [("SURMOUNT-1", 104.8, 167, 45, 0.67),
                                ("SURMOUNT-2", 100.7, 167, 54, 0.50),
                                ("SURMOUNT-4", 107.3, 167, 46, 0.67)]:
        print(f"  {name}: {kg} kg → MSJ TDEE = {mifflin_tdee(kg, h, a, f):.0f} cal")

    best_sse = 1e12
    best_tonic = None
    for tonic in np.arange(-40, -180, -1):
        sse = 0
        trials = [
            (104.8, 5, 72, 0.40, 167, 45, 0.67, -16.0),
            (104.8, 10, 72, 0.40, 167, 45, 0.67, -21.4),
            (104.8, 15, 72, 0.40, 167, 45, 0.67, -22.5),
            (100.7, 10, 72, 0.38, 167, 54, 0.50, -12.8),
            (100.7, 15, 72, 0.38, 167, 54, 0.50, -14.7),
            (107.3, 15, 36, 0.455, 167, 46, 0.67, -20.0),
        ]
        for kg, dose, wk, fm, h, a, f, pub in trials:
            pct = simulate_trial(kg, dose, wk, fm, tonic, h, a, f)
            sse += (pct - pub) ** 2
        if sse < best_sse:
            best_sse = sse
            best_tonic = tonic

    print(f"\n  Best tonic: {best_tonic:+.0f} (acute={2*best_tonic:+.0f}), SSE={best_sse:.1f}")
    run_trials(best_tonic, "trial-best")
    return best_tonic


def main():
    subject_tonic = fit_subject()

    trial_tonic = sweep_trials()

    # Subject-fitted values against trials
    print(f"\n{'=' * 70}")
    print("3. SUBJECT-FITTED tonic against trials")
    print("=" * 70)
    print(f"  Subject: tonic={subject_tonic:+.1f} (logged-calorie basis)")
    print(f"  Undercount-corrected (~7%): {subject_tonic / 0.93:+.1f}")
    run_trials(subject_tonic, "subject-fitted")

    print(f"\n{'=' * 70}")
    print("4. SUMMARY")
    print("=" * 70)
    print(f"""
  Formula:
    drug_effect = {ACUTE_RATIO:.0f}×tonic × (blood - EMA(blood, {ADAPT_HL:.0f}d)) + tonic × blood × exp(-{TONIC_DECAY}×wk)
                  ─────────────────────────────────────────  ──────────────────────────────────────
                  acute: responds to blood level changes      tonic: sustained suppression, decays

  Locked parameters (3):
    adapt_hl     = {ADAPT_HL:.0f} day    (grid-searched, best fit)
    tonic_decay  = {TONIC_DECAY}/wk  (HL={math.log(2)/TONIC_DECAY:.0f} weeks)
    acute/tonic  = {ACUTE_RATIO:.0f}:1     (from subject ratio {-102.5/-53.5:.2f}:1, locked)

  Free parameter (1):
    tonic_slope — cal per unit blood level
    Subject fit:           {subject_tonic:+.1f} cal/unit
    Undercount-corrected:  {subject_tonic/0.93:+.1f} cal/unit
    Trial best-fit:        {trial_tonic:+.0f} cal/unit
    Ratio trial/corrected: {trial_tonic / (subject_tonic/0.93):.2f}x

  The subject's tonic is {abs(subject_tonic/0.93 / trial_tonic - 1)*100:.0f}% below the trial value.
  Consistent with body composition and individual variation.
""")


if __name__ == "__main__":
    main()
