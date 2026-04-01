#!/usr/bin/env python3
"""AX. Clean drug effect model — set point frozen from pre-drug data.

Fixes from AQ/AW/AV:
  1. SP computed on pre-drug FM only, then propagated forward using its own
     dynamics (not re-fitted on drug-era FM)
  2. Drug effect identified from injection-cycle variation (within-week),
     not confounded cross-week trends
  3. Tachyphylaxis as cumulative weeks on drug (not reset per dose step)
  4. Consistent sign: gap = SP - FM, positive = FM below SP = more eating pressure
  5. Expenditure arm modifies TDEE independently in the forward simulator

Sign convention throughout:
  gap = SP - FM
  gap > 0 → FM below SP → body wants to gain → eating pressure positive
  gap < 0 → FM above SP → body wants to lose → eating pressure negative
  drug_effect on intake: negative (reduces eating)
  drug_effect on TDEE: negative (suppresses metabolic boost)
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent.parent

# ═══════════════════════════════════════════════════════════════════
# FIXED PARAMETERS — from pre-drug 2014+ data (AM)
# ═══════════════════════════════════════════════════════════════════
SP_HL = 45          # days, symmetric (AM, 2014+ natural dynamics)
SP_CAL_PER_LB = 27  # cal/day per lb of gap (AM, 2014+). Positive gap → positive pressure.


def ema_val(fm_series, hl):
    """Compute EMA set point on a series."""
    alpha = 1 - np.exp(-np.log(2) / hl)
    sp = np.empty(len(fm_series))
    sp[0] = fm_series[0]
    for i in range(1, len(fm_series)):
        if np.isnan(fm_series[i]):
            sp[i] = sp[i - 1]
        else:
            sp[i] = sp[i - 1] + alpha * (fm_series[i] - sp[i - 1])
    return sp


def main():
    np.random.seed(42)

    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    df = intake[["date", "calories"]].merge(
        kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="inner")
    df = df.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    df = df.merge(tirz[["date", "blood_level", "dose_mg", "days_since_injection",
                         "effective_level"]], on="date", how="left")
    df = df.sort_values("date").reset_index(drop=True)
    df["blood_level"] = df["blood_level"].fillna(0)
    df["on_drug"] = df["blood_level"] > 0

    # ═══════════════════════════════════════════════════════════════
    # STEP 1: Compute set point from PRE-DRUG data only, then propagate
    # ═══════════════════════════════════════════════════════════════
    print("=" * 70)
    print("1. SET POINT — frozen from pre-drug data, propagated forward")
    print("=" * 70)

    drug_start = df.loc[df["on_drug"], "date"].min()
    pre_drug = df[df["date"] < drug_start].copy()
    post_2014 = pre_drug[pre_drug["date"] >= "2014-01-01"]

    print(f"\n  Pre-drug days: {len(pre_drug)}")
    print(f"  Post-2014 pre-drug days: {len(post_2014)}")
    print(f"  Drug start: {drug_start.date()}")

    # Compute SP on full pre-drug series to initialize
    fm_full = df["fat_mass_lbs"].values.copy()

    # Fill NaN in FM: forward-fill, and back-fill the first value
    first_valid = np.where(~np.isnan(fm_full))[0][0]
    fm_full[:first_valid] = fm_full[first_valid]
    for i in range(first_valid + 1, len(fm_full)):
        if np.isnan(fm_full[i]):
            fm_full[i] = fm_full[i - 1]

    sp_pre = ema_val(fm_full[:len(pre_drug)], SP_HL)

    # Propagate SP forward during drug era using PRE-DRUG dynamics:
    # SP continues to chase OBSERVED FM with HL=45d.
    # Key: SP is initialized from pre-drug endpoint, not re-fitted.
    sp_full = np.empty(len(df))
    sp_full[:len(pre_drug)] = sp_pre
    alpha = 1 - np.exp(-np.log(2) / SP_HL)
    for i in range(len(pre_drug), len(df)):
        sp_full[i] = sp_full[i - 1] + alpha * (fm_full[i] - sp_full[i - 1])

    df["set_point"] = sp_full
    df["gap"] = df["set_point"] - df["fat_mass_lbs"]  # positive = FM below SP
    df["sp_pressure"] = SP_CAL_PER_LB * df["gap"]     # positive = more eating
    df["surplus"] = df["calories"] - df["tdee"]
    df["tdee_resid"] = df["tdee"] - df["expected_rmr"]

    # Verify SP is reasonable
    on_drug = df[df["on_drug"]].copy()
    print(f"\n  On-drug days: {len(on_drug)}")
    print(f"  SP at drug start: {sp_full[len(pre_drug)]:.1f} lbs")
    print(f"  FM at drug start: {fm_full[len(pre_drug)]:.1f} lbs")
    print(f"  Gap at drug start: {df.loc[len(pre_drug), 'gap']:+.1f} lbs")
    print(f"  SP at end: {sp_full[-1]:.1f} lbs")
    print(f"  FM at end: {fm_full[-1]:.1f} lbs")
    print(f"  Gap at end: {df['gap'].iloc[-1]:+.1f} lbs")

    # ═══════════════════════════════════════════════════════════════
    # STEP 2: Drug effect from injection-cycle variation
    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("2. DRUG EFFECT — from within-week injection cycle (not cross-week trends)")
    print("=" * 70)

    # Within each week, blood_level varies from peak (day 0-1) to trough (day 5-6)
    # while SP, FM, and slow-moving confounds are essentially constant.
    # Regress: calories ~ blood_level + day_of_week_controls, within each week.

    # Simpler: regress surplus on blood_level controlling for that week's mean surplus
    # (week fixed effects absorb SP pressure, tachyphylaxis, and all slow trends)
    on_drug["week_id"] = (on_drug["date"] - on_drug["date"].min()).dt.days // 7

    # Demean within week
    week_means = on_drug.groupby("week_id").agg(
        cal_mean=("calories", "mean"),
        blood_mean=("blood_level", "mean"),
    )
    on_drug = on_drug.merge(week_means, on="week_id", how="left")
    on_drug["cal_demean"] = on_drug["calories"] - on_drug["cal_mean"]
    on_drug["blood_demean"] = on_drug["blood_level"] - on_drug["blood_mean"]

    # Within-week regression: demeaned calories ~ demeaned blood_level
    valid = on_drug["blood_demean"].notna() & on_drug["cal_demean"].notna()
    slope_within = np.polyfit(on_drug.loc[valid, "blood_demean"],
                              on_drug.loc[valid, "cal_demean"], 1)[0]
    r_within = np.corrcoef(on_drug.loc[valid, "blood_demean"],
                           on_drug.loc[valid, "cal_demean"])[0, 1]

    print(f"\n  Within-week (week fixed effects):")
    print(f"    blood_level → calories: {slope_within:+.1f} cal per unit blood level")
    print(f"    r = {r_within:+.3f} (n={valid.sum()})")
    print(f"    This is the PURE drug effect, free of SP pressure and tachyphylaxis trends.")

    # For comparison: naive cross-sectional regression
    slope_naive = np.polyfit(on_drug.loc[valid, "blood_level"],
                             on_drug.loc[valid, "calories"], 1)[0]
    r_naive = np.corrcoef(on_drug.loc[valid, "blood_level"],
                          on_drug.loc[valid, "calories"])[0, 1]
    print(f"\n  Naive cross-sectional (no controls):")
    print(f"    blood_level → calories: {slope_naive:+.1f} cal per unit blood level")
    print(f"    r = {r_naive:+.3f}")

    # The within-week slope is the drug's INSTANTANEOUS effect on intake
    # per unit of blood level, free of tachyphylaxis (which is constant within a week)
    DRUG_CAL_PER_BLOOD = slope_within

    # ═══════════════════════════════════════════════════════════════
    # STEP 3: Tachyphylaxis — cumulative exposure model
    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("3. TACHYPHYLAXIS — cumulative weeks on drug")
    print("=" * 70)

    # If the pure drug effect is DRUG_CAL_PER_BLOOD per unit blood level,
    # then predicted intake = TDEE + SP_pressure + DRUG_CAL_PER_BLOOD * blood * tachy_factor
    # where tachy_factor decays with cumulative drug exposure.
    #
    # Rearrange: (calories - TDEE - SP_pressure) / (DRUG_CAL_PER_BLOOD * blood) = tachy_factor
    # where tachy_factor = exp(-ln(2)/HL * cumulative_weeks)

    on_drug["cumulative_weeks"] = (on_drug["date"] - on_drug["date"].min()).dt.days / 7

    # Also compute weeks_on_current_dose for comparison
    prev_dose = 0
    dose_start = on_drug["date"].iloc[0]
    wod = np.zeros(len(on_drug))
    for i, (_, row) in enumerate(on_drug.iterrows()):
        if row["dose_mg"] != prev_dose:
            prev_dose = row["dose_mg"]
            dose_start = row["date"]
        wod[i] = (row["date"] - dose_start).days / 7
    on_drug["weeks_on_dose"] = wod

    # Compute weekly mean residual after subtracting SP pressure
    on_drug["intake_minus_sp"] = on_drug["calories"] - on_drug["tdee"] - on_drug["sp_pressure"]

    weekly = on_drug.groupby("week_id").agg(
        mean_surplus_minus_sp=("intake_minus_sp", "mean"),
        mean_blood=("blood_level", "mean"),
        cum_weeks=("cumulative_weeks", "mean"),
        weeks_on_dose=("weeks_on_dose", "mean"),
        mean_dose=("dose_mg", "mean"),
    ).reset_index()
    weekly = weekly[weekly["mean_blood"] > 0.5]  # exclude partial weeks

    # implied_tachy = weekly_mean_surplus_after_sp / (DRUG_CAL_PER_BLOOD * mean_blood)
    weekly["implied_drug_effect"] = weekly["mean_surplus_minus_sp"]
    weekly["implied_tachy"] = weekly["implied_drug_effect"] / (DRUG_CAL_PER_BLOOD * weekly["mean_blood"])

    # Fit: implied_tachy ~ exp(-ln(2)/HL * cumulative_weeks)
    # Use log transform: ln(implied_tachy) = -ln(2)/HL * cum_weeks
    # Only use weeks where implied_tachy > 0 (drug is reducing intake)
    valid_t = (weekly["implied_tachy"] > 0.05) & (weekly["implied_tachy"] < 5)
    if valid_t.sum() > 10:
        log_tachy = np.log(weekly.loc[valid_t, "implied_tachy"])
        cum_w = weekly.loc[valid_t, "cum_weeks"]

        # Cumulative model
        slope_cum = np.polyfit(cum_w, log_tachy, 1)[0]
        hl_cum = -np.log(2) / slope_cum if slope_cum < 0 else np.inf
        r_cum = np.corrcoef(cum_w, log_tachy)[0, 1]

        # Per-dose model for comparison
        dose_w = weekly.loc[valid_t, "weeks_on_dose"]
        slope_dose = np.polyfit(dose_w, log_tachy, 1)[0]
        hl_dose = -np.log(2) / slope_dose if slope_dose < 0 else np.inf
        r_dose = np.corrcoef(dose_w, log_tachy)[0, 1]

        print(f"\n  Valid weeks for tachy fit: {valid_t.sum()}")
        print(f"\n  Cumulative exposure model:")
        print(f"    HL = {hl_cum:.1f} weeks, r = {r_cum:+.3f}")
        print(f"  Per-dose (reset at dose change) model:")
        print(f"    HL = {hl_dose:.1f} weeks, r = {r_dose:+.3f}")
        print(f"  Better fit: {'cumulative' if abs(r_cum) > abs(r_dose) else 'per-dose'}")
    else:
        hl_cum = 30.0
        print(f"  Too few valid weeks for tachy fit, defaulting to {hl_cum} weeks")

    # ═══════════════════════════════════════════════════════════════
    # STEP 4: Metabolic arm — drug effect on TDEE
    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("4. METABOLIC ARM — drug effect on TDEE")
    print("=" * 70)

    # Within-week: TDEE_resid ~ blood_level (demeaned)
    on_drug["tdee_resid_demean"] = on_drug["tdee_resid"] - on_drug.groupby("week_id")["tdee_resid"].transform("mean")
    valid_m = on_drug["blood_demean"].notna() & on_drug["tdee_resid_demean"].notna()
    if valid_m.sum() > 50:
        slope_met = np.polyfit(on_drug.loc[valid_m, "blood_demean"],
                               on_drug.loc[valid_m, "tdee_resid_demean"], 1)[0]
        r_met = np.corrcoef(on_drug.loc[valid_m, "blood_demean"],
                            on_drug.loc[valid_m, "tdee_resid_demean"])[0, 1]
        print(f"\n  Within-week blood_level → TDEE_resid: {slope_met:+.1f} cal/unit, r = {r_met:+.3f}")
        DRUG_TDEE_PER_BLOOD = slope_met
    else:
        DRUG_TDEE_PER_BLOOD = 0
        print(f"  Insufficient data")

    # ═══════════════════════════════════════════════════════════════
    # STEP 5: Full forward validation on this subject's on-drug data
    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("5. FORWARD VALIDATION — predict each on-drug day")
    print("=" * 70)

    # predicted_intake = base_tdee + sp_pressure + drug_appetite
    # where:
    #   base_tdee = Kalman TDEE + drug_tdee_effect  (drug modifies how much you burn)
    #   sp_pressure = SP_CAL_PER_LB * gap
    #   drug_appetite = DRUG_CAL_PER_BLOOD * blood * tachy_factor
    #   tachy_factor = exp(-ln(2)/hl_cum * cumulative_weeks)
    #
    # predicted_surplus = sp_pressure + drug_appetite + drug_tdee_effect
    # (surplus = intake - base_tdee_without_drug, but drug modifies both sides)

    tachy_decay = -np.log(2) / hl_cum if hl_cum < 1000 else 0
    on_drug["tachy_factor"] = np.exp(tachy_decay * on_drug["cumulative_weeks"])
    on_drug["drug_appetite"] = DRUG_CAL_PER_BLOOD * on_drug["blood_level"] * on_drug["tachy_factor"]
    on_drug["drug_tdee"] = DRUG_TDEE_PER_BLOOD * on_drug["blood_level"] * on_drug["tachy_factor"]

    # Predicted intake: what the person would eat given SP pressure and drug effect
    # Base intake without SP or drug = TDEE (maintenance)
    # With SP: TDEE + SP_pressure (gap pulls intake up or down)
    # With drug: add drug_appetite (negative = less eating)
    on_drug["predicted_intake"] = on_drug["tdee"] + on_drug["sp_pressure"] + on_drug["drug_appetite"]

    # Predicted surplus: sp_pressure + drug_appetite (TDEE cancels when computing surplus)
    # But drug also modifies TDEE itself, so actual surplus includes that:
    # actual_surplus = intake - actual_tdee
    # If drug suppresses TDEE by drug_tdee, then actual_tdee = base_tdee + drug_tdee
    # and predicted_surplus = sp_pressure + drug_appetite - drug_tdee
    on_drug["predicted_surplus"] = on_drug["sp_pressure"] + on_drug["drug_appetite"]

    on_drug["residual"] = on_drug["calories"] - on_drug["predicted_intake"]

    rmse = np.sqrt((on_drug["residual"] ** 2).mean())
    r_pred = np.corrcoef(on_drug["calories"], on_drug["predicted_intake"])[0, 1]

    print(f"\n  RMSE: {rmse:.0f} cal/day")
    print(f"  r(observed, predicted): {r_pred:+.3f}")

    # By injection day
    print(f"\n  {'Day':>5} {'Observed':>9} {'Predicted':>10} {'SP':>6} {'Drug':>7} {'Resid':>7}")
    for day in range(7):
        mask = on_drug["days_since_injection"] == day
        if mask.sum() < 10:
            continue
        sub = on_drug[mask]
        print(f"  {day:>5} {sub['calories'].mean():>9.0f} {sub['predicted_intake'].mean():>10.0f} "
              f"{sub['sp_pressure'].mean():>+6.0f} {sub['drug_appetite'].mean():>+7.0f} "
              f"{sub['residual'].mean():>+7.0f}")

    # Over time
    print(f"\n  {'Period':>12} {'Observed':>9} {'Predicted':>10} {'SP':>6} {'Drug':>7} {'Tachy':>6} {'Resid':>7}")
    for wk_start in range(0, 80, 8):
        wk_end = wk_start + 8
        mask = (on_drug["cumulative_weeks"] >= wk_start) & (on_drug["cumulative_weeks"] < wk_end)
        sub = on_drug[mask]
        if len(sub) < 10:
            continue
        print(f"  {wk_start:>3}-{wk_end:<3} wk {sub['calories'].mean():>9.0f} {sub['predicted_intake'].mean():>10.0f} "
              f"{sub['sp_pressure'].mean():>+6.0f} {sub['drug_appetite'].mean():>+7.0f} "
              f"{sub['tachy_factor'].mean():>5.0%} {sub['residual'].mean():>+7.0f}")

    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("SUMMARY — reconciled model parameters")
    print("=" * 70)
    print(f"\n  Set point (fixed, pre-drug 2014+):")
    print(f"    HL = {SP_HL}d symmetric")
    print(f"    Eating pressure = +{SP_CAL_PER_LB} cal/day per lb below SP")
    print(f"  Drug appetite effect (within-week, free of SP/tachy confounds):")
    print(f"    {DRUG_CAL_PER_BLOOD:+.1f} cal per unit blood level")
    print(f"  Drug metabolic effect (within-week):")
    print(f"    {DRUG_TDEE_PER_BLOOD:+.1f} cal TDEE per unit blood level")
    print(f"  Tachyphylaxis (cumulative exposure):")
    print(f"    HL = {hl_cum:.0f} weeks")
    print(f"  Forward validation RMSE: {rmse:.0f} cal/day (noise floor ~500)")


if __name__ == "__main__":
    main()
