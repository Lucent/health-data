#!/usr/bin/env python3
"""AW. Reconcile set point and tachyphylaxis models on tirzepatide data.

Use pre-tirz-derived set point parameters (HL=45d, -27 cal/lb from AM 2014+)
as fixed inputs. Then fit the drug parameters (per-unit appetite effect and
tachyphylaxis half-life) to minimize prediction error on the 529 on-drug days.

Forward model for each on-drug day:
  predicted_intake = TDEE + SP_pressure + drug_effect
where:
  SP_pressure = -27 * (SP - FM)    [from pre-tirz, fixed]
  drug_effect = appetite_per_unit * effective_level
  effective_level = blood_level * exp(-ln(2)/tachy_hl * weeks_on_current_dose)

Free parameters: appetite_per_unit, tachy_hl
Fixed parameters: SP HL=45d (symmetric), SP pressure=-27 cal/lb
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent.parent

SP_HL = 45  # days, symmetric (AM, 2014+ data)
SP_PRESSURE = -27  # cal/day per lb of gap (AM, 2014+ data)


def ema(series, hl):
    alpha = 1 - np.exp(-np.log(2) / hl)
    return series.ewm(alpha=alpha, min_periods=30).mean()


def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    # Build full timeline
    df = intake[["date", "calories"]].merge(
        kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="inner")
    df = df.merge(tirz[["date", "blood_level", "dose_mg", "days_since_injection"]],
                  on="date", how="left")
    df = df.sort_values("date").reset_index(drop=True)

    # Set point: EMA of FM at HL=45d on full data
    df["set_point"] = ema(df["fat_mass_lbs"], SP_HL).values
    df["sp_gap"] = df["set_point"] - df["fat_mass_lbs"]  # positive = FM below SP
    df["sp_pressure"] = SP_PRESSURE * df["sp_gap"]

    # On-drug subset
    on = df[df["blood_level"].notna() & (df["blood_level"] > 0)].copy()
    on = on.dropna(subset=["calories", "tdee", "fat_mass_lbs"])
    print(f"On-drug days: {len(on)}")

    # Track weeks on current dose for tachyphylaxis
    prev_dose = 0
    dose_start_date = on["date"].iloc[0]
    weeks_on_dose = np.zeros(len(on))
    for i, (_, row) in enumerate(on.iterrows()):
        if row["dose_mg"] != prev_dose:
            prev_dose = row["dose_mg"]
            dose_start_date = row["date"]
        weeks_on_dose[i] = (row["date"] - dose_start_date).days / 7
    on["weeks_on_dose"] = weeks_on_dose

    # Observed surplus (what we're predicting)
    observed = on["calories"].values
    tdee = on["tdee"].values
    sp_pressure = on["sp_pressure"].values
    blood = on["blood_level"].values
    wod = on["weeks_on_dose"].values

    # ═══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("1. FIT: appetite_per_unit and tachy_hl")
    print("=" * 70)

    def predict_intake(params, return_components=False):
        appetite_per_unit, tachy_hl = params
        if tachy_hl <= 0:
            return np.full(len(observed), 1e10)
        tachy_factor = np.exp(-np.log(2) / tachy_hl * wod)
        effective = blood * tachy_factor
        drug_effect = appetite_per_unit * effective
        predicted = tdee + sp_pressure + drug_effect
        if return_components:
            return predicted, sp_pressure, drug_effect, effective
        return predicted

    def rmse(params):
        pred = predict_intake(params)
        return np.sqrt(np.mean((observed - pred) ** 2))

    # Grid search then refine
    best_rmse = np.inf
    best_params = (-50, 30)
    for app in range(-120, -10, 5):
        for thl in range(8, 100, 4):
            r = rmse([app, thl])
            if r < best_rmse:
                best_rmse = r
                best_params = (app, thl)

    # Refine with continuous optimization
    result = minimize(lambda p: rmse(p), best_params, method="Nelder-Mead",
                      options={"maxiter": 10000, "xatol": 0.1})
    app_fit, thl_fit = result.x
    rmse_fit = result.fun

    print(f"\n  Appetite effect: {app_fit:.1f} cal per unit effective level")
    print(f"  Tachyphylaxis HL: {thl_fit:.1f} weeks")
    print(f"  RMSE: {rmse_fit:.0f} cal/day")

    # Compare with F's original estimates
    print(f"\n  F's original: -49 cal/unit, 32-week tachy HL")
    print(f"  Reconciled:   {app_fit:.0f} cal/unit, {thl_fit:.0f}-week tachy HL")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("2. COMPONENT DECOMPOSITION — what drives intake on each day?")
    print("=" * 70)

    pred, sp_p, drug_eff, eff_level = predict_intake([app_fit, thl_fit], return_components=True)
    on["predicted"] = pred
    on["sp_component"] = sp_p
    on["drug_component"] = drug_eff
    on["eff_level_fit"] = eff_level
    on["residual"] = observed - pred

    # By injection day
    print(f"\n  {'Day':>5} {'Observed':>9} {'Predicted':>10} {'SP press':>9} {'Drug eff':>9} {'Residual':>9}")
    for day in range(7):
        mask = on["days_since_injection"] == day
        if mask.sum() < 10:
            continue
        sub = on[mask]
        print(f"  {day:>5} {sub['calories'].mean():>9.0f} {sub['predicted'].mean():>10.0f} "
              f"{sub['sp_component'].mean():>+9.0f} {sub['drug_component'].mean():>+9.0f} "
              f"{sub['residual'].mean():>+9.0f}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("3. TEMPORAL VALIDATION — does the model track intake over time?")
    print("=" * 70)

    print(f"\n  {'Period':>12} {'Observed':>9} {'Predicted':>10} {'SP':>6} {'Drug':>7} {'Resid':>7} {'n':>4}")
    on["weeks_on_drug"] = (on["date"] - on["date"].min()).dt.days / 7
    for wk_start in range(0, 80, 8):
        wk_end = wk_start + 8
        mask = (on["weeks_on_drug"] >= wk_start) & (on["weeks_on_drug"] < wk_end)
        sub = on[mask]
        if len(sub) < 10:
            continue
        print(f"  {wk_start:>3}-{wk_end:<3} wk {sub['calories'].mean():>9.0f} {sub['predicted'].mean():>10.0f} "
              f"{sub['sp_component'].mean():>+6.0f} {sub['drug_component'].mean():>+7.0f} "
              f"{sub['residual'].mean():>+7.0f} {len(sub):>4}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("4. TROUGH vs PEAK TRENDS — does the model explain the +12 cal/week?")
    print("=" * 70)

    stable = on[on["dose_mg"] == 12.5].copy()
    stable["weeks_stable"] = (stable["date"] - stable["date"].min()).dt.days / 7
    trough = stable[stable["days_since_injection"].isin([4, 5, 6])]
    peak = stable[stable["days_since_injection"].isin([0, 1])]

    for label, sub in [("Trough (day 4-6)", trough), ("Peak (day 0-1)", peak)]:
        if len(sub) < 20:
            continue
        r_obs = np.corrcoef(sub["weeks_stable"], sub["calories"])[0, 1]
        r_pred = np.corrcoef(sub["weeks_stable"], sub["predicted"])[0, 1]
        r_resid = np.corrcoef(sub["weeks_stable"], sub["residual"])[0, 1]
        slope_obs = np.polyfit(sub["weeks_stable"], sub["calories"], 1)[0]
        slope_pred = np.polyfit(sub["weeks_stable"], sub["predicted"], 1)[0]
        slope_resid = np.polyfit(sub["weeks_stable"], sub["residual"], 1)[0]
        print(f"\n  {label}:")
        print(f"    Observed slope:  {slope_obs:+.1f} cal/week (r = {r_obs:+.3f})")
        print(f"    Predicted slope: {slope_pred:+.1f} cal/week (r = {r_pred:+.3f})")
        print(f"    Residual slope:  {slope_resid:+.1f} cal/week (r = {r_resid:+.3f})")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("5. SENSITIVITY — how do the fitted params change with SP HL?")
    print("=" * 70)

    print(f"\n  {'SP HL':>6} {'App/unit':>9} {'Tachy HL':>9} {'RMSE':>7}")
    for sp_hl_test in [30, 35, 40, 45, 50, 60, 80]:
        sp_test = ema(df["fat_mass_lbs"], sp_hl_test).values
        gap_test = sp_test[on.index] - on["fat_mass_lbs"].values
        sp_p_test = SP_PRESSURE * gap_test

        def rmse_test(params):
            app, thl = params
            if thl <= 0: return 1e10
            eff = blood * np.exp(-np.log(2) / thl * wod)
            pred = tdee + sp_p_test + app * eff
            return np.sqrt(np.mean((observed - pred) ** 2))

        res = minimize(lambda p: rmse_test(p), [-50, 30], method="Nelder-Mead")
        print(f"  {sp_hl_test:>4}d {res.x[0]:>+9.1f} {res.x[1]:>8.1f}wk {res.fun:>7.0f}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)
    print(f"\n  Set point (fixed, from pre-tirz 2014+):")
    print(f"    HL = {SP_HL}d symmetric, pressure = {SP_PRESSURE} cal/lb")
    print(f"  Drug parameters (fitted to on-drug data):")
    print(f"    Appetite: {app_fit:.0f} cal per unit effective level")
    print(f"    Tachyphylaxis HL: {thl_fit:.0f} weeks")
    print(f"    RMSE: {rmse_fit:.0f} cal/day (daily intake noise ~500 cal)")
    print(f"  The model decomposes each on-drug day into:")
    print(f"    TDEE (from Kalman) + SP pressure ({SP_PRESSURE} × gap) + drug ({app_fit:.0f} × eff_level)")


if __name__ == "__main__":
    main()
