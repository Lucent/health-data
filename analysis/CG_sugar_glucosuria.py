"""CG. Does the subject urinate glucose on high-sugar days?

Renal threshold ≈ 180 mg/dL blood glucose. The SURMOUNT-2 HbA1c-stratified
analysis attributes T2D's attenuated weight loss largely to glucosuria —
high HbA1c → urinary glucose loss → passive caloric drain. This subject is
non-diabetic (HbA1c ~5.0%), so glucosuria should be rare. But very large
rapid-sugar meals could spike BG above threshold transiently.

Test: clip sugar above a threshold (per day or per meal), subtract the
clipped excess as urinated calories from intake, re-run the Kalman filter,
and compare fit metrics. If clipping improves the filter — tighter scan
MAE, lower innovation ACF, better energy-balance closure — the subject
genuinely urinates glucose above that threshold and the model treats those
calories as reaching tissue when they shouldn't.

If clipping makes the fit worse or makes no difference: the calories did
reach tissue; no glucosuria; energy balance is unchanged. Either way, the
test is decisive against the post hoc residual regression that confused
gut-content noise with energy loss.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analysis"))

# Import Kalman pieces directly
from P4_kalman_filter import (
    kalman_v4_forward, rts_smoother_v4, diagnostics, build_inputs,
    SIGMA_TDEE_GRID, MAX_DTDEE_GRID,
)


def load_inputs(clip_mode, threshold):
    """Build all Kalman inputs with optional sugar clipping.

    clip_mode in {"none", "day", "meal"}.
    threshold in grams of sugar; clip_mode="day" clips per-day total, "meal" per meal item.
    Returns: dates, calories_clipped, weight_obs, fm_scan, lm_scan, expected_rmr,
             training_inc, workouts_per_wk, daily_clipped_cal_total
    """
    smooth = pd.read_csv(ROOT / "analysis" / "P1_smoothed_weight.csv", parse_dates=["date"])
    intake_daily = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    foods = pd.read_csv(ROOT / "intake" / "intake_foods.csv", parse_dates=["date"])
    strength = pd.read_csv(ROOT / "workout" / "strength.csv", parse_dates=["date"])
    comp_scans = pd.read_csv(ROOT / "composition" / "composition.csv", parse_dates=["date"])
    rmr_df = pd.read_csv(ROOT / "RMR" / "rmr.csv", parse_dates=["date"])
    p2 = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])

    intake_clipped = intake_daily.copy()
    discarded_per_day = pd.Series(0.0, index=intake_daily["date"])
    if clip_mode == "day":
        # All sources used (no source filter). Sum sugars per day, clip excess above T.
        sugars_day = foods.groupby("date").sugars_g.sum()
        excess = (sugars_day - threshold).clip(lower=0)
        discarded = (excess * 4).rename("discarded_cal")
        intake_clipped = intake_clipped.set_index("date")
        intake_clipped["calories"] = intake_clipped["calories"].subtract(discarded, fill_value=0).clip(lower=0)
        intake_clipped = intake_clipped.reset_index()
        discarded_per_day = discarded
    elif clip_mode == "meal":
        meal_sugar = foods.groupby(["date", "meal"]).sugars_g.sum().reset_index()
        meal_sugar["excess"] = (meal_sugar["sugars_g"] - threshold).clip(lower=0)
        discarded = meal_sugar.groupby("date")["excess"].sum() * 4
        discarded.name = "discarded_cal"
        intake_clipped = intake_clipped.set_index("date")
        intake_clipped["calories"] = intake_clipped["calories"].subtract(discarded, fill_value=0).clip(lower=0)
        intake_clipped = intake_clipped.reset_index()
        discarded_per_day = discarded.reindex(intake_daily["date"], fill_value=0)

    dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr, workouts_per_wk = \
        build_inputs(smooth, intake_clipped, strength, comp_scans, p2)

    return (dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr,
            workouts_per_wk, rmr_df, discarded_per_day)


def run_kalman_with_intake(inputs, sigma_tdee=50, max_dtdee=100, temp_df=None):
    (dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr,
     workouts_per_wk, rmr_df, _) = inputs
    first_scan = np.where(~np.isnan(fm_scan))[0][0]
    fm_start = fm_scan[first_scan]
    lm_start = lm_scan[first_scan]
    x_filt, P_filt, x_pred, P_pred, innov, _ = kalman_v4_forward(
        calories, weight_obs, fm_scan, lm_scan, expected_rmr,
        workouts_per_wk, fm_start, lm_start, sigma_tdee, max_dtdee,
    )
    x_smooth, _ = rts_smoother_v4(x_filt, P_filt, x_pred, P_pred, first_scan, workouts_per_wk)
    diag = diagnostics(dates, x_smooth, innov, fm_scan, lm_scan, calories, rmr_df,
                       weight_obs, temp_df)
    return diag, x_smooth, calories


def per_year_eb_residual(dates, calories, x_smooth):
    """Per-year cumulative EB residual: implied dFM (from intake-TDEE) minus observed dFM.
    Negative = model over-predicts loss; positive = model under-predicts loss.
    """
    df = pd.DataFrame({
        "date": pd.to_datetime(dates),
        "cal": calories,
        "tdee": x_smooth[:, 3],
        "fm": x_smooth[:, 0],
    }).dropna(subset=["cal", "tdee", "fm"])
    df["year"] = df.date.dt.year
    out = {}
    for y, g in df.groupby("year"):
        cum_def = (g.cal - g.tdee).sum()
        fm_change = g.fm.iloc[-1] - g.fm.iloc[0]
        implied = cum_def / 3500.0
        out[y] = implied - fm_change
    return out


def main():
    temp_df = pd.read_csv(ROOT / "analysis" / "AS_temperature_baseline.csv", parse_dates=["date"])

    # Baseline: no clipping
    print("=" * 100)
    print("CG. Sugar-glucosuria test: clip sugar above T, re-run Kalman, compare fit")
    print("    All food-log sources used (no MHTML/HTML filter).")
    print("=" * 100)
    print("\nBaseline (no clipping):")
    base_inputs = load_inputs("none", 0)
    base_diag, base_x, base_cal = run_kalman_with_intake(base_inputs, temp_df=temp_df)
    base_dates = base_inputs[0]
    base_yearly = per_year_eb_residual(base_dates, base_cal, base_x)
    fmt = ("FM_MAE={fm_scan_mae:.3f}  LM_MAE={lm_scan_mae:.3f}  "
           "innov_acf={innov_acf1:+.3f}  innov_std={innov_std:.3f}  "
           "TDEE/RMR mean={tdee_rmr_mean:.3f}  EB={eb_residual_lb:+.2f}  "
           "r(T)|FM={tdee_temp_r_partial_fm:+.3f}")
    print(f"  {fmt.format(**base_diag)}")

    # Per-day threshold sweep
    print("\n--- Per-day clipping ---")
    print(f"{'T_day':>6} {'discarded_total':>16} {'days_clipped':>13} {'FM_MAE':>7} "
          f"{'LM_MAE':>7} {'innov_acf':>10} {'innov_std':>10} {'EB(lb)':>8} {'r(T)|FM':>8}")
    print(f"  none {'0':>16} {'0':>13} "
          f"{base_diag['fm_scan_mae']:>7.3f} {base_diag['lm_scan_mae']:>7.3f} "
          f"{base_diag['innov_acf1']:>+10.3f} {base_diag['innov_std']:>10.3f} "
          f"{base_diag['eb_residual_lb']:>+8.2f} "
          f"{base_diag.get('tdee_temp_r_partial_fm', float('nan')):>+8.3f}")
    yearly_runs = {}
    for T in [50, 75, 100, 125, 150, 200]:
        ins = load_inputs("day", T)
        discarded = ins[-1]
        d, x_s, cal = run_kalman_with_intake(ins, temp_df=temp_df)
        n_clipped = (discarded > 0).sum()
        total = discarded.sum()
        print(f"  {T:>4} {total:>16,.0f} {n_clipped:>13} "
              f"{d['fm_scan_mae']:>7.3f} {d['lm_scan_mae']:>7.3f} "
              f"{d['innov_acf1']:>+10.3f} {d['innov_std']:>10.3f} "
              f"{d['eb_residual_lb']:>+8.2f} "
              f"{d.get('tdee_temp_r_partial_fm', float('nan')):>+8.3f}")
        yearly_runs[("day", T)] = per_year_eb_residual(ins[0], cal, x_s)

    # Per-meal threshold sweep
    print("\n--- Per-meal clipping ---")
    print(f"{'T_meal':>6} {'discarded_total':>16} {'days_clipped':>13} {'FM_MAE':>7} "
          f"{'LM_MAE':>7} {'innov_acf':>10} {'innov_std':>10} {'EB(lb)':>8} {'r(T)|FM':>8}")
    for T in [25, 40, 50, 75, 100]:
        ins = load_inputs("meal", T)
        discarded = ins[-1]
        d, x_s, cal = run_kalman_with_intake(ins, temp_df=temp_df)
        n_clipped = (discarded > 0).sum()
        total = discarded.sum()
        print(f"  {T:>4} {total:>16,.0f} {n_clipped:>13} "
              f"{d['fm_scan_mae']:>7.3f} {d['lm_scan_mae']:>7.3f} "
              f"{d['innov_acf1']:>+10.3f} {d['innov_std']:>10.3f} "
              f"{d['eb_residual_lb']:>+8.2f} "
              f"{d.get('tdee_temp_r_partial_fm', float('nan')):>+8.3f}")
        yearly_runs[("meal", T)] = per_year_eb_residual(ins[0], cal, x_s)

    # Year-by-year residual under each clipping
    print(f"\n--- Per-year EB residual (lb): baseline vs clipping (negative = model over-predicts loss) ---")
    years = sorted(base_yearly.keys())
    cols = [("baseline", base_yearly)]
    for T in [50, 75, 100, 125, 150, 200]:
        cols.append((f"d{T}", yearly_runs[("day", T)]))
    for T in [25, 40, 50, 75, 100]:
        cols.append((f"m{T}", yearly_runs[("meal", T)]))
    header = f"  {'year':>5}"
    for label, _ in cols:
        header += f" {label:>9}"
    header += f"  {'best_clip':>10}"
    print(header)
    for y in years:
        row = f"  {y:>5}"
        vals = []
        for label, dct in cols:
            v = dct.get(y, float("nan"))
            row += f" {v:>+9.2f}"
            vals.append((label, v))
        # Find smallest |residual|
        valid = [(l, v) for l, v in vals if not np.isnan(v)]
        best = min(valid, key=lambda kv: abs(kv[1]))
        row += f"  {best[0]:>10}"
        print(row)

    # Year-paired: 2-year rolling windows
    print(f"\n--- 2-year-pair EB residuals (sum across pair) ---")
    pairs = [(years[i], years[i + 1]) for i in range(len(years) - 1)]
    header = f"  {'pair':>10}"
    for label, _ in cols:
        header += f" {label:>9}"
    header += f"  {'best_clip':>10}"
    print(header)
    for y1, y2 in pairs:
        row = f"  {y1}-{y2:>4}"
        vals = []
        for label, dct in cols:
            v = dct.get(y1, 0) + dct.get(y2, 0)
            row += f" {v:>+9.2f}"
            vals.append((label, v))
        best = min(vals, key=lambda kv: abs(kv[1]))
        row += f"  {best[0]:>10}"
        print(row)

    print("\n=== Verdict criteria ===")
    print("  If glucosuria is real: FM_MAE drops, innov_acf drops, EB closure tightens")
    print("  at some non-trivial T (e.g., T_day=100 or T_meal=50). Improvements should")
    print("  be small (sugar > threshold is rare for this subject) but consistent.")
    print()
    print("  If clipping makes fit worse or unchanged: subject does NOT meaningfully")
    print("  urinate glucose at any threshold tested. Energy-balance is intact.")


if __name__ == "__main__":
    main()
