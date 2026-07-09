"""CI. Are the textbook Atwater factors (4 cal/g protein, 4 cal/g carb, 9 cal/g fat)
optimal for this subject's intake conversion, or would alternative values tighten
the Kalman fit?

Same clip-and-rerun protocol as CG, CH:

  alt_cal = logged_cal + (a−4)×protein_g + (b−4)×carbs_g + (c−9)×fat_g

The filler (calories from sources beyond P/C/F — alcohol, sugar alcohols) is
preserved as `logged_cal − 4P − 4C − 9F` and added back automatically. On days
where macros are missing or zero, alt_cal reverts to logged_cal.

Sweep (a, b, c) and re-run the Kalman. If a non-textbook combination produces
a meaningfully tighter fit, it suggests this subject's metabolisable energy per
gram of macronutrient differs from the population average — possibly from gut
absorption efficiency, food-class composition, or a systematic MFP-database bias.

Sanity: per-day Atwater residual (logged − 4P − 4C − 9F) has median −19 cal,
IQR ±60 cal — tight enough that 1 cal/g shifts on macros produce detectable
~50-100 cal/day adjustments on macro-rich days.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analysis"))

from P4_kalman_filter import (
    kalman_v4_forward, rts_smoother_v4, diagnostics, build_inputs, load_data,
)


def load_inputs_with_atwater(a_p, b_c, c_f):
    """Load Kalman inputs with intake recomputed at alternative Atwater factors."""
    smooth, intake_daily, strength, comp_scans, rmr_df, p2, temp_df = load_data()
    # Adjust calories: alt = logged + (a-4)P + (b-4)C + (c-9)F
    intake_alt = intake_daily.copy()
    adj = ((a_p - 4.0) * intake_alt.protein_g.fillna(0) +
           (b_c - 4.0) * intake_alt.carbs_g.fillna(0) +
           (c_f - 9.0) * intake_alt.fat_g.fillna(0))
    intake_alt["calories"] = (intake_alt["calories"] + adj).clip(lower=0)

    dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr, workouts_per_wk = \
        build_inputs(smooth, intake_alt, strength, comp_scans, p2)
    return (dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr,
            workouts_per_wk, rmr_df, temp_df, adj.sum())


def run(inputs):
    (dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr,
     workouts_per_wk, rmr_df, temp_df, _) = inputs
    first_scan = np.where(~np.isnan(fm_scan))[0][0]
    fm_start = fm_scan[first_scan]
    lm_start = lm_scan[first_scan]
    x_filt, P_filt, x_pred, P_pred, innov, _ = kalman_v4_forward(
        calories, weight_obs, fm_scan, lm_scan, expected_rmr,
        workouts_per_wk, fm_start, lm_start, 50, 100,
    )
    x_smooth, _ = rts_smoother_v4(x_filt, P_filt, x_pred, P_pred, first_scan, workouts_per_wk)
    return diagnostics(dates, x_smooth, innov, fm_scan, lm_scan, calories, rmr_df,
                        weight_obs, temp_df)


def fmt_row(label, d, total_adj):
    return (f"  {label:<32} "
            f"{total_adj:>+11,.0f} "
            f"{d['fm_scan_mae']:>7.3f} {d['lm_scan_mae']:>7.3f} "
            f"{d['innov_acf1']:>+10.3f} {d['innov_std']:>8.3f} "
            f"{d['eb_residual_lb']:>+8.2f} "
            f"{d.get('tdee_temp_r_partial_fm', float('nan')):>+8.3f}")


def main():
    print("=" * 100)
    print("CI. Atwater factors test: sweep (cal/g protein, carb, fat) and re-run Kalman.")
    print("    Filler (logged − 4P − 4C − 9F) preserved, so days with no macros are unaffected.")
    print("=" * 100)
    HEADER = (f"  {'config':<32} {'Σadj(cal)':>11} {'FM_MAE':>7} {'LM_MAE':>7} "
              f"{'innov_acf':>10} {'innov_std':>8} {'EB(lb)':>8} {'r(T)|FM':>8}")

    # Reference: textbook Atwater
    print("\n--- Reference: textbook Atwater (4, 4, 9) ---")
    print(HEADER)
    ins_ref = load_inputs_with_atwater(4.0, 4.0, 9.0)
    ref = run(ins_ref)
    print(fmt_row("textbook (4, 4, 9)", ref, 0))

    # 1D sweeps
    print("\n--- 1D sweep: protein (carb=4, fat=9) ---")
    print(HEADER)
    for a in [3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0]:
        ins = load_inputs_with_atwater(a, 4.0, 9.0)
        d = run(ins)
        print(fmt_row(f"P={a}", d, ins[-1]))

    print("\n--- 1D sweep: carbs (protein=4, fat=9) ---")
    print(HEADER)
    for b in [3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0]:
        ins = load_inputs_with_atwater(4.0, b, 9.0)
        d = run(ins)
        print(fmt_row(f"C={b}", d, ins[-1]))

    print("\n--- 1D sweep: fat (protein=4, carbs=4) ---")
    print(HEADER)
    for c in [7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0]:
        ins = load_inputs_with_atwater(4.0, 4.0, c)
        d = run(ins)
        print(fmt_row(f"F={c}", d, ins[-1]))

    print("\n=== Verdict criteria ===")
    print("  If textbook (4, 4, 9) is best: subject metabolises macros at population average.")
    print("  If a different combo wins on FM_MAE + EB closure jointly: subject has a")
    print("  systematic per-gram efficiency offset worth modelling. Then 2D refinement.")


if __name__ == "__main__":
    main()
