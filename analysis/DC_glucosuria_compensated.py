"""DC. Glucosuria retest with calorie compensation.

Why this exists (vs CG): CG clipped sugar above a threshold and subtracted
those calories from logged intake, then ran the Kalman filter with the
canonical INTAKE_MULTIPLIER = 1.09. But 1.09 was calibrated to close
cumulative EB residual to ~0 on the *unclipped* intake — so any clipping
adds to that 9 % undercount and the model gets driven away from the data
purely from magnitude loss. The MAE/EB metrics then drifted in a way that
told us nothing about whether the *shape* of the clipping was right.

Fix: after each clipping scheme, rescale the clipped intake so that
cumulative consumed calories match the baseline. Concretely, replace the
fixed 1.09 multiplier with

    m_clip = 1.09 × Σ(intake_original) / Σ(intake_clipped)

so that `m_clip × intake_clipped` integrates to the same total over the
record as `1.09 × intake_original`. Now the Kalman comparison reflects only
the *timing* benefit of clipping: do the over-threshold sugar spikes line
up with composition residuals such that removing them and redistributing
the calories tightens the fit?

If clipping still improves FM_MAE / innov_acf under this compensation, the
glucosuria signal is timing-real. If it doesn't, the previous CG drift was
just magnitude loss and there's nothing to find.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analysis"))

from P4_kalman_filter import (
    kalman_v4_forward, rts_smoother_v4, diagnostics, build_inputs,
    INTAKE_MULTIPLIER as M_BASE,
)

OUT = ROOT / "analysis" / "DC_glucosuria_compensated.csv"


def clipped_intake(clip_mode, threshold, intake_daily, foods):
    """Return (intake_with_clip_applied, total_discarded_cal)."""
    out = intake_daily.copy()
    if clip_mode == "none":
        return out, 0.0
    if clip_mode == "day":
        sugars_day = foods.groupby("date").sugars_g.sum()
        excess = (sugars_day - threshold).clip(lower=0)
        discarded = (excess * 4).rename("discarded_cal")
    elif clip_mode == "meal":
        meal_sugar = foods.groupby(["date", "meal"]).sugars_g.sum().reset_index()
        meal_sugar["excess"] = (meal_sugar["sugars_g"] - threshold).clip(lower=0)
        discarded = meal_sugar.groupby("date")["excess"].sum() * 4
        discarded.name = "discarded_cal"
    else:
        raise ValueError(clip_mode)
    out = out.set_index("date")
    out["calories"] = out["calories"].subtract(discarded, fill_value=0).clip(lower=0)
    out = out.reset_index()
    return out, float(discarded.sum())


def run_one(label, intake_clipped, total_discarded, total_original,
            smooth, strength, comp_scans, p2, rmr_df, temp_df):
    """Build inputs, override the multiplier via pre-scaling, run Kalman, report."""
    # Compensation factor: scale clipped calories so 1.09 × scaled equals
    # 1.09 × original on the 15-year sum.
    if total_discarded > 0:
        compensation = total_original / (total_original - total_discarded)
        intake_for_kalman = intake_clipped.copy()
        intake_for_kalman["calories"] = intake_for_kalman["calories"] * compensation
    else:
        compensation = 1.0
        intake_for_kalman = intake_clipped

    dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr, workouts_per_wk = (
        build_inputs(smooth, intake_for_kalman, strength, comp_scans, p2)
    )
    first_scan = np.where(~np.isnan(fm_scan))[0][0]
    fm_start = fm_scan[first_scan]; lm_start = lm_scan[first_scan]
    x_filt, P_filt, x_pred, P_pred, innov, _ = kalman_v4_forward(
        calories, weight_obs, fm_scan, lm_scan, expected_rmr,
        workouts_per_wk, fm_start, lm_start, 50, 100,
    )
    x_smooth, _ = rts_smoother_v4(x_filt, P_filt, x_pred, P_pred, first_scan, workouts_per_wk)
    diag = diagnostics(dates, x_smooth, innov, fm_scan, lm_scan, calories, rmr_df,
                       weight_obs, temp_df)

    # Verify EB closure improved (or didn't drift away from 0)
    return dict(
        label=label,
        compensation=compensation,
        m_effective=M_BASE * compensation,
        discarded_total=total_discarded,
        fm_mae=diag["fm_scan_mae"],
        lm_mae=diag["lm_scan_mae"],
        innov_acf=diag["innov_acf1"],
        innov_std=diag["innov_std"],
        eb_residual=diag["eb_residual_lb"],
        tdee_rmr_mean=diag["tdee_rmr_mean"],
        r_temp=diag.get("tdee_temp_r_partial_fm", float("nan")),
    )


def main():
    smooth = pd.read_csv(ROOT / "analysis" / "P1_smoothed_weight.csv", parse_dates=["date"])
    intake_daily = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    foods = pd.read_csv(ROOT / "intake" / "intake_foods.csv", parse_dates=["date"])
    strength = pd.read_csv(ROOT / "workout" / "strength.csv", parse_dates=["date"])
    comp_scans = pd.read_csv(ROOT / "composition" / "composition.csv", parse_dates=["date"])
    rmr_df = pd.read_csv(ROOT / "RMR" / "rmr.csv", parse_dates=["date"])
    p2 = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    temp_df = pd.read_csv(ROOT / "analysis" / "AS_temperature_baseline.csv", parse_dates=["date"])

    total_original_cal = float(intake_daily["calories"].sum())

    rows = []

    # Baseline
    rows.append(run_one("baseline", intake_daily, 0.0, total_original_cal,
                         smooth, strength, comp_scans, p2, rmr_df, temp_df))

    print("=" * 110)
    print("DC. Glucosuria retest with compensated intake multiplier")
    print(f"    Cumulative logged intake (15 yr): {total_original_cal:,.0f} cal")
    print(f"    Baseline multiplier: {M_BASE}  →  effective intake = {M_BASE*total_original_cal:,.0f} cal")
    print("=" * 110)

    fmt_hdr = (f"{'label':>10} {'discarded':>11} {'m_eff':>7} "
               f"{'FM_MAE':>7} {'LM_MAE':>7} {'innov_acf':>10} "
               f"{'innov_std':>10} {'EB(lb)':>8} {'TDEE/RMR':>9} {'r(T)|FM':>8}")
    print(fmt_hdr)

    def emit(r):
        print(f"  {r['label']:>8} {r['discarded_total']:>11,.0f} "
              f"{r['m_effective']:>7.4f} "
              f"{r['fm_mae']:>7.3f} {r['lm_mae']:>7.3f} "
              f"{r['innov_acf']:>+10.3f} {r['innov_std']:>10.3f} "
              f"{r['eb_residual']:>+8.2f} {r['tdee_rmr_mean']:>9.3f} "
              f"{r['r_temp']:>+8.3f}")
    emit(rows[0])

    print("\n--- Per-day clipping (compensated) ---")
    for T in [50, 75, 100, 125, 150, 200]:
        clipped, disc = clipped_intake("day", T, intake_daily, foods)
        rows.append(run_one(f"d{T}", clipped, disc, total_original_cal,
                            smooth, strength, comp_scans, p2, rmr_df, temp_df))
        emit(rows[-1])

    print("\n--- Per-meal clipping (compensated) ---")
    for T in [25, 40, 50, 75, 100]:
        clipped, disc = clipped_intake("meal", T, intake_daily, foods)
        rows.append(run_one(f"m{T}", clipped, disc, total_original_cal,
                            smooth, strength, comp_scans, p2, rmr_df, temp_df))
        emit(rows[-1])

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print(f"\nArtifact: {OUT.relative_to(ROOT)}")

    # Verdict
    base = df.iloc[0]
    improved = df[(df.fm_mae < base.fm_mae) & (df.label != "baseline")]
    if len(improved):
        best = improved.sort_values("fm_mae").iloc[0]
        print(f"\nBest clip vs baseline: {best.label}  "
              f"FM_MAE {base.fm_mae:.4f} → {best.fm_mae:.4f}  "
              f"(Δ={best.fm_mae-base.fm_mae:+.4f} lb)")
        print(f"  This is timing-only since compensation preserves total calories.")
    else:
        print("\nNo clip beats baseline FM_MAE under compensation → "
              "no glucosuria timing signal even with magnitude controlled.")


if __name__ == "__main__":
    main()
