"""CZ2. Kalman tuning, with co-fit intake multiplier and Pareto curve.

For each candidate σ_TDEE / max_dTDEE, refit the intake multiplier so that
cumulative EB closes (otherwise tighter filters look broken just from
miscalibration). Then report:

  FM scan MAE          — tracks composition scans
  LM scan MAE          — tracks LM scans
  cal ratio mean / SD  — TDEE / measured_RMR at the 25 calorimetry days
  TDEE σ_dt            — day-to-day TDEE volatility (the tightness measure)
  intake multiplier    — what multiplier makes EB close
  EB residual          — should be ~0 after refit

We're looking for the Pareto frontier: lowest TDEE σ_dt subject to FM MAE
not blowing up. The user picks the trade-off.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analysis"))

import P4_kalman_filter as P4
from P4_kalman_filter import (
    build_inputs, load_data, kalman_v4_forward, rts_smoother_v4, E_FAT,
)


def run_filter_with_multiplier(intake_multiplier, sig, cap):
    """Run baseline P4 with a specified intake multiplier and TDEE constraint."""
    P4.INTAKE_MULTIPLIER = intake_multiplier
    smooth, intake_daily, strength, comp_scans, rmr_df, p2, temp_df = load_data()
    dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr, workouts_per_wk = \
        build_inputs(smooth, intake_daily, strength, comp_scans, p2)
    fm_start = fm_scan[~np.isnan(fm_scan)][0]
    lm_start = lm_scan[~np.isnan(lm_scan)][0]
    x_filt, P_filt, x_pred, P_pred, innov, fs = kalman_v4_forward(
        calories, weight_obs, fm_scan, lm_scan, expected_rmr,
        workouts_per_wk, fm_start, lm_start, sig, cap)
    x_smooth, _ = rts_smoother_v4(x_filt, P_filt, x_pred, P_pred, fs, workouts_per_wk)
    return dates, calories, fm_scan, lm_scan, x_smooth, rmr_df


def find_eb_closing_multiplier(sig, cap, m_init=1.09, tol=0.5):
    """Bisect intake multiplier so EB residual ≈ 0."""
    def eb_at(m):
        dates, calories, fm_scan, lm_scan, x_smooth, rmr_df = run_filter_with_multiplier(
            m, sig, cap)
        valid = ~np.isnan(calories) & ~np.isnan(x_smooth[:, 3])
        cum_def = np.cumsum((calories - x_smooth[:, 3])[valid])
        fm_first = x_smooth[np.where(valid)[0][0], 0]
        fm_last = x_smooth[np.where(valid)[0][-1], 0]
        return float(cum_def[-1] / E_FAT - (fm_last - fm_first))

    eb1 = eb_at(m_init)
    if abs(eb1) < tol:
        return m_init, eb1
    m2 = m_init - 0.01 * np.sign(eb1)
    eb2 = eb_at(m2)
    if abs(eb2 - eb1) < 1e-6:
        return m_init, eb1
    m_zero = m_init + (m2 - m_init) * (-eb1) / (eb2 - eb1)
    eb3 = eb_at(m_zero)
    if abs(eb3) > tol:
        # one more iteration
        if abs(eb3 - eb1) > 1e-6:
            m_zero = m_init + (m_zero - m_init) * (-eb1) / (eb3 - eb1)
            eb3 = eb_at(m_zero)
    return m_zero, eb3


def metrics(label, dates, calories, fm_scan, lm_scan, x_smooth, rmr_df, m_used):
    out = {"label": label, "intake_mult": m_used}
    fm_errs = [abs(x_smooth[i, 0] - fm_scan[i])
                for i in range(len(fm_scan))
                if not np.isnan(fm_scan[i]) and not np.isnan(x_smooth[i, 0])]
    out["fm_mae"] = float(np.mean(fm_errs)) if fm_errs else float("nan")

    lm_errs = [abs(x_smooth[i, 1] - lm_scan[i])
                for i in range(len(lm_scan))
                if not np.isnan(lm_scan[i]) and not np.isnan(x_smooth[i, 1])]
    out["lm_mae"] = float(np.mean(lm_errs)) if lm_errs else float("nan")

    rmr_dates = pd.to_datetime(rmr_df["date"]).values
    rmr_kcal = rmr_df["rmr_kcal"].values
    ratios = []
    for d, r in zip(rmr_dates, rmr_kcal):
        idx = np.searchsorted(dates, np.datetime64(d))
        if idx < len(x_smooth) and not np.isnan(x_smooth[idx, 3]):
            ratios.append(x_smooth[idx, 3] / r)
    out["cal_n"] = len(ratios)
    out["cal_mean"] = float(np.mean(ratios)) if ratios else float("nan")
    out["cal_sd"] = float(np.std(ratios)) if ratios else float("nan")
    out["cal_in_wide_band"] = float(np.mean([1.05 <= r <= 1.40 for r in ratios])) if ratios else 0.0

    valid = ~np.isnan(calories) & ~np.isnan(x_smooth[:, 3])
    cum_def = np.cumsum((calories - x_smooth[:, 3])[valid])
    fm_first = x_smooth[np.where(valid)[0][0], 0]
    fm_last = x_smooth[np.where(valid)[0][-1], 0]
    out["eb_residual_lb"] = float(cum_def[-1] / E_FAT - (fm_last - fm_first))

    tdee_valid = x_smooth[~np.isnan(x_smooth[:, 3]), 3]
    out["tdee_dt_sd"] = float(np.std(np.diff(tdee_valid))) if len(tdee_valid) > 1 else float("nan")
    out["tdee_total_sd"] = float(np.std(tdee_valid))
    out["tdee_mean"] = float(np.mean(tdee_valid))
    return out


def main():
    print("=" * 110)
    print("CZ2. Kalman tuning — co-fit intake multiplier per σ_TDEE setting")
    print("=" * 110)

    configs = [
        (50, 100, 1.09),  # baseline anchor
        (35, 70, None),
        (25, 50, None),
        (15, 25, None),
        (10, 20, None),
        (5, 10, None),
        (3, 7, None),
        (2, 5, None),
        (1, 3, None),
    ]
    results = []
    for sig, cap, m_init in configs:
        if m_init is None:
            m, eb_after = find_eb_closing_multiplier(sig, cap, m_init=1.09)
        else:
            m = m_init
            _, eb_after = find_eb_closing_multiplier(sig, cap, m_init=m_init, tol=0.1)
        dates, calories, fm_scan, lm_scan, x_smooth, rmr_df = run_filter_with_multiplier(
            m, sig, cap)
        v = metrics(f"σ={sig:>4} cap={cap:>4}", dates, calories, fm_scan, lm_scan,
                    x_smooth, rmr_df, m)
        v["sig"] = sig
        v["cap"] = cap
        results.append(v)
        print(f"  σ={sig:>4} cap={cap:>4}  m={m:.4f}  "
              f"FM_MAE={v['fm_mae']:.4f}  LM_MAE={v['lm_mae']:.4f}  "
              f"cal_mean={v['cal_mean']:.3f}±{v['cal_sd']:.3f} "
              f"in_wide={v['cal_in_wide_band']*100:.0f}%  "
              f"EB={v['eb_residual_lb']:+.2f}  σ_dt={v['tdee_dt_sd']:.2f}")

    print("\n" + "=" * 110)
    print("PARETO TABLE — sorted by tightness (TDEE σ_dt) ascending")
    print("=" * 110)
    print(f"{'σ_TDEE':>7} {'cap':>5} {'m':>7} "
          f"{'FM MAE':>7} {'LM MAE':>7} "
          f"{'cal mean':>9} {'cal SD':>7} {'wide':>6} "
          f"{'EB lb':>7} {'σ_dt':>6} {'σ_total':>8}")
    sorted_r = sorted(results, key=lambda r: r["tdee_dt_sd"])
    for r in sorted_r:
        print(f"{r['sig']:>7d} {r['cap']:>5d} {r['intake_mult']:>7.4f} "
              f"{r['fm_mae']:>7.4f} {r['lm_mae']:>7.4f} "
              f"{r['cal_mean']:>9.3f} {r['cal_sd']:>7.3f} {r['cal_in_wide_band']*100:>5.0f}% "
              f"{r['eb_residual_lb']:>+7.2f} {r['tdee_dt_sd']:>6.2f} {r['tdee_total_sd']:>8.0f}")

    print("\n" + "=" * 110)
    print("ASSESSMENT")
    print("=" * 110)
    baseline = next(r for r in results if r["sig"] == 50)
    bl_fm = baseline["fm_mae"]
    bl_lm = baseline["lm_mae"]
    print(f"Baseline (σ=50): FM_MAE = {bl_fm:.4f}, LM_MAE = {bl_lm:.4f}, "
          f"σ_dt = {baseline['tdee_dt_sd']:.2f}")
    print(f"\nFor each tighter setting, tradeoff vs baseline:")
    print(f"{'σ_TDEE':>7} {'ΔFM':>7} {'ΔLM':>7} {'σ_dt':>6} {'σ_dt reduction':>16}")
    for r in sorted(results, key=lambda r: -r["sig"]):
        if r["sig"] == 50:
            print(f"{r['sig']:>7d} {0:>+7.4f} {0:>+7.4f} {r['tdee_dt_sd']:>6.2f}  "
                  f"baseline")
        else:
            print(f"{r['sig']:>7d} {r['fm_mae']-bl_fm:>+7.4f} {r['lm_mae']-bl_lm:>+7.4f} "
                  f"{r['tdee_dt_sd']:>6.2f}  "
                  f"{(baseline['tdee_dt_sd'] - r['tdee_dt_sd'])/baseline['tdee_dt_sd']*100:>+5.0f}%")


if __name__ == "__main__":
    main()
