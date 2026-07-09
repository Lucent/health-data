"""CZ3. Kalman tuning with step-aware TDEE anchor.

The hybrid filter:

    TDEE_t state drifts as random walk with σ_TDEE,
    but each step is gently pulled toward an external prediction:

        TDEE_external_t = c · expected_rmr_t + d · steps_t

    via an observation update with noise R_anchor.

Tightening σ_TDEE alone broke calorimetry plausibility (CZ2). The hope here:
the step-anchor explains some real TDEE variation, so the filter can
afford a tighter σ without breaking calibration.

Sweep:
  - σ_TDEE ∈ {5, 10, 15, 25, 50}
  - d ∈ {0.00, 0.02, 0.04, 0.06}
  - R_anchor ∈ {50, 200}
  - c held at 1.15 (best from CY)

For each config, bisect intake multiplier so EB closes.
Validate: FM/LM MAE, calorimetry mean and SD, EB residual, TDEE σ_dt.

Goal: find a (σ, d, R) that beats CZ2's σ=50 baseline on tightness without
losing FM tracking or calorimetry plausibility.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analysis"))

import P4_kalman_filter as P4
from P4_kalman_filter import (
    forbes_fat_fraction, build_inputs, load_data,
    E_FAT, E_LEAN, BUF_DECAY, BUF_CAP,
    ALPERT_RATE, DFM_GAIN_CAP, DLM_LOSS_CAP, DLM_GAIN_CAP,
    R_WEIGHT, R_FM_SCAN, R_LM_SCAN, R_BOOT_FM, R_BOOT_LM,
    Q_FM, Q_LM, Q_BUF,
)


def load_steps_arr(dates):
    steps_df = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])
    steps_map = steps_df.set_index("date")["steps"].to_dict()
    return np.array([float(steps_map.get(pd.Timestamp(d), 0.0)) for d in dates])


def kalman_step_anchored(calories, weight_obs, fm_scan, lm_scan, expected_rmr,
                          workouts_per_wk, tdee_external,
                          fm_start, lm_start,
                          sigma_tdee, max_dtdee, R_anchor):
    n = len(calories)
    q_tdee = sigma_tdee ** 2
    first_scan = np.where(~np.isnan(fm_scan))[0][0]
    init_tdee = (tdee_external[first_scan] if not np.isnan(tdee_external[first_scan])
                 else 2200.0)
    x = np.array([fm_start, lm_start, 0.0, init_tdee])
    P = np.diag([4.0, 4.0, 1.0, 200.0 ** 2])
    x_filt = np.full((n, 4), np.nan)
    P_filt = np.full((n, 4, 4), np.nan)
    x_pred_arr = np.full((n, 4), np.nan)
    P_pred_arr = np.full((n, 4, 4), np.nan)
    innovations = np.full(n, np.nan)
    x_filt[first_scan] = x
    P_filt[first_scan] = P
    x_pred_arr[first_scan] = x
    P_pred_arr[first_scan] = P
    Q = np.diag([Q_FM, Q_LM, Q_BUF, q_tdee])

    for i in range(first_scan, n - 1):
        fm, lm, buf, tdee = x
        ff = forbes_fat_fraction(fm, workouts_per_wk[i])
        deficit = (calories[i] if not np.isnan(calories[i]) else 0.0) - tdee
        e_per_lb = ff * E_FAT + (1 - ff) * E_LEAN
        dW_raw = deficit / e_per_lb
        dFM_raw = ff * dW_raw
        dLM_raw = (1 - ff) * dW_raw
        dFM = np.clip(dFM_raw, -fm * ALPERT_RATE, DFM_GAIN_CAP)
        dLM = np.clip(dLM_raw, -DLM_LOSS_CAP, DLM_GAIN_CAP)
        x_pred = np.array([
            max(fm + dFM, 5.0),
            max(lm + dLM, 50.0),
            BUF_DECAY * buf,
            tdee,
        ])
        F = np.array([
            [1.0, 0.0, 0.0, -ff / e_per_lb],
            [0.0, 1.0, 0.0, -(1 - ff) / e_per_lb],
            [0.0, 0.0, BUF_DECAY, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
        P_pred = F @ P @ F.T + Q
        x_pred_arr[i + 1] = x_pred
        P_pred_arr[i + 1] = P_pred
        x = x_pred.copy()
        P = P_pred.copy()

        # TDEE anchor pseudo-observation
        if not np.isnan(tdee_external[i + 1]):
            H = np.array([[0.0, 0.0, 0.0, 1.0]])
            y = tdee_external[i + 1]
            S = (H @ P @ H.T + R_anchor)[0, 0]
            K = (P @ H.T) / S
            x = x + (K * (y - (H @ x)[0])).flatten()
            P = P - np.outer(K.flatten(), H @ P)

        # Weight obs
        if not np.isnan(weight_obs[i + 1]):
            H = np.array([[1.0, 1.0, 1.0, 0.0]])
            y = weight_obs[i + 1]
            innovation = y - (H @ x)[0]
            S = (H @ P @ H.T + R_WEIGHT)[0, 0]
            K = (P @ H.T) / S
            x = x + (K * innovation).flatten()
            P = P - np.outer(K.flatten(), H @ P)
            innovations[i + 1] = innovation

        if not np.isnan(fm_scan[i + 1]):
            H = np.array([[1.0, 0.0, 0.0, 0.0]])
            y = fm_scan[i + 1]
            r_use = R_BOOT_FM if (i + 1 - first_scan) < 30 else R_FM_SCAN
            S = (H @ P @ H.T + r_use)[0, 0]
            K = (P @ H.T) / S
            x = x + (K * (y - (H @ x)[0])).flatten()
            P = P - np.outer(K.flatten(), H @ P)

        if not np.isnan(lm_scan[i + 1]):
            H = np.array([[0.0, 1.0, 0.0, 0.0]])
            y = lm_scan[i + 1]
            r_use = R_BOOT_LM if (i + 1 - first_scan) < 30 else R_LM_SCAN
            S = (H @ P @ H.T + r_use)[0, 0]
            K = (P @ H.T) / S
            x = x + (K * (y - (H @ x)[0])).flatten()
            P = P - np.outer(K.flatten(), H @ P)

        x[0] = max(x[0], 5.0)
        x[1] = max(x[1], 50.0)
        x[2] = np.clip(x[2], -BUF_CAP, BUF_CAP)
        prev_tdee = x_filt[i, 3] if not np.isnan(x_filt[i, 3]) else tdee
        x[3] = prev_tdee + np.clip(x[3] - prev_tdee, -max_dtdee, max_dtdee)

        x_filt[i + 1] = x
        P_filt[i + 1] = P

    return x_filt, P_filt, x_pred_arr, P_pred_arr, innovations, first_scan


def make_anchor(expected_rmr, steps_arr, c, d):
    rmr_safe = np.where(np.isnan(expected_rmr), 1900.0, expected_rmr)
    return c * rmr_safe + d * steps_arr


def run_config(intake_mult, c, d, sigma, cap, R_anchor):
    """Run one config and return SMOOTHED outputs."""
    P4.INTAKE_MULTIPLIER = intake_mult
    smooth, intake_daily, strength, comp_scans, rmr_df, p2, temp_df = load_data()
    dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr, workouts_per_wk = \
        build_inputs(smooth, intake_daily, strength, comp_scans, p2)
    steps_arr = load_steps_arr(dates)
    fm_start = fm_scan[~np.isnan(fm_scan)][0]
    lm_start = lm_scan[~np.isnan(lm_scan)][0]
    if R_anchor > 1e9:
        # Pure free-drift baseline: use P4 directly
        x_filt, P_filt, x_pred, P_pred, innov, fs = P4.kalman_v4_forward(
            calories, weight_obs, fm_scan, lm_scan, expected_rmr, workouts_per_wk,
            fm_start, lm_start, sigma, cap)
    else:
        ext = make_anchor(expected_rmr, steps_arr, c, d)
        x_filt, P_filt, x_pred, P_pred, innov, fs = kalman_step_anchored(
            calories, weight_obs, fm_scan, lm_scan, expected_rmr, workouts_per_wk,
            ext, fm_start, lm_start, sigma, cap, R_anchor)
    x_smooth, _ = P4.rts_smoother_v4(x_filt, P_filt, x_pred, P_pred, fs, workouts_per_wk)
    return dates, calories, fm_scan, lm_scan, x_smooth, rmr_df


def eb_residual(calories, x_smooth):
    valid = ~np.isnan(calories) & ~np.isnan(x_smooth[:, 3])
    cum_def = np.cumsum((calories - x_smooth[:, 3])[valid])
    fm_first = x_smooth[np.where(valid)[0][0], 0]
    fm_last = x_smooth[np.where(valid)[0][-1], 0]
    return float(cum_def[-1] / E_FAT - (fm_last - fm_first))


def find_eb_closing_m(c, d, sigma, cap, R_anchor, m_init=1.09, tol=0.5):
    def eb_at(m):
        _, calories, _, _, x_filt, _ = run_config(m, c, d, sigma, cap, R_anchor)
        return eb_residual(calories, x_filt)
    eb1 = eb_at(m_init)
    if abs(eb1) < tol:
        return m_init, eb1
    m2 = m_init - 0.02 * np.sign(eb1)
    eb2 = eb_at(m2)
    if abs(eb2 - eb1) < 1e-6:
        return m_init, eb1
    m_zero = m_init + (m2 - m_init) * (-eb1) / (eb2 - eb1)
    eb3 = eb_at(m_zero)
    if abs(eb3) > tol and abs(eb3 - eb1) > 1e-6:
        m_zero = m_init + (m_zero - m_init) * (-eb1) / (eb3 - eb1)
        eb3 = eb_at(m_zero)
    return m_zero, eb3


def metrics(label, dates, calories, fm_scan, lm_scan, x_smooth, rmr_df, m, c, d, sigma, R):
    out = {"label": label, "c": c, "d": d, "sigma": sigma, "R": R, "m": m}
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
    for cd, r in zip(rmr_dates, rmr_kcal):
        idx = np.searchsorted(dates, np.datetime64(cd))
        if idx < len(x_smooth) and not np.isnan(x_smooth[idx, 3]):
            ratios.append(x_smooth[idx, 3] / r)
    out["cal_mean"] = float(np.mean(ratios)) if ratios else float("nan")
    out["cal_sd"] = float(np.std(ratios)) if ratios else float("nan")
    out["cal_in_wide"] = float(np.mean([1.05 <= r <= 1.40 for r in ratios])) if ratios else 0.0

    out["eb_residual_lb"] = eb_residual(calories, x_smooth)

    tdee_valid = x_smooth[~np.isnan(x_smooth[:, 3]), 3]
    out["tdee_dt_sd"] = float(np.std(np.diff(tdee_valid))) if len(tdee_valid) > 1 else float("nan")
    out["tdee_mean"] = float(np.mean(tdee_valid))
    out["tdee_total_sd"] = float(np.std(tdee_valid))
    return out


def main():
    print("=" * 110)
    print("CZ3. Step-aware Kalman tuning. Anchor TDEE = c·RMR + d·steps")
    print("=" * 110)

    # Configs to test. (c, d, σ, cap, R_anchor)
    configs = [
        # Reference: free-drift baselines (R = ∞ effectively → use very large R)
        ("BASELINE σ=50 no anchor", 1.0, 0.0, 50, 100, 1e10),
        ("σ=25 no anchor",          1.0, 0.0, 25, 50,  1e10),
        ("σ=10 no anchor",          1.0, 0.0, 10, 20,  1e10),
        # Step-anchored at various tightness
        ("σ=25 c=1.15 d=0.04 R=200", 1.15, 0.04, 25, 50, 200),
        ("σ=15 c=1.15 d=0.04 R=200", 1.15, 0.04, 15, 25, 200),
        ("σ=10 c=1.15 d=0.04 R=200", 1.15, 0.04, 10, 20, 200),
        ("σ=10 c=1.15 d=0.04 R=50",  1.15, 0.04, 10, 20, 50),
        ("σ=5  c=1.15 d=0.04 R=50",  1.15, 0.04, 5,  10, 50),
        # Different d values at fixed σ=10, R=100
        ("σ=10 c=1.15 d=0.00 R=100", 1.15, 0.00, 10, 20, 100),
        ("σ=10 c=1.15 d=0.02 R=100", 1.15, 0.02, 10, 20, 100),
        ("σ=10 c=1.15 d=0.04 R=100", 1.15, 0.04, 10, 20, 100),
        ("σ=10 c=1.15 d=0.06 R=100", 1.15, 0.06, 10, 20, 100),
        # Tighter σ with stronger d
        ("σ=5  c=1.15 d=0.06 R=100", 1.15, 0.06, 5,  10, 100),
        ("σ=3  c=1.15 d=0.06 R=50",  1.15, 0.06, 3,  7,  50),
    ]

    results = []
    for label, c, d, sigma, cap, R in configs:
        try:
            m, eb_after = find_eb_closing_m(c, d, sigma, cap, R, m_init=1.09)
            dates, calories, fm_scan, lm_scan, x_filt, rmr_df = run_config(
                m, c, d, sigma, cap, R)
            v = metrics(label, dates, calories, fm_scan, lm_scan, x_filt, rmr_df,
                        m, c, d, sigma, R)
            results.append(v)
            print(f"  {label:<32}  m={m:.4f}  FM={v['fm_mae']:.4f}  LM={v['lm_mae']:.4f}  "
                  f"cal={v['cal_mean']:.3f}±{v['cal_sd']:.3f} wide={v['cal_in_wide']*100:.0f}%  "
                  f"EB={v['eb_residual_lb']:+.2f}  σ_dt={v['tdee_dt_sd']:.2f}")
        except Exception as ex:
            print(f"  {label}: FAILED {ex}")

    print("\n" + "=" * 130)
    print("PARETO TABLE — sorted by σ_dt ascending")
    print("=" * 130)
    print(f"{'config':<32} {'m':>7} {'FM MAE':>7} {'LM MAE':>7} "
          f"{'cal mean':>9} {'cal SD':>7} {'wide':>6} {'EB':>7} {'σ_dt':>6}")
    for r in sorted(results, key=lambda r: r["tdee_dt_sd"]):
        print(f"{r['label']:<32} {r['m']:>7.4f} "
              f"{r['fm_mae']:>7.4f} {r['lm_mae']:>7.4f} "
              f"{r['cal_mean']:>9.3f} {r['cal_sd']:>7.3f} {r['cal_in_wide']*100:>5.0f}% "
              f"{r['eb_residual_lb']:>+7.2f} {r['tdee_dt_sd']:>6.2f}")

    # ── Best step-anchored filter that beats baseline on tightness ──
    print("\n" + "=" * 110)
    print("BEST STEP-ANCHORED FILTER (lower σ_dt than baseline AND FM MAE ≤ baseline + 0.10)")
    print("=" * 110)
    baseline = next(r for r in results if "BASELINE" in r["label"])
    bl_fm = baseline["fm_mae"]
    bl_dt = baseline["tdee_dt_sd"]
    bl_cal = baseline["cal_mean"]
    bl_wide = baseline["cal_in_wide"]
    print(f"Baseline: FM={bl_fm:.4f}, σ_dt={bl_dt:.2f}, cal_mean={bl_cal:.3f}, "
          f"in_wide={bl_wide*100:.0f}%")
    candidates = [r for r in results
                  if r["tdee_dt_sd"] < bl_dt
                  and r["fm_mae"] <= bl_fm + 0.10
                  and abs(r["eb_residual_lb"]) < 5.0]
    if candidates:
        best = min(candidates, key=lambda r: r["tdee_dt_sd"])
        print(f"\nBest: {best['label']}")
        print(f"  m = {best['m']:.4f}, FM MAE = {best['fm_mae']:.4f} (baseline {bl_fm:.4f})")
        print(f"  σ_dt = {best['tdee_dt_sd']:.2f} (baseline {bl_dt:.2f}), "
              f"reduction = {(bl_dt-best['tdee_dt_sd'])/bl_dt*100:.0f}%")
        print(f"  cal_mean = {best['cal_mean']:.3f}, in_wide = {best['cal_in_wide']*100:.0f}%")
    else:
        print("\nNo step-anchored filter beats baseline tightness with FM MAE ≤ +0.10")


if __name__ == "__main__":
    main()
