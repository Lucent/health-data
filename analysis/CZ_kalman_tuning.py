"""CZ. Tune the Kalman filter as tight as possible while still tracking
FM, LM, and RMR properly.

Validation criteria (all must pass for a filter to be acceptable):
  1. FM scan MAE   — should not degrade much from baseline (~1.08 lb)
  2. LM scan MAE   — should not degrade much from baseline
  3. TDEE/RMR at calorimetry days — ratio should fall in 0.95-1.30 band (held-out)
  4. EB residual   — cumulative over 15 years, |residual| ≤ ~6 lb
  5. TDEE day-to-day SD — should be low (free-drift baseline ≈ 12 cal/day)

Filter variants tested:
  A. Baseline P4: σ_TDEE = 50, max_dTDEE = 100  (free drift, current canonical)
  B-F. Free-drift with tighter σ
  G. Deterministic: TDEE_t = c · expected_rmr_t  (no activity)
  H. Deterministic + steps: TDEE = c · RMR + d · steps
  I. Deterministic + steps + workouts
  J. Deterministic + steps + workouts + TEF
  K. Deterministic + steps + workouts + TEF, with small process noise added

For each variant: report all metrics. Choose the tightest filter that
satisfies all validation criteria.
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


def load_extra_inputs(dates):
    """Load steps and workout indicator arrays aligned with `dates`."""
    steps_df = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])
    workouts_df = pd.read_csv(ROOT / "workout" / "strength.csv", parse_dates=["date"])
    steps_map = steps_df.set_index("date")["steps"].to_dict()
    workout_set = set(pd.to_datetime(workouts_df["date"]).dt.normalize())
    steps_arr = np.array([float(steps_map.get(pd.Timestamp(d), 0.0)) for d in dates])
    workout_ind = np.array(
        [1.0 if pd.Timestamp(d).normalize() in workout_set else 0.0 for d in dates])
    return steps_arr, workout_ind


def kalman_with_tdee_anchor(calories, weight_obs, fm_scan, lm_scan, expected_rmr,
                              workouts_per_wk, tdee_external,
                              fm_start, lm_start,
                              sigma_tdee, max_dtdee, R_tdee_anchor):
    """Forward Kalman where TDEE is pulled toward `tdee_external[t]` at each
    step via an observation update with noise R_tdee_anchor.

    R_tdee_anchor very small (e.g. 1.0) => TDEE locked to external prediction.
    R_tdee_anchor very large => behaves like baseline P4 (free drift).
    """
    n = len(calories)
    q_tdee = sigma_tdee ** 2
    first_scan = np.where(~np.isnan(fm_scan))[0][0]

    init_tdee = (tdee_external[first_scan] if not np.isnan(tdee_external[first_scan])
                 else expected_rmr[first_scan] * 1.15 if not np.isnan(expected_rmr[first_scan])
                 else 2100.0)

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
            tdee,  # random walk
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

        # TDEE anchor observation (the new bit)
        if not np.isnan(tdee_external[i + 1]):
            H = np.array([[0.0, 0.0, 0.0, 1.0]])
            y = tdee_external[i + 1]
            S = (H @ P @ H.T + R_tdee_anchor)[0, 0]
            K = (P @ H.T) / S
            x = x + (K * (y - (H @ x)[0])).flatten()
            P = P - np.outer(K.flatten(), H @ P)

        # Weight obs
        if not np.isnan(weight_obs[i + 1]):
            H = np.array([[1.0, 1.0, 1.0, 0.0]])
            y = weight_obs[i + 1]
            y_pred = (H @ x)[0]
            innovation = y - y_pred
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


def validate(label, x_smooth, dates, fm_scan, lm_scan, calories, rmr_df):
    """Compute the validation metrics for a filter trajectory."""
    out = {"label": label}

    # 1. FM MAE
    fm_errs = [abs(x_smooth[i, 0] - fm_scan[i])
                for i in range(len(fm_scan))
                if not np.isnan(fm_scan[i]) and not np.isnan(x_smooth[i, 0])]
    out["fm_mae"] = float(np.mean(fm_errs)) if fm_errs else float("nan")

    # 2. LM MAE
    lm_errs = [abs(x_smooth[i, 1] - lm_scan[i])
                for i in range(len(lm_scan))
                if not np.isnan(lm_scan[i]) and not np.isnan(x_smooth[i, 1])]
    out["lm_mae"] = float(np.mean(lm_errs)) if lm_errs else float("nan")

    # 3. TDEE/RMR at calorimetry days
    rmr_dates = pd.to_datetime(rmr_df["date"]).values
    rmr_kcal = rmr_df["rmr_kcal"].values
    ratios = []
    for d, r in zip(rmr_dates, rmr_kcal):
        idx = np.searchsorted(dates, np.datetime64(d))
        if idx < len(x_smooth) and not np.isnan(x_smooth[idx, 3]):
            ratios.append(x_smooth[idx, 3] / r)
    out["cal_n"] = len(ratios)
    out["cal_ratio_mean"] = float(np.mean(ratios)) if ratios else float("nan")
    out["cal_ratio_sd"] = float(np.std(ratios)) if ratios else float("nan")
    out["cal_ratio_in_band"] = float(np.mean([0.95 <= r <= 1.30 for r in ratios])) if ratios else 0.0

    # 4. Cumulative EB
    valid_cal = ~np.isnan(calories) & ~np.isnan(x_smooth[:, 3])
    if valid_cal.sum() > 0:
        cum_deficit = np.cumsum((calories - x_smooth[:, 3])[valid_cal])
        fm_first = x_smooth[np.where(valid_cal)[0][0], 0]
        fm_last = x_smooth[np.where(valid_cal)[0][-1], 0]
        out["eb_residual_lb"] = float(cum_deficit[-1] / E_FAT - (fm_last - fm_first))
    else:
        out["eb_residual_lb"] = float("nan")

    # 5. TDEE smoothness
    tdee_valid = x_smooth[~np.isnan(x_smooth[:, 3]), 3]
    if len(tdee_valid) > 1:
        out["tdee_dt_sd"] = float(np.std(np.diff(tdee_valid)))
        out["tdee_mean"] = float(np.mean(tdee_valid))
        out["tdee_total_sd"] = float(np.std(tdee_valid))
    else:
        out["tdee_dt_sd"] = out["tdee_mean"] = out["tdee_total_sd"] = float("nan")

    return out


def main():
    print("=" * 100)
    print("CZ. Kalman tuning sweep — find tightest filter that still tracks")
    print("=" * 100)

    smooth, intake_daily, strength, comp_scans, rmr_df, p2, temp_df = load_data()
    dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr, workouts_per_wk = \
        build_inputs(smooth, intake_daily, strength, comp_scans, p2)
    steps_arr, workout_ind = load_extra_inputs(dates)

    fm_start = fm_scan[~np.isnan(fm_scan)][0]
    lm_start = lm_scan[~np.isnan(lm_scan)][0]

    results = []

    # ─── Free-drift sweep (varying σ_TDEE, no anchor) ───
    print("\nFamily 1: free-drift, varying σ_TDEE (no anchor)")
    for sig, cap in [(50, 100), (25, 50), (10, 20), (5, 10), (2, 5)]:
        x_filt, P_filt, x_pred, P_pred, innov, fs = P4.kalman_v4_forward(
            calories, weight_obs, fm_scan, lm_scan, expected_rmr,
            workouts_per_wk, fm_start, lm_start, sig, cap)
        x_smooth, _ = P4.rts_smoother_v4(x_filt, P_filt, x_pred, P_pred, fs, workouts_per_wk)
        v = validate(f"FreeDrift σ={sig}, cap={cap}", x_smooth, dates,
                      fm_scan, lm_scan, calories, rmr_df)
        results.append(v)

    # ─── Anchor-to-deterministic sweep ───
    # tdee_external_t = c·RMR + d·steps + e·workouts + f·TEF (for various combos)
    print("\nFamily 2: anchored to deterministic prediction")

    def make_external(c, d=0.0, e=0.0, f=0.0):
        rmr_safe = np.where(np.isnan(expected_rmr), 1900.0, expected_rmr)
        cal_safe = np.where(np.isnan(calories), 0.0, calories)
        return c * rmr_safe + d * steps_arr + e * workout_ind + f * cal_safe

    # Use σ_TDEE=2 and small max_dTDEE for tight anchoring; vary R_tdee_anchor
    # to allow some slack
    configs = [
        # (label, c, d, e, f, R_anchor, sig, cap)
        ("AnchRMR only c=1.20",          1.20, 0.0,  0.0,  0.0, 25,    5, 10),
        ("AnchRMR+steps c=1.15 d=0.04",  1.15, 0.04, 0.0,  0.0, 25,    5, 10),
        ("AnchRMR+steps c=1.18 d=0.04",  1.18, 0.04, 0.0,  0.0, 25,    5, 10),
        ("AnchRMR+steps+TEF c=1.15 d=0.04 f=0.10", 1.15, 0.04, 0.0, 0.10, 25, 5, 10),
        ("AnchRMR+steps+wo c=1.15 d=0.04 e=100",   1.15, 0.04, 100.0, 0.0, 25, 5, 10),
        ("LooseAnch c=1.18 d=0.04 R=200", 1.18, 0.04, 0.0, 0.0, 200,   25, 50),
        ("LooseAnch c=1.18 d=0.04 R=2500", 1.18, 0.04, 0.0, 0.0, 2500, 50, 100),
        ("TightAnch c=1.18 d=0.04 R=5",   1.18, 0.04, 0.0, 0.0, 5,     2, 5),
    ]
    for label, c, d, e, f, R_anchor, sig, cap in configs:
        ext = make_external(c, d, e, f)
        x_filt, P_filt, x_pred, P_pred, innov, fs = kalman_with_tdee_anchor(
            calories, weight_obs, fm_scan, lm_scan, expected_rmr,
            workouts_per_wk, ext,
            fm_start, lm_start, sig, cap, R_anchor)
        # No smoother for anchored (different F structure); use filtered
        x_smooth = x_filt
        v = validate(label, x_smooth, dates, fm_scan, lm_scan, calories, rmr_df)
        v["c"] = c
        v["d"] = d
        v["e"] = e
        v["f"] = f
        v["R_anchor"] = R_anchor
        v["sig"] = sig
        results.append(v)

    # ─── Report ───
    print("\n" + "=" * 130)
    print(f"{'Filter':<46} {'FM MAE':>7} {'LM MAE':>7} {'cal mean':>9} "
          f"{'cal SD':>7} {'in band':>8} {'EB lb':>7} {'TDEE σ_dt':>10} {'TDEE mean':>10}")
    print("=" * 130)
    for r in results:
        print(f"{r['label']:<46} "
              f"{r['fm_mae']:>7.4f} {r['lm_mae']:>7.4f} "
              f"{r['cal_ratio_mean']:>9.3f} {r['cal_ratio_sd']:>7.3f} "
              f"{r['cal_ratio_in_band']*100:>7.0f}% {r['eb_residual_lb']:>+7.2f} "
              f"{r['tdee_dt_sd']:>10.2f} {r['tdee_mean']:>10.0f}")

    # ─── Verdict ───
    print("\n" + "=" * 100)
    print("ACCEPTABILITY (FM MAE ≤ 1.20, LM MAE ≤ 2.0, cal in_band ≥ 0.80, |EB| ≤ 8 lb)")
    print("=" * 100)
    for r in results:
        ok_fm = r["fm_mae"] <= 1.20
        ok_lm = r["lm_mae"] <= 2.0
        ok_cal = r["cal_ratio_in_band"] >= 0.80
        ok_eb = abs(r["eb_residual_lb"]) <= 8.0 if not np.isnan(r["eb_residual_lb"]) else False
        passes = ok_fm and ok_lm and ok_cal and ok_eb
        marks = ("F" if not ok_fm else "·") + ("L" if not ok_lm else "·") + \
                ("C" if not ok_cal else "·") + ("E" if not ok_eb else "·")
        verdict = "PASS" if passes else f"FAIL [{marks}]"
        print(f"  {r['label']:<46} {verdict:<14} "
              f"σ_dt = {r['tdee_dt_sd']:>6.2f}")


if __name__ == "__main__":
    main()
