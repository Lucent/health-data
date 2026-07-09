"""CY. Constrain TDEE drift; see if Kalman still locks on; add activity if needed.

Step 1 — Sweep TDEE process noise σ_TDEE and rate cap. Compare FM-scan MAE
   to the baseline filter (σ_TDEE=50, max_dTDEE=100). If FM MAE stays low
   under tight constraint, free TDEE drift was bookkeeping, not signal.

Step 2 — Rebuild the filter so TDEE is DETERMINISTIC each day:

       TDEE_t = c · expected_rmr_t  +  d · steps_t  +  e · workout_indicator_t

   Sweep / fit (c, d, e). Check whether explicit activity captures the TDEE
   variation that the free-drift filter was using.

Step 3 — Under the tightest filter that still tracks FM, re-run K hysteresis
   test. If K-style fall>rise TDEE asymmetry survives, it's genuine; if it
   collapses, it was a Kalman artifact of free-drift attribution.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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


def load_steps_workouts():
    """Return per-date arrays of steps and workout indicator aligned with intake dates."""
    steps = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])
    workouts = pd.read_csv(ROOT / "workout" / "strength.csv", parse_dates=["date"])
    return steps, workouts


def kalman_constrained(calories, weight_obs, fm_scan, lm_scan, expected_rmr,
                        workouts_per_wk, steps_arr, workout_ind_arr,
                        fm_start, lm_start,
                        c_rmr, d_steps, e_workout,
                        sigma_tdee=2.0, max_dtdee=5.0):
    """Modified forward Kalman where TDEE is deterministic each day:

       TDEE_t = c_rmr · expected_rmr_t + d_steps · steps_t + e_workout · workout_t

    Random-walk noise on TDEE is set tiny; rate cap is tight. The deterministic
    function dominates."""
    n = len(calories)
    q_tdee = sigma_tdee ** 2

    first_scan = np.where(~np.isnan(fm_scan))[0][0]

    def predicted_tdee(i):
        rmr_i = expected_rmr[i] if not np.isnan(expected_rmr[i]) else 1900.0
        steps_i = steps_arr[i] if not np.isnan(steps_arr[i]) else 0.0
        wo_i = workout_ind_arr[i] if not np.isnan(workout_ind_arr[i]) else 0.0
        return c_rmr * rmr_i + d_steps * steps_i + e_workout * wo_i

    init_tdee = predicted_tdee(first_scan)

    x = np.array([fm_start, lm_start, 0.0, init_tdee])
    P = np.diag([4.0, 4.0, 1.0, 50.0 ** 2])

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

        # OVERRIDE TDEE with deterministic prediction (forces activity-driven TDEE)
        tdee_target = predicted_tdee(i + 1)
        # Pull TDEE strongly toward target (within rate cap)
        tdee_used = tdee_target

        deficit = (calories[i] if not np.isnan(calories[i]) else 0.0) - tdee_used
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
            tdee_target,  # deterministic
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
        # Force TDEE to deterministic prediction (ignore drift)
        x[3] = predicted_tdee(i + 1)

        x_filt[i + 1] = x
        P_filt[i + 1] = P

    return x_filt, P_filt, x_pred_arr, P_pred_arr, innovations, first_scan


def fm_mae(x_smooth_or_filt, fm_scan):
    """Mean absolute error of FM state at scan dates."""
    err = []
    for i in range(len(fm_scan)):
        if not np.isnan(fm_scan[i]) and not np.isnan(x_smooth_or_filt[i, 0]):
            err.append(abs(x_smooth_or_filt[i, 0] - fm_scan[i]))
    return float(np.mean(err)) if err else float("nan")


def main():
    print("=" * 90)
    print("CY. Constrained-TDEE Kalman test")
    print("=" * 90)

    smooth, intake_daily, strength, comp_scans, rmr_df, p2, temp_df = load_data()
    steps_df, workouts_df = load_steps_workouts()

    dates, calories, weight_obs, fm_scan, lm_scan, expected_rmr, workouts_per_wk = \
        build_inputs(smooth, intake_daily, strength, comp_scans, p2)

    # Build steps_arr aligned with dates
    steps_map = steps_df.set_index("date")["steps"].to_dict()
    steps_arr = np.array([float(steps_map.get(pd.Timestamp(d), 0.0))
                           for d in dates])
    # workout indicator (1 if a strength session that day)
    workout_set = set(pd.to_datetime(workouts_df["date"]).dt.normalize())
    workout_ind_arr = np.array([1.0 if pd.Timestamp(d).normalize() in workout_set else 0.0
                                 for d in dates])

    print(f"\nDays: {len(dates)}")
    print(f"Steps mean = {np.nanmean(steps_arr[steps_arr > 0]):.0f}, "
          f"median = {np.nanmedian(steps_arr[steps_arr > 0]):.0f}")
    print(f"Workout days: {int(workout_ind_arr.sum())}")

    fm_start = fm_scan[~np.isnan(fm_scan)][0]
    lm_start = lm_scan[~np.isnan(lm_scan)][0]

    # ── Baseline P4 ──
    print("\nBaseline P4 (free-drift TDEE σ=50, max_dTDEE=100):")
    x_filt0, P_filt0, x_pred0, P_pred0, innov0, fs0 = P4.kalman_v4_forward(
        calories, weight_obs, fm_scan, lm_scan, expected_rmr,
        workouts_per_wk, fm_start, lm_start, 50, 100)
    x_smooth0, _ = P4.rts_smoother_v4(x_filt0, P_filt0, x_pred0, P_pred0, fs0,
                                        workouts_per_wk)
    mae_baseline = fm_mae(x_smooth0, fm_scan)
    tdee_mean_baseline = float(np.nanmean(x_smooth0[:, 3]))
    tdee_std_baseline = float(np.nanstd(x_smooth0[:, 3]))
    print(f"  FM scan MAE: {mae_baseline:.4f} lb")
    print(f"  TDEE mean:   {tdee_mean_baseline:.0f} cal/day, SD = {tdee_std_baseline:.0f}")

    # ── Step 1: Sweep σ_TDEE under free-drift architecture ──
    print("\n" + "=" * 90)
    print("STEP 1 — sweep σ_TDEE / max_dTDEE under existing free-drift architecture")
    print("=" * 90)
    print(f"  {'σ_TDEE':>7} {'max_dTDEE':>10} {'FM MAE':>9} {'TDEE SD':>9}")
    for sig, cap in [(50, 100), (30, 50), (15, 25), (5, 10), (2, 5), (0.5, 2)]:
        x_filt, P_filt, x_pred, P_pred, innov, fs = P4.kalman_v4_forward(
            calories, weight_obs, fm_scan, lm_scan, expected_rmr,
            workouts_per_wk, fm_start, lm_start, sig, cap)
        x_smooth, _ = P4.rts_smoother_v4(x_filt, P_filt, x_pred, P_pred, fs,
                                          workouts_per_wk)
        mae = fm_mae(x_smooth, fm_scan)
        tdee_sd = float(np.nanstd(x_smooth[:, 3]))
        print(f"  {sig:>7.1f} {cap:>10.1f} {mae:>9.4f} {tdee_sd:>9.0f}")

    # ── Step 2: Constrained Kalman with deterministic TDEE = c·RMR + d·steps + e·workouts ──
    print("\n" + "=" * 90)
    print("STEP 2 — sweep (c, d, e) for deterministic TDEE")
    print("=" * 90)
    # Coarse grid
    print(f"  {'c (RMR)':>8} {'d (steps)':>10} {'e (wo)':>8} {'FM MAE':>9} {'TDEE mean':>10} {'TDEE SD':>9}")
    best = {"mae": float("inf"), "c": None, "d": None, "e": None}
    grid_c = [1.10, 1.15, 1.20, 1.25, 1.30]
    grid_d = [0.0, 0.02, 0.04, 0.06]
    grid_e = [0.0, 50, 100, 150]
    for c in grid_c:
        for d in grid_d:
            for e in grid_e:
                try:
                    x_filt, P_filt, x_pred, P_pred, innov, fs = kalman_constrained(
                        calories, weight_obs, fm_scan, lm_scan, expected_rmr,
                        workouts_per_wk, steps_arr, workout_ind_arr,
                        fm_start, lm_start, c, d, e)
                    mae = fm_mae(x_filt, fm_scan)
                    tdee_mean = float(np.nanmean(x_filt[:, 3]))
                    tdee_sd = float(np.nanstd(x_filt[:, 3]))
                    print(f"  {c:>8.2f} {d:>10.3f} {e:>8.0f} {mae:>9.4f} "
                          f"{tdee_mean:>10.0f} {tdee_sd:>9.0f}")
                    if mae < best["mae"]:
                        best.update({"mae": mae, "c": c, "d": d, "e": e})
                except Exception as ex:
                    print(f"  {c:>8.2f} {d:>10.3f} {e:>8.0f}  FAILED: {ex}")

    print(f"\nBest constrained fit:")
    print(f"  c = {best['c']:.2f}, d = {best['d']:.3f}, e = {best['e']:.0f}")
    print(f"  FM MAE = {best['mae']:.4f} lb (baseline {mae_baseline:.4f})")
    print(f"  ΔMAE vs baseline: {best['mae'] - mae_baseline:+.4f} lb")

    # ── Step 3: Re-test K hysteresis under best constrained model ──
    if best["c"] is not None:
        print("\n" + "=" * 90)
        print("STEP 3 — K hysteresis under constrained TDEE")
        print("=" * 90)
        x_filt, P_filt, x_pred, P_pred, innov, fs = kalman_constrained(
            calories, weight_obs, fm_scan, lm_scan, expected_rmr,
            workouts_per_wk, steps_arr, workout_ind_arr,
            fm_start, lm_start, best["c"], best["d"], best["e"])

        date_idx = pd.to_datetime(dates)
        df = pd.DataFrame({
            "date": date_idx,
            "fm": x_filt[:, 0],
            "lm": x_filt[:, 1],
            "tdee": x_filt[:, 3],
            "expected_rmr": expected_rmr,
        })
        # Only look at pre-tirz for the cycle test
        tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
        df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
        df["effective_level"] = df["effective_level"].fillna(0)
        pre = df[df["effective_level"] == 0].copy()

        cycles = pd.read_csv(ROOT / "analysis" / "CT_cycles.csv",
                              parse_dates=["peak_a_date", "trough_date", "peak_b_date"])

        diffs = []
        for _, c in cycles.iterrows():
            mfall = (pre["date"] >= c["peak_a_date"]) & (pre["date"] <= c["trough_date"])
            mrise = (pre["date"] >= c["trough_date"]) & (pre["date"] <= c["peak_b_date"])
            gf, gr = pre[mfall], pre[mrise]
            if len(gf) > 7 and len(gr) > 7:
                fall_ratio = float((gf["tdee"] / gf["expected_rmr"]).mean())
                rise_ratio = float((gr["tdee"] / gr["expected_rmr"]).mean())
                diffs.append(fall_ratio - rise_ratio)
        diffs = np.array(diffs)
        n_pos = int(np.sum(diffs > 0))
        print(f"\n  Constrained K hysteresis test (CX-style):")
        print(f"  Cycles with fall TDEE/RMR > rise: {n_pos}/{len(diffs)} "
              f"({100*n_pos/len(diffs):.0f}%)")
        print(f"  Median Δ ratio: {np.median(diffs):+.4f}")
        print(f"  Compare to baseline P4: 14/19 (74%), median +0.0324")

        if n_pos / len(diffs) > 0.65 or abs(np.median(diffs)) > 0.02:
            print("  → K signal SURVIVES under constrained TDEE — possibly real")
        else:
            print("  → K signal COLLAPSES under constrained TDEE — was bookkeeping artifact")


if __name__ == "__main__":
    main()
