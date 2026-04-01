#!/usr/bin/env python3
"""BQ. Adaptive Freeze Latch: reconcile subject and trial latch behavior.

Motivation:
  The corrected full-window hold models fit the subject better.
  The short-lookback near-freeze models fit withdrawal trials better.

Treat those not as contradictions, but as clues:
  - adaptation should increase with stability / residence
  - adaptation should be suppressed by recent movement ("freeze")
  - downward adaptation toward lower FM may be harder than upward adaptation

Shared biology:
  pressure = 55 * (SP - FM)

Shared SP update:
  SP_{t+1} = SP_t + rate_t * (target_t - SP_t)

where:
  target_t = trailing mean FM over a residence window
  rate_t = base_rate * stability_score * freeze_release_score * direction_factor

This script searches that rule on:
  1. subject 2014+ pre-tirzepatide mean surplus fit
  2. trial-only withdrawal transfer
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PRESSURE_PER_LB = 55.0
KCAL_PER_LB_FAT = 3500.0
REGAIN_FAT_FRACTIONS = [1.00, 0.90, 0.80, 0.70]

TRIALS = [
    {"name": "SURMOUNT-4", "baseline_weight": 107.3, "stop_weight": 85.8, "target_regain_pct": 14.0, "treatment_days": 36 * 7, "regain_days": 52 * 7, "baseline_bf": 46.6 / 102.5, "loss_fat_fraction": 0.75},
    {"name": "STEP-1 extension", "baseline_weight": 105.6, "stop_weight": 87.5, "target_regain_pct": (99.0 - 87.5) / 87.5 * 100.0, "treatment_days": 68 * 7, "regain_days": 52 * 7, "baseline_bf": 0.434, "loss_fat_fraction": 0.56},
    {"name": "STEP-4", "baseline_weight": 107.2, "stop_weight": 96.1, "target_regain_pct": 6.9, "treatment_days": 20 * 7, "regain_days": 48 * 7, "baseline_bf": 0.434, "loss_fat_fraction": 0.56},
]


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def load_subject():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    df = intake[["date", "calories"]].merge(kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="inner")
    df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df = df.sort_values("date").reset_index(drop=True)

    fm = df["fat_mass_lbs"].values.copy()
    first_valid = np.where(~np.isnan(fm))[0][0]
    fm[:first_valid] = fm[first_valid]
    for i in range(first_valid + 1, len(fm)):
        if np.isnan(fm[i]):
            fm[i] = fm[i - 1]
    df["fm"] = fm
    df["surplus"] = df["calories"] - df["tdee"]
    df["surplus_90"] = df["surplus"].rolling(90, min_periods=90).mean()
    df["surplus_180"] = df["surplus"].rolling(180, min_periods=180).mean()
    mask = ((df["date"] >= "2014-01-01") & (df["date"] < "2024-01-01") & (df["effective_level"] == 0)).values
    return df, fm, mask


def trial_anchor(trial):
    baseline_fm = trial["baseline_weight"] * trial["baseline_bf"]
    lost_weight = trial["baseline_weight"] - trial["stop_weight"]
    stop_fm = baseline_fm - lost_weight * trial["loss_fat_fraction"]
    return baseline_fm, stop_fm


def sp_adaptive_freeze_series(fm, residence_days, tol, move_window, move_thresh, move_scale, base_rate, down_factor):
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        start = max(0, i - residence_days + 1)
        window = fm[start:i + 1]
        target = window.mean()
        stability = max(0.0, 1.0 - np.max(np.abs(window - target)) / tol)

        j = max(0, i - move_window)
        vel = abs((fm[i] - fm[j]) / max(1, i - j))
        freeze_release = sigmoid((move_thresh - vel) / move_scale)

        direction = down_factor if target < sp[i - 1] else 1.0
        rate = base_rate * stability * freeze_release * direction
        sp[i] = sp[i - 1] + rate * (target - sp[i - 1])
    return sp


def sp_adaptive_freeze_step(hist, sp_prev, residence_days, tol, move_window, move_thresh, move_scale, base_rate, down_factor):
    i = len(hist) - 1
    start = max(0, i - residence_days + 1)
    window = np.asarray(hist[start:i + 1])
    target = window.mean()
    stability = max(0.0, 1.0 - np.max(np.abs(window - target)) / tol)

    j = max(0, i - move_window)
    vel = abs((hist[i] - hist[j]) / max(1, i - j))
    freeze_release = sigmoid((move_thresh - vel) / move_scale)

    direction = down_factor if target < sp_prev else 1.0
    rate = base_rate * stability * freeze_release * direction
    return sp_prev + rate * (target - sp_prev)


def eval_subject(sp, fm, surplus90, surplus180, mask):
    dist = sp - fm
    out = {}
    for name, target in [("r90", surplus90), ("r180", surplus180)]:
        valid = mask & ~np.isnan(target) & ~np.isnan(dist)
        out[name] = np.corrcoef(dist[valid], target[valid])[0, 1] if valid.sum() > 200 else np.nan
    return out


def simulate_trial(params, trial, regain_fat_fraction):
    start_fm, stop_fm = trial_anchor(trial)
    fm_treat = np.linspace(start_fm, stop_fm, trial["treatment_days"])
    sp_treat = sp_adaptive_freeze_series(
        fm_treat,
        params["residence_days"],
        params["tol"],
        params["move_window"],
        params["move_thresh"],
        params["move_scale"],
        params["base_rate"],
        params["down_factor"],
    )
    fm_hist = list(fm_treat)
    fm = stop_fm
    sp = sp_treat[-1]
    for _ in range(1, trial["regain_days"]):
        gap = sp - fm
        fm = fm + PRESSURE_PER_LB * gap / KCAL_PER_LB_FAT
        fm_hist.append(fm)
        sp = sp_adaptive_freeze_step(
            fm_hist, sp,
            params["residence_days"], params["tol"],
            params["move_window"], params["move_thresh"], params["move_scale"],
            params["base_rate"], params["down_factor"],
        )
    regained_weight = (fm - stop_fm) / regain_fat_fraction
    return regained_weight / trial["stop_weight"] * 100.0


def eval_trials(params):
    totals = []
    mids = {}
    for frac in REGAIN_FAT_FRACTIONS:
        errs = []
        for trial in TRIALS:
            pred = simulate_trial(params, trial, frac)
            errs.append(abs(pred - trial["target_regain_pct"]))
            mids.setdefault(trial["name"], []).append(pred)
        totals.append(sum(errs))
    out = {
        "trial_mean_abs_err": float(np.mean(totals)),
        "trial_worst_abs_err": float(np.max(totals)),
    }
    for trial in TRIALS:
        out[f"{trial['name']}_mid"] = np.asarray(mids[trial["name"]])[1]
    return out


def main():
    df, fm, mask = load_subject()
    rows = []

    for residence_days in [14, 21, 28, 42]:
        for tol in [2.5, 3.0, 4.0]:
            for move_window in [7, 14, 21]:
                for move_thresh in [0.03, 0.05, 0.08, 0.10]:
                    for base_rate in [0.005, 0.0075, 0.01]:
                        for down_factor in [0.1, 0.2, 0.3, 0.5, 1.0]:
                            params = {
                                "residence_days": residence_days,
                                "tol": tol,
                                "move_window": move_window,
                                "move_thresh": move_thresh,
                                "move_scale": 0.02,
                                "base_rate": base_rate,
                                "down_factor": down_factor,
                            }
                            sp = sp_adaptive_freeze_series(
                                fm, residence_days, tol, move_window, move_thresh, 0.02, base_rate, down_factor
                            )
                            fit = eval_subject(sp, fm, df["surplus_90"].values, df["surplus_180"].values, mask)
                            trials = eval_trials(params)
                            score = abs(fit["r90"]) - 0.02 * trials["trial_mean_abs_err"] - 0.01 * trials["trial_worst_abs_err"]
                            rows.append({
                                **params,
                                "subject_r90": fit["r90"],
                                "subject_r180": fit["r180"],
                                "score": score,
                                **trials,
                            })

    out = pd.DataFrame(rows).sort_values(["score", "subject_r90"], ascending=False).reset_index(drop=True)
    artifact_path = ROOT / "analysis" / "BQ_adaptive_freeze_latch_search.csv"
    out.to_csv(artifact_path, index=False)

    print("=" * 100)
    print("ADAPTIVE FREEZE LATCH SEARCH")
    print("=" * 100)
    print(f"\nTotal models tested: {len(out)}")

    cols = [
        "residence_days", "tol", "move_window", "move_thresh", "base_rate", "down_factor",
        "subject_r90", "subject_r180",
        "trial_mean_abs_err", "trial_worst_abs_err",
        "SURMOUNT-4_mid", "STEP-1 extension_mid", "STEP-4_mid",
    ]

    print("\nTop 30 overall:")
    print(out[cols].head(30).to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    good = out[out["subject_r90"].abs() > 0.75].copy()
    print("\nClosest to trials among |subject_r90| > 0.75:")
    print(good.sort_values(["trial_mean_abs_err", "subject_r90"], ascending=[True, False])[cols]
          .head(30)
          .to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print(f"\nArtifact: {artifact_path}")


if __name__ == "__main__":
    main()
