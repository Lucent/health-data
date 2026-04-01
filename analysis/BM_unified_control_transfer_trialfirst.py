#!/usr/bin/env python3
"""BM. Unified subject-control + trial-withdrawal model, scored trial-first.

Same shared lipostat law for both subject and trials:
  pressure = 55 * (SP - FM)

Subject only:
  exhaustible control stock opposes pressure

Trials:
  zero control

Unlike BL, this script searches the known stronger trial-regain candidates and
prioritizes trial transfer first, then subject-side coherence.
"""

from math import exp, log
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
    mask = (df["date"] >= "2011-04-21") & (df["date"] < "2017-01-01") & (df["effective_level"] == 0)
    return df.loc[mask].copy().reset_index(drop=True)


def trial_anchor(trial):
    baseline_fm = trial["baseline_weight"] * trial["baseline_bf"]
    lost_weight = trial["baseline_weight"] - trial["stop_weight"]
    stop_fm = baseline_fm - lost_weight * trial["loss_fat_fraction"]
    return baseline_fm, stop_fm


def sp_ema_series(fm, hl):
    alpha = 1 - exp(-log(2) / hl)
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        sp[i] = sp[i - 1] + alpha * (fm[i] - sp[i - 1])
    return sp


def sp_lookback_series(fm, tol, hold, rate):
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        if i >= hold and abs(fm[i] - fm[i - hold]) <= tol:
            sp[i] = sp[i - 1] + rate * (fm[i] - sp[i - 1])
        else:
            sp[i] = sp[i - 1]
    return sp


def sp_hold_mean_series(fm, tol, hold, rate):
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        if i >= hold:
            window = fm[i - hold + 1:i + 1]
            center = window.mean()
            if np.max(np.abs(window - center)) <= tol:
                sp[i] = sp[i - 1] + rate * (center - sp[i - 1])
                continue
        sp[i] = sp[i - 1]
    return sp


def sp_hold_current_series(fm, tol, hold, rate):
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        if i >= hold:
            window = fm[i - hold + 1:i + 1]
            if np.max(np.abs(window - fm[i])) <= tol:
                sp[i] = sp[i - 1] + rate * (fm[i] - sp[i - 1])
                continue
        sp[i] = sp[i - 1]
    return sp


def step_rule(kind, params, hist, sp_prev):
    if kind == "EMA":
        hl = params["hl"]
        alpha = 1 - exp(-log(2) / hl)
        return sp_prev + alpha * (hist[-1] - sp_prev)
    tol, hold, rate = params["tol"], params["hold"], params["rate"]
    i = len(hist) - 1
    if kind == "Lookback":
        if i >= hold and abs(hist[i] - hist[i - hold]) <= tol:
            return sp_prev + rate * (hist[i] - sp_prev)
        return sp_prev
    if i < hold:
        return sp_prev
    window = np.asarray(hist[i - hold + 1:i + 1])
    if kind == "HoldMean":
        center = window.mean()
        if np.max(np.abs(window - center)) <= tol:
            return sp_prev + rate * (center - sp_prev)
        return sp_prev
    if np.max(np.abs(window - hist[i])) <= tol:
        return sp_prev + rate * (hist[i] - sp_prev)
    return sp_prev


def infer_control_stock(required, cmax, recovery, depletion, rest_gain):
    n = len(required)
    stock = np.empty(n)
    exerted = np.empty(n)
    shortfall = np.empty(n)
    stock[0] = cmax
    exerted[0] = min(stock[0], required[0])
    shortfall[0] = required[0] - exerted[0]
    for i in range(1, n):
        recover = recovery * (cmax - stock[i - 1])
        if required[i - 1] < 100:
            recover += rest_gain
        stock_now = np.clip(stock[i - 1] + recover - depletion * exerted[i - 1], 0, cmax)
        stock[i] = stock_now
        exerted[i] = min(stock_now, required[i])
        shortfall[i] = required[i] - exerted[i]
    return stock, exerted, shortfall


def eval_subject_with_control(df, sp, baseline, cmax, recovery, depletion, rest_gain):
    pressure = PRESSURE_PER_LB * (sp - df["fm"].values)
    required = np.maximum(0.0, baseline + pressure - df["surplus"].values)
    stock, exerted, shortfall = infer_control_stock(required, cmax, recovery, depletion, rest_gain)
    explained = baseline + pressure - exerted
    residual = df["surplus"].values - explained
    return {
        "pressure": pressure,
        "required": required,
        "stock": stock,
        "exerted": exerted,
        "shortfall": shortfall,
        "rmse": float(np.sqrt(np.mean(residual ** 2))),
        "impossible_days": float((shortfall > 100).mean()),
        "mean_shortfall": float(shortfall.mean()),
        "total_control": float(exerted.sum()),
    }


def simulate_trial(kind, params, trial, regain_fat_fraction):
    start_fm, stop_fm = trial_anchor(trial)
    fm_treat = np.linspace(start_fm, stop_fm, trial["treatment_days"])
    if kind == "EMA":
        sp_treat = sp_ema_series(fm_treat, params["hl"])
    elif kind == "Lookback":
        sp_treat = sp_lookback_series(fm_treat, params["tol"], params["hold"], params["rate"])
    elif kind == "HoldMean":
        sp_treat = sp_hold_mean_series(fm_treat, params["tol"], params["hold"], params["rate"])
    else:
        sp_treat = sp_hold_current_series(fm_treat, params["tol"], params["hold"], params["rate"])
    fm_hist = list(fm_treat)
    fm = stop_fm
    sp = sp_treat[-1]
    for _ in range(1, trial["regain_days"]):
        gap = sp - fm
        fm = fm + PRESSURE_PER_LB * gap / KCAL_PER_LB_FAT
        fm_hist.append(fm)
        sp = step_rule(kind, params, fm_hist, sp)
    regained_weight = (fm - stop_fm) / regain_fat_fraction
    return regained_weight / trial["stop_weight"] * 100.0


def eval_trials(kind, params):
    totals = []
    mids = {}
    for frac in REGAIN_FAT_FRACTIONS:
        errs = []
        for trial in TRIALS:
            pred = simulate_trial(kind, params, trial, frac)
            errs.append(abs(pred - trial["target_regain_pct"]))
            mids.setdefault(trial["name"], []).append(pred)
        totals.append(sum(errs))
    out = {"trial_mean_abs_err": float(np.mean(totals)), "trial_worst_abs_err": float(np.max(totals))}
    for trial in TRIALS:
        arr = np.asarray(mids[trial["name"]])
        out[f"{trial['name']}_mid"] = arr[1]
    return out


def main():
    df = load_subject()
    model_specs = []
    for hl in [45, 50]:
        model_specs.append(("EMA", {"hl": hl}, f"hl={hl}", sp_ema_series(df["fm"].values, hl)))
    for tol in [1.5, 2.5, 3.0]:
        for hold in [7, 14]:
            for rate in [0.001]:
                model_specs.append(("Lookback", {"tol": tol, "hold": hold, "rate": rate}, f"tol={tol},hold={hold},rate={rate}", sp_lookback_series(df["fm"].values, tol, hold, rate)))
    for tol in [2.5, 3.0, 4.0]:
        for hold in [28, 42, 56]:
            for rate in [0.0075, 0.01]:
                model_specs.append(("HoldMean", {"tol": tol, "hold": hold, "rate": rate}, f"tol={tol},hold={hold},rate={rate}", sp_hold_mean_series(df["fm"].values, tol, hold, rate)))
    for tol in [6.0]:
        for hold in [14]:
            for rate in [0.0075]:
                model_specs.append(("HoldCurrent", {"tol": tol, "hold": hold, "rate": rate}, f"tol={tol},hold={hold},rate={rate}", sp_hold_current_series(df["fm"].values, tol, hold, rate)))

    rows = []
    for kind, params, label, sp in model_specs:
        trial_eval = eval_trials(kind, params)
        for baseline in [-400, -300]:
            for cmax in [900, 1200]:
                for recovery in [0.02]:
                    for depletion in [0.4, 0.7]:
                        for rest_gain in [20]:
                            subj = eval_subject_with_control(df, sp, baseline, cmax, recovery, depletion, rest_gain)
                            score = (
                                -0.08 * trial_eval["trial_mean_abs_err"]
                                -0.03 * trial_eval["trial_worst_abs_err"]
                                -0.001 * subj["rmse"]
                                -1.5 * subj["impossible_days"]
                                -0.002 * subj["mean_shortfall"]
                            )
                            rows.append({
                                "family": kind,
                                "sp_params": label,
                                "baseline": baseline,
                                "cmax": cmax,
                                "recovery": recovery,
                                "depletion": depletion,
                                "rest_gain": rest_gain,
                                "rmse": subj["rmse"],
                                "impossible_days": subj["impossible_days"],
                                "mean_shortfall": subj["mean_shortfall"],
                                "total_control": subj["total_control"],
                                "score": score,
                                **trial_eval,
                            })
    out = pd.DataFrame(rows).sort_values(["trial_mean_abs_err", "trial_worst_abs_err", "rmse"]).reset_index(drop=True)
    print("=" * 100)
    print("UNIFIED CONTROL-TRANSFER MODEL, TRIAL-FIRST")
    print("=" * 100)
    print("\nTop 30 unified fits:")
    print(out[[
        "family", "sp_params", "baseline", "cmax", "recovery", "depletion", "rest_gain",
        "trial_mean_abs_err", "trial_worst_abs_err", "SURMOUNT-4_mid", "STEP-1 extension_mid", "STEP-4_mid",
        "rmse", "impossible_days", "mean_shortfall", "total_control"
    ]].head(30).to_string(index=False, float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()
