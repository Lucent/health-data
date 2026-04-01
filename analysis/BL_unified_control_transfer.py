#!/usr/bin/env python3
"""BL. Unified set-point + control-transfer model.

Goal:
  use one underlying lipostat formula for both:
  1. published GLP-1 withdrawal trials (assume zero willpower/control)
  2. this subject (allow an exhaustible latent control stock)

Shared biology:
  pressure_t = k * (SP_t - FM_t)
  SP_t updates from FM_t by one chosen rule

Subject-only modifier:
  observed_surplus_t = baseline + pressure_t - control_exertion_t + noise_t
  control comes from a finite stock that depletes and recovers

Trials:
  same pressure and same SP rule, but control_exertion_t = 0

This does not claim a full causal behavioral model. It asks a narrower question:
how much of the subject-vs-trial gap can be closed by adding an exhaustible control
stock on top of the same lipostat law?
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
    {
        "name": "SURMOUNT-4",
        "baseline_weight": 107.3,
        "stop_weight": 85.8,
        "target_regain_pct": 14.0,
        "treatment_days": 36 * 7,
        "regain_days": 52 * 7,
        "baseline_bf": 46.6 / 102.5,
        "loss_fat_fraction": 0.75,
    },
    {
        "name": "STEP-1 extension",
        "baseline_weight": 105.6,
        "stop_weight": 87.5,
        "target_regain_pct": (99.0 - 87.5) / 87.5 * 100.0,
        "treatment_days": 68 * 7,
        "regain_days": 52 * 7,
        "baseline_bf": 0.434,
        "loss_fat_fraction": 0.56,
    },
    {
        "name": "STEP-4",
        "baseline_weight": 107.2,
        "stop_weight": 96.1,
        "target_regain_pct": 6.9,
        "treatment_days": 20 * 7,
        "regain_days": 48 * 7,
        "baseline_bf": 0.434,
        "loss_fat_fraction": 0.56,
    },
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


def sp_ema_step(sp_prev, fm, hl):
    alpha = 1 - exp(-log(2) / hl)
    return sp_prev + alpha * (fm - sp_prev)


def sp_hold_mean_step(hist, sp_prev, tol, hold, rate):
    i = len(hist) - 1
    if i >= hold:
        window = np.asarray(hist[i - hold + 1:i + 1])
        center = window.mean()
        if np.max(np.abs(window - center)) <= tol:
            return sp_prev + rate * (center - sp_prev)
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
    rmse = float(np.sqrt(np.mean(residual ** 2)))
    impossible_days = float((shortfall > 100).mean())
    total_control = float(exerted.sum())
    return {
        "pressure": pressure,
        "required": required,
        "stock": stock,
        "exerted": exerted,
        "shortfall": shortfall,
        "rmse": rmse,
        "impossible_days": impossible_days,
        "total_control": total_control,
        "mean_required": float(required.mean()),
        "mean_shortfall": float(shortfall.mean()),
    }


def simulate_trial_ema(trial, hl, regain_fat_fraction):
    start_fm, stop_fm = trial_anchor(trial)
    fm_treat = np.linspace(start_fm, stop_fm, trial["treatment_days"])
    sp_treat = sp_ema_series(fm_treat, hl)
    fm = stop_fm
    sp = sp_treat[-1]
    for _ in range(1, trial["regain_days"]):
        gap = sp - fm
        fm = fm + PRESSURE_PER_LB * gap / KCAL_PER_LB_FAT
        sp = sp_ema_step(sp, fm, hl)
    regained_fm = fm - stop_fm
    regained_weight = regained_fm / regain_fat_fraction
    return regained_weight / trial["stop_weight"] * 100.0


def simulate_trial_hold(trial, tol, hold, rate, regain_fat_fraction):
    start_fm, stop_fm = trial_anchor(trial)
    fm_treat = np.linspace(start_fm, stop_fm, trial["treatment_days"])
    sp_treat = sp_hold_mean_series(fm_treat, tol, hold, rate)
    fm_hist = list(fm_treat)
    fm = stop_fm
    sp = sp_treat[-1]
    for _ in range(1, trial["regain_days"]):
        gap = sp - fm
        fm = fm + PRESSURE_PER_LB * gap / KCAL_PER_LB_FAT
        fm_hist.append(fm)
        sp = sp_hold_mean_step(fm_hist, sp, tol, hold, rate)
    regained_fm = fm - stop_fm
    regained_weight = regained_fm / regain_fat_fraction
    return regained_weight / trial["stop_weight"] * 100.0


def eval_trials(sim_fn):
    total_errs = []
    mids = {}
    for frac in REGAIN_FAT_FRACTIONS:
        errs = []
        for trial in TRIALS:
            pred = sim_fn(trial, frac)
            errs.append(abs(pred - trial["target_regain_pct"]))
            mids.setdefault(trial["name"], []).append(pred)
        total_errs.append(sum(errs))
    out = {
        "trial_mean_abs_err": float(np.mean(total_errs)),
        "trial_worst_abs_err": float(np.max(total_errs)),
    }
    for trial in TRIALS:
        arr = np.asarray(mids[trial["name"]])
        out[f"{trial['name']}_mid"] = arr[1]
    return out


def episode_summary(df, start, end):
    sub = df[(df["date"] >= start) & (df["date"] <= end)]
    return {
        "period": f"{start}..{end}",
        "surplus": float(sub["surplus"].mean()),
        "pressure": float(sub["pressure"].mean()),
        "required": float(sub["required_control"].mean()),
        "stock": float(sub["control_stock"].mean()),
        "shortfall": float(sub["control_shortfall"].mean()),
        "cum_control": float(sub["control_exerted"].sum()),
        "fm_start": float(sub["fm"].iloc[0]),
        "fm_end": float(sub["fm"].iloc[-1]),
    }


def main():
    df = load_subject()
    rows = []

    model_specs = []
    for hl in [45, 50]:
        model_specs.append(("EMA", f"hl={hl}", sp_ema_series(df["fm"].values, hl), lambda trial, frac, hl=hl: simulate_trial_ema(trial, hl, frac)))
    for tol in [2.5, 3.0]:
        for hold in [28, 42]:
            for rate in [0.005, 0.0075]:
                model_specs.append((
                    "HoldMean",
                    f"tol={tol},hold={hold},rate={rate}",
                    sp_hold_mean_series(df["fm"].values, tol, hold, rate),
                    lambda trial, frac, tol=tol, hold=hold, rate=rate: simulate_trial_hold(trial, tol, hold, rate, frac),
                ))

    for family, params, sp, trial_sim in model_specs:
        trial_eval = eval_trials(trial_sim)
        for baseline in [-400, -300, -200]:
            for cmax in [900, 1200]:
                for recovery in [0.01, 0.02]:
                    for depletion in [0.4, 0.7]:
                        for rest_gain in [10, 20]:
                            subj = eval_subject_with_control(df, sp, baseline, cmax, recovery, depletion, rest_gain)
                            # Combined score: better subject explanatory power with fewer impossible days,
                            # but still reward models that transfer to trials without willpower.
                            score = (
                                -0.002 * subj["rmse"]
                                - 2.5 * subj["impossible_days"]
                                - 0.004 * subj["mean_shortfall"]
                                - 0.03 * trial_eval["trial_mean_abs_err"]
                                - 0.01 * trial_eval["trial_worst_abs_err"]
                            )
                            rows.append({
                                "family": family,
                                "sp_params": params,
                                "baseline": baseline,
                                "cmax": cmax,
                                "recovery": recovery,
                                "depletion": depletion,
                                "rest_gain": rest_gain,
                                "rmse": subj["rmse"],
                                "impossible_days": subj["impossible_days"],
                                "mean_required": subj["mean_required"],
                                "mean_shortfall": subj["mean_shortfall"],
                                "total_control": subj["total_control"],
                                "score": score,
                                **trial_eval,
                            })

    out = pd.DataFrame(rows).sort_values(["score", "rmse"], ascending=[False, True]).reset_index(drop=True)
    best = out.iloc[0]

    # Reconstruct best model series and save artifact.
    if best["family"] == "EMA":
        hl = int(best["sp_params"].split("=")[1])
        sp = sp_ema_series(df["fm"].values, hl)
    else:
        parts = dict(item.split("=") for item in best["sp_params"].split(","))
        tol = float(parts["tol"])
        hold = int(parts["hold"])
        rate = float(parts["rate"])
        sp = sp_hold_mean_series(df["fm"].values, tol, hold, rate)

    subj = eval_subject_with_control(
        df, sp,
        best["baseline"], best["cmax"], best["recovery"], best["depletion"], best["rest_gain"]
    )

    df["sp"] = sp
    df["pressure"] = subj["pressure"]
    df["required_control"] = subj["required"]
    df["control_stock"] = subj["stock"]
    df["control_exerted"] = subj["exerted"]
    df["control_shortfall"] = subj["shortfall"]
    df["surplus_pred_no_control"] = best["baseline"] + df["pressure"]
    df["surplus_pred_with_control"] = best["baseline"] + df["pressure"] - df["control_exerted"]

    artifact_path = ROOT / "analysis" / "BL_unified_control_daily.csv"
    df[
        [
            "date", "fm", "surplus", "sp", "pressure",
            "required_control", "control_stock", "control_exerted",
            "control_shortfall", "surplus_pred_no_control", "surplus_pred_with_control",
        ]
    ].to_csv(artifact_path, index=False)

    print("=" * 100)
    print("UNIFIED CONTROL-TRANSFER MODEL")
    print("=" * 100)
    print("\nBest model:")
    print(best.to_string())

    print("\nTop 20 unified fits:")
    print(
        out[
            [
                "family", "sp_params", "baseline", "cmax", "recovery", "depletion", "rest_gain",
                "rmse", "impossible_days", "mean_shortfall", "total_control",
                "trial_mean_abs_err", "trial_worst_abs_err",
                "SURMOUNT-4_mid", "STEP-1 extension_mid", "STEP-4_mid",
            ]
        ].head(20).to_string(index=False, float_format=lambda x: f"{x:.3f}")
    )

    print("\nSelected subject episodes:")
    for start, end in [
        ("2011-05-01", "2012-12-31"),
        ("2013-01-01", "2013-12-31"),
        ("2014-01-01", "2014-12-31"),
        ("2015-01-01", "2015-12-31"),
        ("2016-01-01", "2016-12-31"),
    ]:
        s = episode_summary(df, start, end)
        print(
            f"  {s['period']}  surplus={s['surplus']:+.0f}  pressure={s['pressure']:+.0f}"
            f"  required={s['required']:.0f}  stock={s['stock']:.0f}"
            f"  shortfall={s['shortfall']:.0f}  cum_control={s['cum_control']:.0f}"
            f"  FM {s['fm_start']:.1f}->{s['fm_end']:.1f}"
        )

    print(f"\nArtifact: {artifact_path}")


if __name__ == "__main__":
    main()
