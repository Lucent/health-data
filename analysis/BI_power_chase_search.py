#!/usr/bin/env python3
"""BI. Search sublinear power-law chase rules for the defended fat-mass level.

Keeps the appetite pressure law fixed:
  pressure (cal/day) = 55 * (SP - FM)

Changes only the update rule:
  SP_{t+1} = SP_t + a * sign(FM - SP) * |FM - SP|^p

This interpolates between:
  - near-fixed-speed chase when p is small
  - EMA-like proportional chase when p = 1

We test:
  - symmetric power chase
  - asymmetric power chase

against both subject fit and trial withdrawal transfer.
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PRESSURE = 55.0
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
    df["surplus_90"] = df["surplus"].rolling(90, min_periods=90).mean()
    df["surplus_180"] = df["surplus"].rolling(180, min_periods=180).mean()
    mask = ((df["date"] >= "2014-01-01") & (df["effective_level"] == 0)).values
    return df, fm, mask


def eval_subject(sp, fm, surplus90, surplus180, mask):
    dist = sp - fm
    out = {}
    for name, target in [("r90", surplus90), ("r180", surplus180)]:
        valid = mask & ~np.isnan(target) & ~np.isnan(dist)
        out[name] = np.corrcoef(dist[valid], target[valid])[0, 1] if valid.sum() > 200 else np.nan
    return out


def trial_anchor(trial):
    baseline_fm = trial["baseline_weight"] * trial["baseline_bf"]
    lost_weight = trial["baseline_weight"] - trial["stop_weight"]
    stop_fm = baseline_fm - lost_weight * trial["loss_fat_fraction"]
    return baseline_fm, stop_fm


def apply_rule_series(fm, rule):
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        sp[i] = rule(sp[i - 1], fm[i])
    return sp


def make_power_rule(scale, power, cap=None):
    def rule(sp, fm):
        gap = fm - sp
        step = scale * np.sign(gap) * (abs(gap) ** power)
        if cap is not None:
            step = np.clip(step, -cap, cap)
        return sp + step
    return rule


def make_asym_power_rule(scale_up, power_up, scale_down, power_down, cap_up=None, cap_down=None):
    def rule(sp, fm):
        gap = fm - sp
        if gap >= 0:
            step = scale_up * (gap ** power_up)
            if cap_up is not None:
                step = min(step, cap_up)
        else:
            step = -scale_down * ((-gap) ** power_down)
            if cap_down is not None:
                step = max(step, -cap_down)
        return sp + step
    return rule


def simulate_trial(trial, rule, regain_fat_fraction):
    start_fm, stop_fm = trial_anchor(trial)
    fm_treat = np.linspace(start_fm, stop_fm, trial["treatment_days"])
    sp_treat = apply_rule_series(fm_treat, rule)

    fm = stop_fm
    sp = sp_treat[-1]
    for _ in range(1, trial["regain_days"]):
        gap = sp - fm
        fm = fm + PRESSURE * gap / KCAL_PER_LB_FAT
        sp = rule(sp, fm)

    regained_fm = fm - stop_fm
    regained_weight = regained_fm / regain_fat_fraction
    return regained_weight / trial["stop_weight"] * 100.0


def eval_trials(rule):
    errs_total = []
    preds = {}
    for regain_fat_fraction in REGAIN_FAT_FRACTIONS:
        errs = []
        for trial in TRIALS:
            pred = simulate_trial(trial, rule, regain_fat_fraction)
            errs.append(abs(pred - trial["target_regain_pct"]))
            preds.setdefault(trial["name"], []).append(pred)
        errs_total.append(sum(errs))
    row = {
        "mean_abs_err": float(np.mean(errs_total)),
        "worst_abs_err": float(np.max(errs_total)),
    }
    for trial in TRIALS:
        arr = np.asarray(preds[trial["name"]])
        row[f"{trial['name']}_mid"] = arr[1]
        row[f"{trial['name']}_lo"] = arr.min()
        row[f"{trial['name']}_hi"] = arr.max()
    return row


def main():
    df, fm, mask = load_subject()
    rows = []

    for power in [0.25, 0.33, 0.50, 0.67, 0.80, 1.00]:
        for scale in [0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.10]:
            for cap in [None, 0.03, 0.05, 0.08, 0.10]:
                rule = make_power_rule(scale, power, cap)
                sp = apply_rule_series(fm, rule)
                fit = eval_subject(sp, fm, df["surplus_90"].values, df["surplus_180"].values, mask)
                label = f"scale={scale},p={power}" + ("" if cap is None else f",cap={cap}")
                row = {
                    "model": "Power",
                    "params": label,
                    "subject_r90": fit["r90"],
                    "subject_r180": fit["r180"],
                }
                row.update(eval_trials(rule))
                rows.append(row)

    for power_up in [0.33, 0.50, 0.67, 0.80]:
        for power_down in [0.33, 0.50, 0.67, 0.80, 1.00]:
            for scale_up in [0.01, 0.02, 0.03, 0.05]:
                for scale_down in [0.003, 0.005, 0.01, 0.02]:
                    rule = make_asym_power_rule(scale_up, power_up, scale_down, power_down)
                    sp = apply_rule_series(fm, rule)
                    fit = eval_subject(sp, fm, df["surplus_90"].values, df["surplus_180"].values, mask)
                    row = {
                        "model": "AsymPower",
                        "params": f"up={scale_up}^{power_up},down={scale_down}^{power_down}",
                        "subject_r90": fit["r90"],
                        "subject_r180": fit["r180"],
                    }
                    row.update(eval_trials(rule))
                    rows.append(row)

    out = pd.DataFrame(rows)
    out["score"] = out["subject_r90"].abs() - 0.03 * out["mean_abs_err"] - 0.01 * out["worst_abs_err"]
    out = out.sort_values(["score", "subject_r90"], ascending=False).reset_index(drop=True)

    cols = [
        "model", "params", "subject_r90", "subject_r180",
        "mean_abs_err", "worst_abs_err",
        "SURMOUNT-4_mid", "STEP-1 extension_mid", "STEP-4_mid",
    ]
    print("=" * 100)
    print("POWER-CHASE SEARCH")
    print("=" * 100)
    print(f"Total models tested: {len(out)}")
    print("\nTop 30 overall:")
    print(out[cols].head(30).to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    good = out[out["subject_r90"].abs() > 0.75].copy()
    print("\nClosest to trials among |subject_r90| > 0.75:")
    print(good.sort_values(["mean_abs_err", "subject_r90"], ascending=[True, False])[cols]
          .head(30)
          .to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print("\nBest by family:")
    for family in out["model"].unique():
        sub = out[out["model"] == family].head(1)
        print(sub[cols].to_string(index=False, float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()
