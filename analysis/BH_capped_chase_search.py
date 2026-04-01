#!/usr/bin/env python3
"""BH. Search simple chase rules that preserve linear calorie pressure.

Keeps the same appetite law:
  pressure (cal/day) = k * (SP - FM)

But tries simpler alternatives to pure half-life or latch updates:
  - EMA
  - fixed-speed chase
  - capped EMA
  - asymmetric fixed-speed chase
  - asymmetric capped EMA

Evaluates both:
  1. subject fit on 2014+ pre-tirzepatide mean surplus
  2. trial-only withdrawal transfer using published weight anchors and uncertain
     regain fat fractions
"""

from math import exp, log
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


def make_ema_rule(hl):
    alpha = 1 - exp(-log(2) / hl)
    return lambda sp, fm: sp + alpha * (fm - sp)


def make_fixed_rule(speed):
    return lambda sp, fm: sp + np.clip(fm - sp, -speed, speed)


def make_capped_ema_rule(hl, cap):
    alpha = 1 - exp(-log(2) / hl)
    return lambda sp, fm: sp + np.clip(alpha * (fm - sp), -cap, cap)


def make_asym_fixed_rule(speed_up, speed_down):
    def rule(sp, fm):
        gap = fm - sp
        speed = speed_up if gap > 0 else speed_down
        return sp + np.clip(gap, -speed, speed)
    return rule


def make_asym_capped_ema_rule(hl_up, hl_down, cap_up, cap_down):
    alpha_up = 1 - exp(-log(2) / hl_up)
    alpha_down = 1 - exp(-log(2) / hl_down)

    def rule(sp, fm):
        gap = fm - sp
        if gap > 0:
            return sp + np.clip(alpha_up * gap, -cap_up, cap_up)
        return sp + np.clip(alpha_down * gap, -cap_down, cap_down)
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
    vals = []
    mids = {}
    for regain_fat_fraction in REGAIN_FAT_FRACTIONS:
        errs = []
        for trial in TRIALS:
            pred = simulate_trial(trial, rule, regain_fat_fraction)
            errs.append(abs(pred - trial["target_regain_pct"]))
            mids.setdefault(trial["name"], []).append(pred)
        vals.append(sum(errs))
    row = {
        "mean_abs_err": float(np.mean(vals)),
        "worst_abs_err": float(np.max(vals)),
    }
    for trial in TRIALS:
        arr = np.asarray(mids[trial["name"]])
        row[f"{trial['name']}_mid"] = arr[1]
        row[f"{trial['name']}_lo"] = arr.min()
        row[f"{trial['name']}_hi"] = arr.max()
    return row


def main():
    df, fm, mask = load_subject()
    rows = []

    for hl in [30, 45, 60, 80, 100, 120, 140, 160, 200]:
        rule = make_ema_rule(hl)
        sp = apply_rule_series(fm, rule)
        fit = eval_subject(sp, fm, df["surplus_90"].values, df["surplus_180"].values, mask)
        row = {"model": "EMA", "params": f"hl={hl}", "subject_r90": fit["r90"], "subject_r180": fit["r180"]}
        row.update(eval_trials(rule))
        rows.append(row)

    for speed in [0.01, 0.02, 0.03, 0.05, 0.08, 0.10, 0.15, 0.20]:
        rule = make_fixed_rule(speed)
        sp = apply_rule_series(fm, rule)
        fit = eval_subject(sp, fm, df["surplus_90"].values, df["surplus_180"].values, mask)
        row = {"model": "Fixed", "params": f"speed={speed}", "subject_r90": fit["r90"], "subject_r180": fit["r180"]}
        row.update(eval_trials(rule))
        rows.append(row)

    for hl in [30, 45, 60, 80, 100, 120, 160]:
        for cap in [0.01, 0.02, 0.03, 0.05, 0.08, 0.10, 0.15]:
            rule = make_capped_ema_rule(hl, cap)
            sp = apply_rule_series(fm, rule)
            fit = eval_subject(sp, fm, df["surplus_90"].values, df["surplus_180"].values, mask)
            row = {"model": "CappedEMA", "params": f"hl={hl},cap={cap}", "subject_r90": fit["r90"], "subject_r180": fit["r180"]}
            row.update(eval_trials(rule))
            rows.append(row)

    for speed_up in [0.03, 0.05, 0.08, 0.10]:
        for speed_down in [0.005, 0.01, 0.02, 0.03, 0.05]:
            rule = make_asym_fixed_rule(speed_up, speed_down)
            sp = apply_rule_series(fm, rule)
            fit = eval_subject(sp, fm, df["surplus_90"].values, df["surplus_180"].values, mask)
            row = {
                "model": "AsymFixed",
                "params": f"up={speed_up},down={speed_down}",
                "subject_r90": fit["r90"],
                "subject_r180": fit["r180"],
            }
            row.update(eval_trials(rule))
            rows.append(row)

    for hl_up in [45, 60, 80, 100]:
        for hl_down in [45, 60, 80, 120]:
            for cap_up in [0.03, 0.05, 0.08, 0.10]:
                for cap_down in [0.01, 0.02, 0.03, 0.05]:
                    rule = make_asym_capped_ema_rule(hl_up, hl_down, cap_up, cap_down)
                    sp = apply_rule_series(fm, rule)
                    fit = eval_subject(sp, fm, df["surplus_90"].values, df["surplus_180"].values, mask)
                    row = {
                        "model": "AsymCappedEMA",
                        "params": f"uphl={hl_up},downhl={hl_down},upcap={cap_up},downcap={cap_down}",
                        "subject_r90": fit["r90"],
                        "subject_r180": fit["r180"],
                    }
                    row.update(eval_trials(rule))
                    rows.append(row)

    out = pd.DataFrame(rows)
    out["joint_err"] = out["mean_abs_err"]
    out["score"] = out["subject_r90"].abs() - 0.03 * out["mean_abs_err"] - 0.01 * out["worst_abs_err"]
    out = out.sort_values(["score", "subject_r90"], ascending=False).reset_index(drop=True)

    cols = [
        "model", "params", "subject_r90", "subject_r180",
        "mean_abs_err", "worst_abs_err",
        "SURMOUNT-4_mid", "STEP-1 extension_mid", "STEP-4_mid",
    ]
    print("=" * 96)
    print("CAPPED/FIXED CHASE SEARCH")
    print("=" * 96)
    print(f"Total models tested: {len(out)}")
    print("\nTop 30 overall:")
    print(out[cols].head(30).to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    good = out[out["subject_r90"].abs() > 0.75].copy()
    print("\nClosest to trials among |subject_r90| > 0.75:")
    print(good.sort_values(["joint_err", "subject_r90"], ascending=[True, False])[cols]
          .head(30)
          .to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print("\nBest by family:")
    for family in out["model"].unique():
        sub = out[out["model"] == family].head(1)
        print(sub[cols].to_string(index=False, float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()
