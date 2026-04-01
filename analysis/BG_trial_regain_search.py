#!/usr/bin/env python3
"""BG. Trial-only search using published withdrawal anchors and uncertain regain partitions.

This script asks a narrower question than BF:
can one simple set-point law fit the major GLP-1/GIP withdrawal trials
without using the subject data at all?

The trials provide accurate body-weight anchors but not off-drug DXA.
So we:
  1. infer baseline and stop fat mass from drug-specific body-composition substudies
  2. simulate latent FM/SP dynamics
  3. map FM regain back to observed weight regain under uncertain regain partitions

This makes the uncertainty explicit instead of hiding it in a single hard-coded
fat fraction.

Published anchors used here:
  - SURMOUNT-4 withdrawal trial:
      baseline 107.3 kg, week-36 85.8 kg, placebo regain +14.0% by week 88
      Aronne et al. JAMA 2024, doi:10.1001/jama.2023.24945, PMID:38078870
  - STEP-1 extension:
      semaglutide extension subset 105.6 -> 87.5 -> 99.0 kg over weeks 0 -> 68 -> 120
      Wilding et al. Diabetes Obes Metab 2022, doi:10.1111/dom.14725, PMID:35441470
  - STEP-4:
      baseline 107.2 kg, week-20 96.1 kg, placebo regain +6.9% by week 68
      Rubino et al. JAMA 2021, doi:10.1001/jama.2021.3224, PMID:33755728

Body-composition priors used to infer FM:
  - SURMOUNT-1 DXA substudy:
      baseline FM 46.6 kg at body weight 102.5 kg, about 75% of lost weight as fat
      Heise et al. Obesity 2025, doi:10.1002/oby.24278, PMID:39996356
  - STEP-1 DXA substudy:
      baseline body-fat proportion 43.4%; FM -19.3%, BW -15.0%
      implies about 56% of lost weight as fat in that substudy
      Bergmann et al. Diabetes Obes Metab 2021, doi:10.1111/dom.14475, PMID:34036786
"""

from math import exp, log

import numpy as np
import pandas as pd

PRESSURE = 55.0
KCAL_PER_LB_FAT = 3500.0


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def sp_ema_series(fm, hl):
    alpha = 1 - exp(-log(2) / hl)
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        sp[i] = sp[i - 1] + alpha * (fm[i] - sp[i - 1])
    return sp


def sp_asym_series(fm, hl_up, hl_down):
    alpha_up = 1 - exp(-log(2) / hl_up)
    alpha_down = 1 - exp(-log(2) / hl_down)
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        alpha = alpha_up if fm[i] > sp[i - 1] else alpha_down
        sp[i] = sp[i - 1] + alpha * (fm[i] - sp[i - 1])
    return sp


def sp_latch_lookback_series(fm, tol, hold, rate):
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        if i >= hold and abs(fm[i] - fm[i - hold]) <= tol:
            sp[i] = sp[i - 1] + rate * (fm[i] - sp[i - 1])
        else:
            sp[i] = sp[i - 1]
    return sp


def sp_latch_hold_mean_series(fm, tol, hold, rate):
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


def sp_latch_hold_current_series(fm, tol, hold, rate):
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


def sp_vel_ema_series(fm, hl, vel_thresh, vel_scale):
    alpha = 1 - exp(-log(2) / hl)
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        j = max(0, i - 14)
        vel = abs((fm[i] - fm[j]) / max(1, i - j))
        gate = sigmoid((vel_thresh - vel) / vel_scale)
        sp[i] = sp[i - 1] + alpha * gate * (fm[i] - sp[i - 1])
    return sp


TRIALS = [
    {
        "name": "SURMOUNT-4",
        "drug": "tirzepatide",
        "baseline_weight": 107.3,
        "stop_weight": 85.8,
        "target_regain_pct": 14.0,
        "treatment_days": 36 * 7,
        "regain_days": 52 * 7,
        # Trial anchor:
        #   SURMOUNT-4 placebo-randomized cohort, mean body weight 107.3 kg at week 0
        #   and 85.8 kg at week 36; placebo regained +14.0% from week 36 to week 88.
        #   Source: doi:10.1001/jama.2023.24945
        # FM prior:
        #   SURMOUNT-1 DXA substudy reported 46.6 kg fat mass at 102.5 kg baseline body
        #   weight, so baseline BF proportion proxy = 46.6 / 102.5 = 45.5%.
        #   The same substudy reported about 75% of lost weight as fat mass, used here
        #   to infer stop-state FM in the absence of SURMOUNT-4 DXA.
        #   Source: doi:10.1002/oby.24278
        "baseline_bf": 46.6 / 102.5,
        "loss_fat_fraction": 0.75,
    },
    {
        "name": "STEP-1 extension",
        "drug": "semaglutide",
        "baseline_weight": 105.6,
        "stop_weight": 87.5,
        "target_regain_pct": (99.0 - 87.5) / 87.5 * 100.0,
        "treatment_days": 68 * 7,
        "regain_days": 52 * 7,
        # Trial anchor:
        #   STEP-1 extension semaglutide subset went 105.6 -> 87.5 -> 99.0 kg across
        #   weeks 0 -> 68 -> 120, so off-drug regain is measured over 52 weeks.
        #   Source: doi:10.1111/dom.14725
        # FM prior:
        #   STEP-1 DXA substudy reported baseline body-fat proportion 43.4%.
        #   It also reported FM -19.3% and BW -15.0%; using
        #     fat_fraction_of_loss ~= (0.434 * 0.193) / 0.150
        #   gives about 0.56 of lost weight as fat mass.
        #   Source: doi:10.1111/dom.14475
        "baseline_bf": 0.434,
        "loss_fat_fraction": 0.56,
    },
    {
        "name": "STEP-4",
        "drug": "semaglutide",
        "baseline_weight": 107.2,
        "stop_weight": 96.1,
        "target_regain_pct": 6.9,
        "treatment_days": 20 * 7,
        "regain_days": 48 * 7,
        # Trial anchor:
        #   STEP-4 randomized withdrawal cohort went 107.2 -> 96.1 kg during the
        #   20-week semaglutide run-in, then placebo regained +6.9% by week 68.
        #   Source: doi:10.1001/jama.2021.3224
        # FM prior:
        #   No STEP-4 DXA was located, so reuse the STEP-1 semaglutide composition priors.
        #   Source: doi:10.1111/dom.14475
        "baseline_bf": 0.434,
        "loss_fat_fraction": 0.56,
    },
]

# Off-drug DXA was not found for the major withdrawal cohorts. These scenarios treat
# regained body weight as 70-100% fat mass, spanning "fatter than the original loss"
# through "essentially all fat".
REGAIN_FAT_FRACTIONS = [1.00, 0.90, 0.80, 0.70]


def trial_anchor(trial):
    # Inferred stop-state FM = baseline FM minus the fat component of on-treatment loss.
    baseline_fm = trial["baseline_weight"] * trial["baseline_bf"]
    lost_weight = trial["baseline_weight"] - trial["stop_weight"]
    stop_fm = baseline_fm - lost_weight * trial["loss_fat_fraction"]
    return baseline_fm, stop_fm


def simulate_trial(trial, updater_factory, regain_fat_fraction):
    start_fm, stop_fm = trial_anchor(trial)
    fm_treat = np.linspace(start_fm, stop_fm, trial["treatment_days"])
    sp_treat = updater_factory("series", fm_treat)

    fm_hist = list(fm_treat)
    fm = stop_fm
    sp = sp_treat[-1]
    for _ in range(1, trial["regain_days"]):
        gap = sp - fm
        fm = fm + PRESSURE * gap / KCAL_PER_LB_FAT
        fm_hist.append(fm)
        sp = updater_factory("step", np.asarray(fm_hist), sp)

    # The latent model evolves fat mass. Observed weight regain depends on how much of
    # regained body weight is fat, which is uncertain because off-drug DXA was not found.
    regained_fm = fm - stop_fm
    regained_weight = regained_fm / regain_fat_fraction
    return regained_weight / trial["stop_weight"] * 100.0


def make_ema(hl):
    alpha = 1 - exp(-log(2) / hl)

    def updater(mode, fm_or_hist, sp_prev=None):
        if mode == "series":
            return sp_ema_series(fm_or_hist, hl)
        return sp_prev + alpha * (fm_or_hist[-1] - sp_prev)

    return updater


def make_asym(hl_up, hl_down):
    alpha_up = 1 - exp(-log(2) / hl_up)
    alpha_down = 1 - exp(-log(2) / hl_down)

    def updater(mode, fm_or_hist, sp_prev=None):
        if mode == "series":
            return sp_asym_series(fm_or_hist, hl_up, hl_down)
        alpha = alpha_up if fm_or_hist[-1] > sp_prev else alpha_down
        return sp_prev + alpha * (fm_or_hist[-1] - sp_prev)

    return updater


def make_latch(tol, hold, rate, kind):
    def updater(mode, fm_or_hist, sp_prev=None):
        if mode == "series":
            if kind == "lookback":
                return sp_latch_lookback_series(fm_or_hist, tol, hold, rate)
            if kind == "hold_mean":
                return sp_latch_hold_mean_series(fm_or_hist, tol, hold, rate)
            return sp_latch_hold_current_series(fm_or_hist, tol, hold, rate)

        hist = fm_or_hist
        i = len(hist) - 1
        if kind == "lookback":
            if i >= hold and abs(hist[i] - hist[i - hold]) <= tol:
                return sp_prev + rate * (hist[i] - sp_prev)
            return sp_prev
        if i < hold:
            return sp_prev
        window = np.asarray(hist[i - hold + 1:i + 1])
        if kind == "hold_mean":
            center = window.mean()
            if np.max(np.abs(window - center)) <= tol:
                return sp_prev + rate * (center - sp_prev)
            return sp_prev
        if np.max(np.abs(window - hist[i])) <= tol:
            return sp_prev + rate * (hist[i] - sp_prev)
        return sp_prev

    return updater


def make_vel_ema(hl, vel_thresh, vel_scale):
    alpha = 1 - exp(-log(2) / hl)

    def updater(mode, fm_or_hist, sp_prev=None):
        if mode == "series":
            return sp_vel_ema_series(fm_or_hist, hl, vel_thresh, vel_scale)
        i = len(fm_or_hist) - 1
        j = max(0, i - 14)
        vel = abs((fm_or_hist[i] - fm_or_hist[j]) / max(1, i - j))
        gate = sigmoid((vel_thresh - vel) / vel_scale)
        return sp_prev + alpha * gate * (fm_or_hist[-1] - sp_prev)

    return updater


def eval_model(model, params, updater_factory):
    vals = []
    by_trial = {}
    for regain_fat_fraction in REGAIN_FAT_FRACTIONS:
        errs = []
        for trial in TRIALS:
            pred = simulate_trial(trial, updater_factory, regain_fat_fraction)
            errs.append(abs(pred - trial["target_regain_pct"]))
            by_trial.setdefault(trial["name"], []).append(pred)
        vals.append(sum(errs))

    row = {
        "model": model,
        "params": params,
        "mean_abs_err": float(np.mean(vals)),
        "worst_abs_err": float(np.max(vals)),
    }
    for trial in TRIALS:
        preds = np.asarray(by_trial[trial["name"]])
        row[f"{trial['name']}_mid"] = preds[1]
        row[f"{trial['name']}_lo"] = preds.min()
        row[f"{trial['name']}_hi"] = preds.max()
        row[f"{trial['name']}_range"] = preds.max() - preds.min()
    return row


def main():
    rows = []

    for hl in [30, 45, 60, 80, 100, 120, 140, 160, 200]:
        rows.append(eval_model("EMA", f"hl={hl}", make_ema(hl)))

    for hl_up in [45, 60, 80, 100, 120]:
        for hl_down in [20, 30, 45, 60]:
            rows.append(eval_model("AsymEMA", f"up={hl_up},down={hl_down}", make_asym(hl_up, hl_down)))

    tol_grid = [1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]
    hold_grid = [7, 10, 14, 21, 28, 42, 56]
    rate_grid = [0.001, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02]
    for kind, label in [("lookback", "Lookback"), ("hold_mean", "HoldMean"), ("hold_current", "HoldCurrent")]:
        for tol in tol_grid:
            for hold in hold_grid:
                for rate in rate_grid:
                    rows.append(
                        eval_model(label, f"tol={tol},hold={hold},rate={rate}", make_latch(tol, hold, rate, kind))
                    )

    for hl in [30, 45, 60, 80]:
        for vel_thresh in [0.03, 0.05, 0.08, 0.10, 0.12]:
            rows.append(eval_model("VelEMA", f"hl={hl},v={vel_thresh}", make_vel_ema(hl, vel_thresh, 0.02)))

    out = pd.DataFrame(rows).sort_values(["mean_abs_err", "worst_abs_err"]).reset_index(drop=True)

    print("=" * 88)
    print("TRIAL-ONLY WITHDRAWAL SEARCH")
    print("=" * 88)
    print("Regain fat-fraction scenarios:", REGAIN_FAT_FRACTIONS)
    print("\nTrial anchors:")
    for trial in TRIALS:
        start_fm, stop_fm = trial_anchor(trial)
        print(
            f"  {trial['name']:<17} W {trial['baseline_weight']:.1f}->{trial['stop_weight']:.1f} kg,"
            f" inferred FM {start_fm:.1f}->{stop_fm:.1f} kg,"
            f" target regain {trial['target_regain_pct']:.2f}%"
        )

    cols = [
        "model",
        "params",
        "mean_abs_err",
        "worst_abs_err",
        "SURMOUNT-4_mid",
        "STEP-1 extension_mid",
        "STEP-4_mid",
    ]
    print("\nTop 30 by average trial error across regain partitions:")
    print(out[cols].head(30).to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print("\nMost robust candidates (lowest worst-case error):")
    print(out.sort_values(["worst_abs_err", "mean_abs_err"])[cols].head(30).to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print("\nBest by family:")
    for family in out["model"].unique():
        sub = out[out["model"] == family].head(1)
        print(sub[cols].to_string(index=False, float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()
