#!/usr/bin/env python3
"""BD. Simple set-point models with corrected hold logic and trial transfer checks.

This script keeps the simple-model family and fixes earlier implementation issues:
  - "Hold" can mean the full trailing window stayed within a ±tol band, not just
    that today matches the value `hold` days ago.
  - Subject fit and regain transfer are reported side by side for several nearby
    simple models rather than hard-coding one preferred narrative.
  - SURMOUNT-4 and STEP-1 regain are simulated from the same latent rule.

Models compared:
  - EMA
  - Asymmetric EMA
  - Lookback latch
  - True-hold latch around trailing mean
  - True-hold latch around current FM
  - Velocity-gated EMA
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

PRESSURE = 55.0
TOLERANCE = 3.0
HOLD_DAYS = 14
RATE_GRID = [0.0025, 0.005, 0.0075, 0.01, 0.015]
EMA_HL_GRID = [30, 45, 60, 80, 100, 120, 160]
ASYM_GRID = [(45, 45), (60, 30), (80, 30), (100, 30), (120, 45)]
VEL_THRESH_GRID = [0.04, 0.06, 0.08, 0.10]
VEL_SCALE = 0.02
FAT_FRACTION = 0.85


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def sp_ema(fm, hl):
    alpha = 1 - np.exp(-np.log(2) / hl)
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        sp[i] = sp[i - 1] + alpha * (fm[i] - sp[i - 1])
    return sp


def sp_asymmetric_ema(fm, hl_up, hl_down):
    alpha_up = 1 - np.exp(-np.log(2) / hl_up)
    alpha_down = 1 - np.exp(-np.log(2) / hl_down)
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        alpha = alpha_up if fm[i] > sp[i - 1] else alpha_down
        sp[i] = sp[i - 1] + alpha * (fm[i] - sp[i - 1])
    return sp


def sp_latch_lookback(fm, tol=TOLERANCE, hold=HOLD_DAYS, rate=0.005):
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        if i >= hold and abs(fm[i] - fm[i - hold]) <= tol:
            sp[i] = sp[i - 1] + rate * (fm[i] - sp[i - 1])
        else:
            sp[i] = sp[i - 1]
    return sp


def sp_latch_true_hold_mean(fm, tol=TOLERANCE, hold=HOLD_DAYS, rate=0.005):
    """Adapt only if the entire trailing window sits within ±tol of its mean."""
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        if i >= hold:
            window = fm[i - hold + 1:i + 1]
            center = window.mean()
            stable = np.max(np.abs(window - center)) <= tol
            if stable:
                sp[i] = sp[i - 1] + rate * (center - sp[i - 1])
                continue
        sp[i] = sp[i - 1]
    return sp


def sp_latch_true_hold_current(fm, tol=TOLERANCE, hold=HOLD_DAYS, rate=0.005):
    """Adapt only if the entire trailing window sits within ±tol of current FM."""
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        if i >= hold:
            window = fm[i - hold + 1:i + 1]
            stable = np.max(np.abs(window - fm[i])) <= tol
            if stable:
                sp[i] = sp[i - 1] + rate * (fm[i] - sp[i - 1])
                continue
        sp[i] = sp[i - 1]
    return sp


def sp_velocity_gated_ema(fm, hl=45, vel_thresh=0.08, vel_scale=VEL_SCALE):
    """Continuous simple model: EMA gain shrinks when FM velocity is high."""
    alpha_base = 1 - np.exp(-np.log(2) / hl)
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        j = max(0, i - HOLD_DAYS)
        vel = abs((fm[i] - fm[j]) / max(1, i - j))
        gate = sigmoid((vel_thresh - vel) / vel_scale)
        sp[i] = sp[i - 1] + alpha_base * gate * (fm[i] - sp[i - 1])
    return sp


def simulate_regain_from_path(fm_treat, sp_treat, update_fn, weeks=52):
    days = weeks * 7
    fm = np.empty(days)
    sp = np.empty(days)
    fm[0] = fm_treat[-1]
    sp[0] = sp_treat[-1]
    total_weight_start = fm[0] / 0.40

    full_fm = list(fm_treat)
    full_sp = list(sp_treat)
    for d in range(1, days):
        gap = sp[d - 1] - fm[d - 1]
        surplus = PRESSURE * gap
        fm[d] = fm[d - 1] + surplus * FAT_FRACTION / 3500.0
        full_fm.append(fm[d])
        next_sp = update_fn(np.asarray(full_fm), np.asarray(full_sp))
        sp[d] = next_sp[-1]
        full_sp.append(sp[d])

    regain_lbs = fm[-1] - fm[0]
    regain_pct = regain_lbs / total_weight_start * 100
    return regain_pct, fm, sp


def evaluate_subject(df, sp):
    rows = []
    dist = sp - df["fat_mass_lbs_filled"].values
    for label, mask in [
        ("Full", np.ones(len(df), dtype=bool)),
        ("2014+ pre-tirz", ((df["date"] >= "2014-01-01") & (df["effective_level"] == 0)).values),
    ]:
        for win_name, target in [("90d", df["surplus_90"].values), ("180d", df["surplus_180"].values)]:
            valid = mask & ~np.isnan(dist) & ~np.isnan(target)
            if valid.sum() < 200:
                continue
            rows.append((label, win_name, np.corrcoef(dist[valid], target[valid])[0, 1]))
    return rows


def treatment_path(sp_original, fm_start, loss_lbs, treatment_days):
    fm_treat = np.linspace(fm_start, fm_start - loss_lbs, treatment_days)
    sp_init = np.empty(len(fm_treat))
    sp_init[0] = sp_original
    return fm_treat, sp_init


def make_recursive_updater(kind, params):
    if kind == "ema":
        hl = params["hl"]

        def updater(full_fm, full_sp):
            return sp_ema(full_fm, hl)
        return updater

    if kind == "asym":
        hl_up, hl_down = params["hl_up"], params["hl_down"]

        def updater(full_fm, full_sp):
            return sp_asymmetric_ema(full_fm, hl_up, hl_down)
        return updater

    if kind == "lookback":
        rate = params["rate"]

        def updater(full_fm, full_sp):
            return sp_latch_lookback(full_fm, rate=rate)
        return updater

    if kind == "hold_mean":
        rate = params["rate"]

        def updater(full_fm, full_sp):
            return sp_latch_true_hold_mean(full_fm, rate=rate)
        return updater

    if kind == "hold_current":
        rate = params["rate"]

        def updater(full_fm, full_sp):
            return sp_latch_true_hold_current(full_fm, rate=rate)
        return updater

    if kind == "vel_ema":
        hl = params["hl"]
        vel_thresh = params["vel_thresh"]

        def updater(full_fm, full_sp):
            return sp_velocity_gated_ema(full_fm, hl=hl, vel_thresh=vel_thresh)
        return updater

    raise ValueError(kind)


def fit_initial_treatment_sp(fm_treat, sp0, updater):
    sp = np.empty(len(fm_treat))
    sp[0] = sp0
    full_fm = [fm_treat[0]]
    for i in range(1, len(fm_treat)):
        full_fm.append(fm_treat[i])
        sp_path = updater(np.asarray(full_fm), np.zeros(len(full_fm)))
        sp[i] = sp0 if len(sp_path) == 0 else sp_path[-1]
    return sp


def run_trial_transfer(kind, params, loss_lbs, treatment_days, regain_weeks=52, sp_original=90.0, fm_start=90.0):
    updater = make_recursive_updater(kind, params)
    fm_treat, _ = treatment_path(sp_original, fm_start, loss_lbs, treatment_days)
    sp_treat = updater(fm_treat, np.zeros(len(fm_treat)))
    return simulate_regain_from_path(fm_treat, sp_treat, updater, weeks=regain_weeks)[0]


def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    df = intake[["date", "calories"]].merge(
        kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="inner"
    )
    df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df = df.sort_values("date").reset_index(drop=True)
    fm = df["fat_mass_lbs"].values.copy()
    first_valid = np.where(~np.isnan(fm))[0][0]
    fm[:first_valid] = fm[first_valid]
    for i in range(first_valid + 1, len(fm)):
        if np.isnan(fm[i]):
            fm[i] = fm[i - 1]
    df["fat_mass_lbs_filled"] = fm
    df["surplus"] = df["calories"] - df["tdee"]
    df["surplus_90"] = df["surplus"].rolling(90, min_periods=90).mean()
    df["surplus_180"] = df["surplus"].rolling(180, min_periods=180).mean()

    candidates = []

    for hl in EMA_HL_GRID:
        sp = sp_ema(fm, hl)
        fits = dict(((label, win), r) for label, win, r in evaluate_subject(df, sp))
        candidates.append({
            "model": "EMA",
            "params": f"hl={hl}",
            "subject_r90": fits.get(("2014+ pre-tirz", "90d"), np.nan),
            "subject_r180": fits.get(("2014+ pre-tirz", "180d"), np.nan),
            "surmount4": run_trial_transfer("ema", {"hl": hl}, loss_lbs=42, treatment_days=252),
            "step1": run_trial_transfer("ema", {"hl": hl}, loss_lbs=32, treatment_days=476, sp_original=92.0, fm_start=92.0),
        })

    for hl_up, hl_down in ASYM_GRID:
        sp = sp_asymmetric_ema(fm, hl_up, hl_down)
        fits = dict(((label, win), r) for label, win, r in evaluate_subject(df, sp))
        candidates.append({
            "model": "AsymEMA",
            "params": f"up={hl_up},down={hl_down}",
            "subject_r90": fits.get(("2014+ pre-tirz", "90d"), np.nan),
            "subject_r180": fits.get(("2014+ pre-tirz", "180d"), np.nan),
            "surmount4": run_trial_transfer("asym", {"hl_up": hl_up, "hl_down": hl_down}, loss_lbs=42, treatment_days=252),
            "step1": run_trial_transfer("asym", {"hl_up": hl_up, "hl_down": hl_down}, loss_lbs=32, treatment_days=476, sp_original=92.0, fm_start=92.0),
        })

    for rate in RATE_GRID:
        for model_name, fn, kind in [
            ("LatchLookback", sp_latch_lookback, "lookback"),
            ("LatchHoldMean", sp_latch_true_hold_mean, "hold_mean"),
            ("LatchHoldCurrent", sp_latch_true_hold_current, "hold_current"),
        ]:
            sp = fn(fm, rate=rate)
            fits = dict(((label, win), r) for label, win, r in evaluate_subject(df, sp))
            candidates.append({
                "model": model_name,
                "params": f"rate={rate:.4f}",
                "subject_r90": fits.get(("2014+ pre-tirz", "90d"), np.nan),
                "subject_r180": fits.get(("2014+ pre-tirz", "180d"), np.nan),
                "surmount4": run_trial_transfer(kind, {"rate": rate}, loss_lbs=42, treatment_days=252),
                "step1": run_trial_transfer(kind, {"rate": rate}, loss_lbs=32, treatment_days=476, sp_original=92.0, fm_start=92.0),
            })

    for hl in [45, 60, 80]:
        for vel_thresh in VEL_THRESH_GRID:
            sp = sp_velocity_gated_ema(fm, hl=hl, vel_thresh=vel_thresh)
            fits = dict(((label, win), r) for label, win, r in evaluate_subject(df, sp))
            candidates.append({
                "model": "VelEMA",
                "params": f"hl={hl},v={vel_thresh:.2f}",
                "subject_r90": fits.get(("2014+ pre-tirz", "90d"), np.nan),
                "subject_r180": fits.get(("2014+ pre-tirz", "180d"), np.nan),
                "surmount4": run_trial_transfer("vel_ema", {"hl": hl, "vel_thresh": vel_thresh}, loss_lbs=42, treatment_days=252),
                "step1": run_trial_transfer("vel_ema", {"hl": hl, "vel_thresh": vel_thresh}, loss_lbs=32, treatment_days=476, sp_original=92.0, fm_start=92.0),
            })

    out = pd.DataFrame(candidates)
    out["s4_err"] = (out["surmount4"] - 14.0).abs()
    out["step1_err"] = (out["step1"] - 10.0).abs()
    out["score"] = out["subject_r90"].abs() - 0.02 * out["s4_err"] - 0.01 * out["step1_err"]
    out = out.sort_values(["score", "subject_r90"], ascending=False).reset_index(drop=True)

    print("=" * 70)
    print("SIMPLE MODEL COMPARISON WITH CORRECTED HOLD LOGIC")
    print("=" * 70)
    print("\nTop 20 overall:")
    print(
        out[["model", "params", "subject_r90", "subject_r180", "surmount4", "step1", "s4_err", "step1_err"]]
        .head(20)
        .to_string(index=False, float_format=lambda x: f"{x:.3f}")
    )

    print("\nBest by model family:")
    for family in out["model"].unique():
        sub = out[out["model"] == family].head(1)
        print(sub[["model", "params", "subject_r90", "subject_r180", "surmount4", "step1"]]
              .to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print("\nReference:")
    print("  SURMOUNT-4 published regain: +14.0%")
    print("  STEP-1 published regain: ~+10.0%")
    print("  Better simple models should keep subject fit high while pushing both trial regains upward.")


if __name__ == "__main__":
    main()
