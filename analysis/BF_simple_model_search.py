#!/usr/bin/env python3
"""BF. Broad parameter search over corrected simple set-point models.

Searches simple latent-state families against three targets:
  1. Subject 2014+ pre-tirzepatide fit on 90d / 180d mean surplus
  2. SURMOUNT-4 post-discontinuation regain (~+14%)
  3. STEP-1 post-discontinuation regain (~+10%)

Families:
  - EMA
  - Asymmetric EMA
  - Lookback latch
  - Full-window hold latch around trailing mean
  - Full-window hold latch around current FM
  - Velocity-gated EMA
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PRESSURE = 55.0
FAT_FRACTION = 0.85


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def sp_ema_series(fm, hl):
    alpha = 1 - np.exp(-np.log(2) / hl)
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        sp[i] = sp[i - 1] + alpha * (fm[i] - sp[i - 1])
    return sp


def sp_asym_series(fm, hl_up, hl_down):
    alpha_up = 1 - np.exp(-np.log(2) / hl_up)
    alpha_down = 1 - np.exp(-np.log(2) / hl_down)
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
    alpha = 1 - np.exp(-np.log(2) / hl)
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        j = max(0, i - 14)
        vel = abs((fm[i] - fm[j]) / max(1, i - j))
        gate = sigmoid((vel_thresh - vel) / vel_scale)
        sp[i] = sp[i - 1] + alpha * gate * (fm[i] - sp[i - 1])
    return sp


def eval_subject(sp, fm, surplus90, surplus180, mask):
    dist = sp - fm
    out = {}
    for name, target in [("r90", surplus90), ("r180", surplus180)]:
        valid = mask & ~np.isnan(target) & ~np.isnan(dist)
        out[name] = np.corrcoef(dist[valid], target[valid])[0, 1] if valid.sum() > 200 else np.nan
    return out


def simulate_regain_ema(start_fm, stop_fm, treatment_days, regain_days, hl):
    fm_treat = np.linspace(start_fm, stop_fm, treatment_days)
    sp_treat = sp_ema_series(fm_treat, hl)
    fm = np.empty(regain_days)
    sp = np.empty(regain_days)
    fm[0] = stop_fm
    sp[0] = sp_treat[-1]
    total_weight_start = fm[0] / 0.40
    for d in range(1, regain_days):
        gap = sp[d - 1] - fm[d - 1]
        fm[d] = fm[d - 1] + PRESSURE * gap * FAT_FRACTION / 3500.0
        sp[d] = sp[d - 1] + (1 - np.exp(-np.log(2) / hl)) * (fm[d] - sp[d - 1])
    return (fm[-1] - fm[0]) / total_weight_start * 100.0


def simulate_regain_asym(start_fm, stop_fm, treatment_days, regain_days, hl_up, hl_down):
    fm_treat = np.linspace(start_fm, stop_fm, treatment_days)
    sp_treat = sp_asym_series(fm_treat, hl_up, hl_down)
    fm = np.empty(regain_days)
    sp = np.empty(regain_days)
    fm[0] = stop_fm
    sp[0] = sp_treat[-1]
    total_weight_start = fm[0] / 0.40
    alpha_up = 1 - np.exp(-np.log(2) / hl_up)
    alpha_down = 1 - np.exp(-np.log(2) / hl_down)
    for d in range(1, regain_days):
        gap = sp[d - 1] - fm[d - 1]
        fm[d] = fm[d - 1] + PRESSURE * gap * FAT_FRACTION / 3500.0
        alpha = alpha_up if fm[d] > sp[d - 1] else alpha_down
        sp[d] = sp[d - 1] + alpha * (fm[d] - sp[d - 1])
    return (fm[-1] - fm[0]) / total_weight_start * 100.0


def simulate_regain_latch(start_fm, stop_fm, treatment_days, regain_days, tol, hold, rate, mode):
    fm_treat = np.linspace(start_fm, stop_fm, treatment_days)
    if mode == "lookback":
        sp_treat = sp_latch_lookback_series(fm_treat, tol, hold, rate)
    elif mode == "hold_mean":
        sp_treat = sp_latch_hold_mean_series(fm_treat, tol, hold, rate)
    else:
        sp_treat = sp_latch_hold_current_series(fm_treat, tol, hold, rate)

    fm_hist = list(fm_treat)
    sp = sp_treat[-1]
    fm = stop_fm
    total_weight_start = fm / 0.40

    for _ in range(1, regain_days):
        gap = sp - fm
        fm = fm + PRESSURE * gap * FAT_FRACTION / 3500.0
        fm_hist.append(fm)
        i = len(fm_hist) - 1
        if mode == "lookback":
            if i >= hold and abs(fm_hist[i] - fm_hist[i - hold]) <= tol:
                sp = sp + rate * (fm_hist[i] - sp)
        elif mode == "hold_mean":
            if i >= hold:
                window = np.asarray(fm_hist[i - hold + 1:i + 1])
                center = window.mean()
                if np.max(np.abs(window - center)) <= tol:
                    sp = sp + rate * (center - sp)
        else:
            if i >= hold:
                window = np.asarray(fm_hist[i - hold + 1:i + 1])
                if np.max(np.abs(window - fm_hist[i])) <= tol:
                    sp = sp + rate * (fm_hist[i] - sp)

    return (fm - stop_fm) / total_weight_start * 100.0


def simulate_regain_vel_ema(start_fm, stop_fm, treatment_days, regain_days, hl, vel_thresh, vel_scale):
    fm_treat = np.linspace(start_fm, stop_fm, treatment_days)
    sp_treat = sp_vel_ema_series(fm_treat, hl, vel_thresh, vel_scale)
    fm_hist = list(fm_treat)
    sp = sp_treat[-1]
    fm = stop_fm
    total_weight_start = fm / 0.40
    alpha = 1 - np.exp(-np.log(2) / hl)
    for _ in range(1, regain_days):
        gap = sp - fm
        fm = fm + PRESSURE * gap * FAT_FRACTION / 3500.0
        fm_hist.append(fm)
        i = len(fm_hist) - 1
        j = max(0, i - 14)
        vel = abs((fm_hist[i] - fm_hist[j]) / max(1, i - j))
        gate = sigmoid((vel_thresh - vel) / vel_scale)
        sp = sp + alpha * gate * (fm - sp)
    return (fm - stop_fm) / total_weight_start * 100.0


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


def main():
    df, fm, mask = load_subject()
    rows = []

    for hl in [30, 45, 60, 80, 100, 120, 140, 160, 200]:
        sp = sp_ema_series(fm, hl)
        fit = eval_subject(sp, fm, df["surplus_90"].values, df["surplus_180"].values, mask)
        rows.append({
            "model": "EMA",
            "params": f"hl={hl}",
            "subject_r90": fit["r90"],
            "subject_r180": fit["r180"],
            "surmount4": simulate_regain_ema(90.0, 48.0, 252, 52 * 7, hl),
            "step1": simulate_regain_ema(92.0, 60.0, 476, 52 * 7, hl),
        })

    for hl_up in [45, 60, 80, 100, 120]:
        for hl_down in [20, 30, 45, 60]:
            sp = sp_asym_series(fm, hl_up, hl_down)
            fit = eval_subject(sp, fm, df["surplus_90"].values, df["surplus_180"].values, mask)
            rows.append({
                "model": "AsymEMA",
                "params": f"up={hl_up},down={hl_down}",
                "subject_r90": fit["r90"],
                "subject_r180": fit["r180"],
                "surmount4": simulate_regain_asym(90.0, 48.0, 252, 52 * 7, hl_up, hl_down),
                "step1": simulate_regain_asym(92.0, 60.0, 476, 52 * 7, hl_up, hl_down),
            })

    tol_grid = [1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]
    hold_grid = [7, 10, 14, 21, 28, 42, 56]
    rate_grid = [0.001, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02]
    for mode, fn in [
        ("Lookback", sp_latch_lookback_series),
        ("HoldMean", sp_latch_hold_mean_series),
        ("HoldCurrent", sp_latch_hold_current_series),
    ]:
        sim_mode = {"Lookback": "lookback", "HoldMean": "hold_mean", "HoldCurrent": "hold_current"}[mode]
        for tol in tol_grid:
            for hold in hold_grid:
                for rate in rate_grid:
                    sp = fn(fm, tol, hold, rate)
                    fit = eval_subject(sp, fm, df["surplus_90"].values, df["surplus_180"].values, mask)
                    rows.append({
                        "model": mode,
                        "params": f"tol={tol},hold={hold},rate={rate}",
                        "subject_r90": fit["r90"],
                        "subject_r180": fit["r180"],
                        "surmount4": simulate_regain_latch(90.0, 48.0, 252, 52 * 7, tol, hold, rate, sim_mode),
                        "step1": simulate_regain_latch(92.0, 60.0, 476, 52 * 7, tol, hold, rate, sim_mode),
                    })

    for hl in [30, 45, 60, 80]:
        for vel_thresh in [0.03, 0.05, 0.08, 0.10, 0.12]:
            sp = sp_vel_ema_series(fm, hl, vel_thresh, 0.02)
            fit = eval_subject(sp, fm, df["surplus_90"].values, df["surplus_180"].values, mask)
            rows.append({
                "model": "VelEMA",
                "params": f"hl={hl},v={vel_thresh}",
                "subject_r90": fit["r90"],
                "subject_r180": fit["r180"],
                "surmount4": simulate_regain_vel_ema(90.0, 48.0, 252, 52 * 7, hl, vel_thresh, 0.02),
                "step1": simulate_regain_vel_ema(92.0, 60.0, 476, 52 * 7, hl, vel_thresh, 0.02),
            })

    out = pd.DataFrame(rows)
    out["s4_err"] = (out["surmount4"] - 14.0).abs()
    out["step1_err"] = (out["step1"] - 10.0).abs()
    out["joint_err"] = out["s4_err"] + out["step1_err"]
    out["score"] = out["subject_r90"].abs() - 0.02 * out["s4_err"] - 0.01 * out["step1_err"]
    out = out.sort_values(["score", "subject_r90"], ascending=False).reset_index(drop=True)

    print("=" * 70)
    print("BROAD SIMPLE-MODEL SEARCH")
    print("=" * 70)
    print(f"\nTotal models tested: {len(out)}")

    print("\nTop 30 overall:")
    print(out[["model", "params", "subject_r90", "subject_r180", "surmount4", "step1", "s4_err", "step1_err"]]
          .head(30)
          .to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    good = out[out["subject_r90"].abs() > 0.75].copy()
    print("\nClosest to both trials among |subject_r90| > 0.75:")
    print(good.sort_values(["joint_err", "subject_r90"], ascending=[True, False])
          [["model", "params", "subject_r90", "subject_r180", "surmount4", "step1", "joint_err"]]
          .head(30)
          .to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print("\nBest by family:")
    for family in out["model"].unique():
        sub = out[out["model"] == family].head(1)
        print(sub[["model", "params", "subject_r90", "subject_r180", "surmount4", "step1"]]
              .to_string(index=False, float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()
