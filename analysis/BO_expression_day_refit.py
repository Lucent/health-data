#!/usr/bin/env python3
"""BO. Refit shared set-point biology on likely expression days only.

Hypothesis:
  the subject's underlying biology may match trial-like regain speed, but many days
  are censored by willpower/restraint. If so, the right comparison is:

    trials: mostly expression days
    subject: mixed expression and suppression days

This script:
  1. tries a small set of shared SP laws
  2. computes daily pressure = 55 * (SP - FM)
  3. identifies likely "expression days" where observed surplus is not strongly below
     baseline + pressure
  4. refits the subject pressure relation on those days only
  5. compares the same SP law's zero-willpower trial regain to the subject expression fit
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


def trial_sim(kind, params, trial, frac):
    start_fm, stop_fm = trial_anchor(trial)
    fm_treat = np.linspace(start_fm, stop_fm, trial["treatment_days"])
    if kind == "EMA":
        sp_treat = sp_ema_series(fm_treat, params["hl"])
    elif kind == "Lookback":
        sp_treat = sp_lookback_series(fm_treat, params["tol"], params["hold"], params["rate"])
    else:
        sp_treat = sp_hold_mean_series(fm_treat, params["tol"], params["hold"], params["rate"])

    fm_hist = list(fm_treat)
    fm = stop_fm
    sp = sp_treat[-1]
    for _ in range(1, trial["regain_days"]):
        gap = sp - fm
        fm = fm + PRESSURE_PER_LB * gap / KCAL_PER_LB_FAT
        fm_hist.append(fm)
        if kind == "EMA":
            sp = sp_ema_series(np.array([sp, fm]), params["hl"])[-1]
        elif kind == "Lookback":
            i = len(fm_hist) - 1
            if i >= params["hold"] and abs(fm_hist[i] - fm_hist[i - params["hold"]]) <= params["tol"]:
                sp = sp + params["rate"] * (fm_hist[i] - sp)
        else:
            i = len(fm_hist) - 1
            if i >= params["hold"]:
                window = np.asarray(fm_hist[i - params["hold"] + 1:i + 1])
                center = window.mean()
                if np.max(np.abs(window - center)) <= params["tol"]:
                    sp = sp + params["rate"] * (center - sp)

    regained_weight = (fm - stop_fm) / frac
    return regained_weight / trial["stop_weight"] * 100.0


def eval_trials(kind, params):
    totals = []
    mids = {}
    for frac in REGAIN_FAT_FRACTIONS:
        errs = []
        for trial in TRIALS:
            pred = trial_sim(kind, params, trial, frac)
            errs.append(abs(pred - trial["target_regain_pct"]))
            mids.setdefault(trial["name"], []).append(pred)
        totals.append(sum(errs))
    out = {"trial_mean_abs_err": float(np.mean(totals))}
    for trial in TRIALS:
        out[f"{trial['name']}_mid"] = np.asarray(mids[trial["name"]])[1]
    return out


def fit_expression_days(df, sp, baseline, suppress_thresh, require_positive_pressure, min_run):
    pressure = PRESSURE_PER_LB * (sp - df["fm"].values)
    suppression = np.maximum(0.0, baseline + pressure - df["surplus"].values)
    expr = suppression <= suppress_thresh
    if require_positive_pressure:
        expr &= pressure > 0

    if min_run > 1:
        run = np.zeros(len(expr), dtype=int)
        cur = 0
        for i, val in enumerate(expr):
            cur = cur + 1 if val else 0
            run[i] = cur
        expr = run >= min_run

    n = int(expr.sum())
    if n < 50:
        return None

    x = pressure[expr]
    y = df["surplus"].values[expr]
    x_mean = x.mean()
    y_mean = y.mean()
    denom = ((x - x_mean) ** 2).sum()
    slope = np.nan if denom == 0 else ((x - x_mean) * (y - y_mean)).sum() / denom
    intercept = y_mean - slope * x_mean if not np.isnan(slope) else np.nan
    corr = np.corrcoef(x, y)[0, 1] if n > 2 else np.nan
    rmse = float(np.sqrt(np.mean((y - (intercept + slope * x)) ** 2)))
    return {
        "pressure": pressure,
        "suppression": suppression,
        "expr_mask": expr,
        "n_expr": n,
        "expr_frac": n / len(df),
        "slope": float(slope),
        "intercept": float(intercept),
        "corr": float(corr),
        "rmse_expr": rmse,
    }


def main():
    df = load_subject()
    candidates = [
        ("EMA", {"hl": 45}, "hl=45", sp_ema_series(df["fm"].values, 45)),
        ("EMA", {"hl": 50}, "hl=50", sp_ema_series(df["fm"].values, 50)),
        ("Lookback", {"tol": 2.5, "hold": 7, "rate": 0.001}, "tol=2.5,hold=7,rate=0.001", sp_lookback_series(df["fm"].values, 2.5, 7, 0.001)),
        ("Lookback", {"tol": 3.0, "hold": 7, "rate": 0.001}, "tol=3,hold=7,rate=0.001", sp_lookback_series(df["fm"].values, 3.0, 7, 0.001)),
        ("HoldMean", {"tol": 2.5, "hold": 28, "rate": 0.0075}, "tol=2.5,hold=28,rate=0.0075", sp_hold_mean_series(df["fm"].values, 2.5, 28, 0.0075)),
        ("HoldMean", {"tol": 3.0, "hold": 42, "rate": 0.0075}, "tol=3,hold=42,rate=0.0075", sp_hold_mean_series(df["fm"].values, 3.0, 42, 0.0075)),
    ]

    rows = []
    for kind, params, label, sp in candidates:
        trial_eval = eval_trials(kind, params)
        for baseline in [-400, -300, -200, -100, 0]:
            for suppress_thresh in [0, 50, 100, 200]:
                for require_positive_pressure in [False, True]:
                    for min_run in [1, 3, 7]:
                        fit = fit_expression_days(df, sp, baseline, suppress_thresh, require_positive_pressure, min_run)
                        if fit is None:
                            continue
                        # Prefer trial-good models whose expression-day slope is near 1 and correlation is high.
                        score = (
                            -0.10 * trial_eval["trial_mean_abs_err"]
                            -0.60 * abs(fit["slope"] - 1.0)
                            +0.40 * abs(fit["corr"])
                            -0.0005 * fit["rmse_expr"]
                            -0.20 * abs(fit["expr_frac"] - 0.25)
                        )
                        rows.append({
                            "family": kind,
                            "sp_params": label,
                            "baseline": baseline,
                            "suppress_thresh": suppress_thresh,
                            "require_positive_pressure": require_positive_pressure,
                            "min_run": min_run,
                            "n_expr": fit["n_expr"],
                            "expr_frac": fit["expr_frac"],
                            "slope": fit["slope"],
                            "intercept": fit["intercept"],
                            "corr": fit["corr"],
                            "rmse_expr": fit["rmse_expr"],
                            "score": score,
                            **trial_eval,
                        })

    out = pd.DataFrame(rows).sort_values(["score", "trial_mean_abs_err"], ascending=[False, True]).reset_index(drop=True)
    best = out.iloc[0]

    # Save best expression mask artifact.
    best_sp = None
    best_kind = None
    for kind, params, label, sp in candidates:
        if kind == best["family"] and label == best["sp_params"]:
            best_sp = sp
            best_kind = kind
            break
    fit = fit_expression_days(
        df,
        best_sp,
        best["baseline"],
        best["suppress_thresh"],
        bool(best["require_positive_pressure"]),
        int(best["min_run"]),
    )
    artifact = df[["date", "fm", "surplus"]].copy()
    artifact["sp"] = best_sp
    artifact["pressure"] = fit["pressure"]
    artifact["suppression"] = fit["suppression"]
    artifact["expression_day"] = fit["expr_mask"].astype(int)
    artifact_path = ROOT / "analysis" / "BO_expression_day_daily.csv"
    artifact.to_csv(artifact_path, index=False)

    search_path = ROOT / "analysis" / "BO_expression_day_search.csv"
    out.to_csv(search_path, index=False)

    print("=" * 100)
    print("EXPRESSION-DAY REFIT")
    print("=" * 100)
    print("\nBest model:")
    print(best.to_string())

    print("\nTop 30 models:")
    print(out.head(30).to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print("\nTop by family:")
    for family in out["family"].unique():
        sub = out[out["family"] == family].head(1)
        print(sub.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print(f"\nArtifact: {search_path}")
    print(f"Artifact: {artifact_path}")


if __name__ == "__main__":
    main()
