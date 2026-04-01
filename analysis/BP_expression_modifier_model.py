#!/usr/bin/env python3
"""BP. Shared lookback biology with a continuous expression modifier.

Same shared biology for subject and trials:
  pressure_t = 55 * (SP_t - FM_t)

Observation law:
  surplus_t = baseline + m_t * pressure_t + noise_t

Trials:
  m_t = 1

Subject:
  m_t evolves by one bounded overload/recovery law rather than hand-picked
  "uncensored" days.
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


def sp_lookback_series(fm, tol, hold, rate):
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        if i >= hold and abs(fm[i] - fm[i - hold]) <= tol:
            sp[i] = sp[i - 1] + rate * (fm[i] - sp[i - 1])
        else:
            sp[i] = sp[i - 1]
    return sp


def sp_lookback_step(hist, sp_prev, tol, hold, rate):
    i = len(hist) - 1
    if i >= hold and abs(hist[i] - hist[i - hold]) <= tol:
        return sp_prev + rate * (hist[i] - sp_prev)
    return sp_prev


def eval_trials(params):
    tol, hold, rate = params["tol"], params["hold"], params["rate"]
    totals = []
    mids = {}
    for frac in REGAIN_FAT_FRACTIONS:
        errs = []
        for trial in TRIALS:
            start_fm, stop_fm = trial_anchor(trial)
            fm_treat = np.linspace(start_fm, stop_fm, trial["treatment_days"])
            sp_treat = sp_lookback_series(fm_treat, tol, hold, rate)
            fm_hist = list(fm_treat)
            fm = stop_fm
            sp = sp_treat[-1]
            for _ in range(1, trial["regain_days"]):
                gap = sp - fm
                fm = fm + PRESSURE_PER_LB * gap / KCAL_PER_LB_FAT
                fm_hist.append(fm)
                sp = sp_lookback_step(fm_hist, sp, tol, hold, rate)
            pred = ((fm - stop_fm) / frac) / trial["stop_weight"] * 100.0
            errs.append(abs(pred - trial["target_regain_pct"]))
            mids.setdefault(trial["name"], []).append(pred)
        totals.append(sum(errs))
    out = {"trial_mean_abs_err": float(np.mean(totals))}
    for trial in TRIALS:
        out[f"{trial['name']}_mid"] = np.asarray(mids[trial["name"]])[1]
    return out


def evolve_modifier(pressure, surplus, baseline, m0, a, b, c):
    """
    m in [0, 1]: fraction of pressure expressed in behavior.

    overload grows when observed surplus sits far below baseline + full pressure.
    recovery grows when there is little or no suppressed pressure.
    """
    n = len(pressure)
    m = np.empty(n)
    pred = np.empty(n)
    overload = np.empty(n)
    recovery = np.empty(n)
    m[0] = m0
    pred[0] = baseline + m[0] * pressure[0]
    raw_gap0 = baseline + pressure[0] - surplus[0]
    overload[0] = max(0.0, raw_gap0 - c)
    recovery[0] = max(0.0, c - raw_gap0)
    for i in range(1, n):
        raw_gap = baseline + pressure[i - 1] - surplus[i - 1]
        overload[i] = max(0.0, raw_gap - c)
        recovery[i] = max(0.0, c - raw_gap)
        m[i] = m[i - 1] + a * overload[i] * (1 - m[i - 1]) - b * recovery[i] * m[i - 1]
        m[i] = float(np.clip(m[i], 0.0, 1.0))
        pred[i] = baseline + m[i] * pressure[i]
    return m, pred, overload, recovery


def main():
    df = load_subject()
    rows = []

    lookback_params = [
        {"tol": 2.5, "hold": 7, "rate": 0.001},
        {"tol": 3.0, "hold": 7, "rate": 0.001},
        {"tol": 3.0, "hold": 14, "rate": 0.001},
    ]

    for params in lookback_params:
        sp = sp_lookback_series(df["fm"].values, params["tol"], params["hold"], params["rate"])
        pressure = PRESSURE_PER_LB * (sp - df["fm"].values)
        trial_eval = eval_trials(params)
        for baseline in [-400, -300, -200, -100, 0]:
            for m0 in [0.0, 0.1, 0.2, 0.3]:
                for a in [0.0005, 0.001, 0.002, 0.005]:
                    for b in [0.0005, 0.001, 0.002, 0.005]:
                        for c in [0, 50, 100, 200]:
                            m, pred, overload, recovery = evolve_modifier(
                                pressure, df["surplus"].values, baseline, m0, a, b, c
                            )
                            rmse = float(np.sqrt(np.mean((df["surplus"].values - pred) ** 2)))
                            corr = float(np.corrcoef(pred, df["surplus"].values)[0, 1])
                            mean_m = float(m.mean())
                            frac_high = float((m > 0.8).mean())
                            frac_low = float((m < 0.2).mean())
                            score = (
                                -0.08 * trial_eval["trial_mean_abs_err"]
                                -0.001 * rmse
                                +0.50 * abs(corr)
                                -0.50 * abs(mean_m - 0.25)
                            )
                            rows.append({
                                "sp_params": f"tol={params['tol']},hold={params['hold']},rate={params['rate']}",
                                "baseline": baseline,
                                "m0": m0,
                                "a": a,
                                "b": b,
                                "c": c,
                                "rmse": rmse,
                                "corr": corr,
                                "mean_m": mean_m,
                                "frac_high": frac_high,
                                "frac_low": frac_low,
                                "score": score,
                                **trial_eval,
                            })

    out = pd.DataFrame(rows).sort_values(["score", "trial_mean_abs_err"], ascending=[False, True]).reset_index(drop=True)
    best = out.iloc[0]

    parts = dict(item.split("=") for item in best["sp_params"].split(","))
    tol = float(parts["tol"])
    hold = int(parts["hold"])
    rate = float(parts["rate"])
    sp = sp_lookback_series(df["fm"].values, tol, hold, rate)
    pressure = PRESSURE_PER_LB * (sp - df["fm"].values)
    m, pred, overload, recovery = evolve_modifier(
        pressure, df["surplus"].values, best["baseline"], best["m0"], best["a"], best["b"], best["c"]
    )

    artifact = df[["date", "fm", "surplus"]].copy()
    artifact["sp"] = sp
    artifact["pressure"] = pressure
    artifact["m_expr"] = m
    artifact["pred_surplus"] = pred
    artifact["overload"] = overload
    artifact["recovery"] = recovery

    artifact_path = ROOT / "analysis" / "BP_expression_modifier_daily.csv"
    search_path = ROOT / "analysis" / "BP_expression_modifier_search.csv"
    artifact.to_csv(artifact_path, index=False)
    out.to_csv(search_path, index=False)

    print("=" * 100)
    print("EXPRESSION MODIFIER MODEL")
    print("=" * 100)
    print("\nBest model:")
    print(best.to_string())

    print("\nTop 30 models:")
    print(out.head(30).to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print(f"\nArtifact: {search_path}")
    print(f"Artifact: {artifact_path}")


if __name__ == "__main__":
    main()
