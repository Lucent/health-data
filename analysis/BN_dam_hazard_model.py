#!/usr/bin/env python3
"""BN. Reverse-model restraint failure as a hazard under suppressed lipostat pressure.

Idea:
  We do not try to predict the exact day "willpower runs out" from a hidden stock.
  Instead, we treat observed restraint as a dam holding back biological pressure.

Shared biological drive:
  pressure_t = 55 * (SP_t - FM_t)

Observed restraint:
  suppression_t = max(0, baseline + pressure_t - observed_surplus_t)

Dam load:
  a rolling or EMA accumulation of prior suppression, optionally weighted by pressure

Target:
  future overage event in the next H days

This script searches simple SP rules, baselines, suppression-load definitions, and
overage targets to see whether failure is predicted better by:
  - current pressure alone
  - accumulated suppression alone
  - a combination of both

Updated extension:
  also tests a true backlog state where suppressed pressure both accumulates and
  decays over time, optionally with saturation. This is closer to the "dam under
  sustained load" idea than raw rolling sums.
"""

from math import exp, log
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PRESSURE_PER_LB = 55.0


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


def ema_accumulate(x, hl):
    alpha = 1 - exp(-log(2) / hl)
    out = np.empty(len(x))
    out[0] = x[0]
    for i in range(1, len(x)):
        out[i] = out[i - 1] + alpha * (x[i] - out[i - 1])
    return out


def backlog_state(inp, hl, gain=1.0, cap=None):
    """Decaying backlog: old load relaxes while new suppressed pressure adds to it."""
    alpha = 1 - exp(-log(2) / hl)
    out = np.empty(len(inp))
    out[0] = gain * inp[0]
    for i in range(1, len(inp)):
        prev = (1 - alpha) * out[i - 1]
        add = gain * inp[i]
        if cap is None:
            out[i] = prev + add
        else:
            # Saturating accumulation: harder to build more when backlog is already high.
            sat = max(0.0, 1.0 - out[i - 1] / cap)
            out[i] = prev + add * sat
            out[i] = min(out[i], cap)
    return out


def local_reset_state(inp, pressure, mode):
    """
    No-carry / instant-reset variants.

    mode = "instant"
      load is just today's suppressed input

    mode = "pressure_only"
      load is just today's positive pressure

    mode = "supp_if_pressure"
      today's suppression only when pressure is positive
    """
    if mode == "instant":
        return inp.copy()
    if mode == "pressure_only":
        return np.maximum(pressure, 0).copy()
    return inp * (np.maximum(pressure, 0) > 0)


def make_future_event(surplus, threshold, horizon):
    event = np.zeros(len(surplus), dtype=float)
    for i in range(len(surplus)):
        j = min(len(surplus), i + horizon + 1)
        event[i] = float(np.any(surplus[i + 1:j] >= threshold)) if i + 1 < j else 0.0
    return event


def auc_rank(y, score):
    y = np.asarray(y).astype(bool)
    score = np.asarray(score)
    pos = score[y]
    neg = score[~y]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    # Mann-Whitney U / AUC
    ranks = pd.Series(score).rank(method="average").values
    rank_sum = ranks[y].sum()
    auc = (rank_sum - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))
    return float(auc)


def logloss(y, score):
    p = 1.0 / (1.0 + np.exp(-score))
    p = np.clip(p, 1e-6, 1 - 1e-6)
    y = np.asarray(y)
    return float(-(y * np.log(p) + (1 - y) * np.log(1 - p)).mean())


def standardize(x):
    x = np.asarray(x, dtype=float)
    sd = x.std()
    if sd == 0 or np.isnan(sd):
        return np.zeros_like(x)
    return (x - x.mean()) / sd


def main():
    df = load_subject()
    rows = []

    sp_models = [
        ("EMA", "hl=45", sp_ema_series(df["fm"].values, 45)),
        ("EMA", "hl=50", sp_ema_series(df["fm"].values, 50)),
        ("Lookback", "tol=3,hold=7,rate=0.001", sp_lookback_series(df["fm"].values, 3.0, 7, 0.001)),
        ("Lookback", "tol=2.5,hold=7,rate=0.001", sp_lookback_series(df["fm"].values, 2.5, 7, 0.001)),
        ("HoldMean", "tol=2.5,hold=28,rate=0.0075", sp_hold_mean_series(df["fm"].values, 2.5, 28, 0.0075)),
        ("HoldMean", "tol=3,hold=42,rate=0.0075", sp_hold_mean_series(df["fm"].values, 3.0, 42, 0.0075)),
    ]

    for sp_family, sp_params, sp in sp_models:
        pressure = PRESSURE_PER_LB * (sp - df["fm"].values)
        for baseline in [-400, -300, -200, -100, 0]:
            suppression = np.maximum(0.0, baseline + pressure - df["surplus"].values)
            for load_kind in [
                "ema_supp", "ema_weighted", "rollsum", "rollsum_weighted",
                "backlog", "backlog_weighted", "backlog_sat", "backlog_weighted_sat",
                "instant", "pressure_only", "supp_if_pressure",
            ]:
                for load_scale in [14, 30, 45, 60, 90]:
                    if load_kind == "ema_supp":
                        load = ema_accumulate(suppression, load_scale)
                    elif load_kind == "ema_weighted":
                        load = ema_accumulate(suppression * np.maximum(pressure, 0), load_scale)
                    elif load_kind == "rollsum":
                        load = pd.Series(suppression).rolling(load_scale, min_periods=1).sum().values
                    elif load_kind == "rollsum_weighted":
                        load = pd.Series(suppression * np.maximum(pressure, 0)).rolling(load_scale, min_periods=1).sum().values
                    elif load_kind == "backlog":
                        load = backlog_state(suppression, load_scale, gain=1.0, cap=None)
                    elif load_kind == "backlog_weighted":
                        load = backlog_state(suppression * np.maximum(pressure, 0), load_scale, gain=1.0, cap=None)
                    elif load_kind == "backlog_sat":
                        load = backlog_state(suppression, load_scale, gain=1.0, cap=5000.0)
                    elif load_kind == "instant":
                        load = local_reset_state(suppression, pressure, "instant")
                    elif load_kind == "pressure_only":
                        load = local_reset_state(suppression, pressure, "pressure_only")
                    elif load_kind == "supp_if_pressure":
                        load = local_reset_state(suppression, pressure, "supp_if_pressure")
                    else:
                        load = backlog_state(suppression * np.maximum(pressure, 0), load_scale, gain=1.0, cap=50000.0)

                    pz = standardize(np.maximum(pressure, 0))
                    lz = standardize(load)
                    sz = standardize(suppression)

                    for threshold in [300, 500, 750, 1000]:
                        for horizon in [3, 7, 14]:
                            y = make_future_event(df["surplus"].values, threshold, horizon)
                            event_rate = y.mean()
                            if event_rate < 0.03 or event_rate > 0.97:
                                continue

                            candidate_scores = [
                                ("pressure", pz),
                                ("load", lz),
                                ("suppression", sz),
                                ("pressure+load", pz + lz),
                                ("pressure+2load", pz + 2 * lz),
                                ("2pressure+load", 2 * pz + lz),
                                ("pressure+load+supp", pz + lz + 0.5 * sz),
                            ]
                            for score_name, score in candidate_scores:
                                auc = auc_rank(y, score)
                                ll = logloss(y, score)
                                rows.append(
                                    {
                                        "sp_family": sp_family,
                                        "sp_params": sp_params,
                                        "baseline": baseline,
                                        "load_kind": load_kind,
                                        "load_scale": load_scale,
                                        "threshold": threshold,
                                        "horizon": horizon,
                                        "score_name": score_name,
                                        "event_rate": event_rate,
                                        "auc": auc,
                                        "logloss": ll,
                                    }
                                )

    out = pd.DataFrame(rows).sort_values(["auc", "logloss"], ascending=[False, True]).reset_index(drop=True)
    best = out.iloc[0]

    # Reconstruct best artifact series.
    best_sp = None
    for sp_family, sp_params, sp in sp_models:
        if sp_family == best["sp_family"] and sp_params == best["sp_params"]:
            best_sp = sp
            break
    pressure = PRESSURE_PER_LB * (best_sp - df["fm"].values)
    suppression = np.maximum(0.0, best["baseline"] + pressure - df["surplus"].values)
    if best["load_kind"] == "ema_supp":
        load = ema_accumulate(suppression, int(best["load_scale"]))
    elif best["load_kind"] == "ema_weighted":
        load = ema_accumulate(suppression * np.maximum(pressure, 0), int(best["load_scale"]))
    elif best["load_kind"] == "rollsum":
        load = pd.Series(suppression).rolling(int(best["load_scale"]), min_periods=1).sum().values
    elif best["load_kind"] == "rollsum_weighted":
        load = pd.Series(suppression * np.maximum(pressure, 0)).rolling(int(best["load_scale"]), min_periods=1).sum().values
    elif best["load_kind"] == "backlog":
        load = backlog_state(suppression, int(best["load_scale"]), gain=1.0, cap=None)
    elif best["load_kind"] == "backlog_weighted":
        load = backlog_state(suppression * np.maximum(pressure, 0), int(best["load_scale"]), gain=1.0, cap=None)
    elif best["load_kind"] == "backlog_sat":
        load = backlog_state(suppression, int(best["load_scale"]), gain=1.0, cap=5000.0)
    elif best["load_kind"] == "instant":
        load = local_reset_state(suppression, pressure, "instant")
    elif best["load_kind"] == "pressure_only":
        load = local_reset_state(suppression, pressure, "pressure_only")
    elif best["load_kind"] == "supp_if_pressure":
        load = local_reset_state(suppression, pressure, "supp_if_pressure")
    else:
        load = backlog_state(suppression * np.maximum(pressure, 0), int(best["load_scale"]), gain=1.0, cap=50000.0)
    pz = standardize(np.maximum(pressure, 0))
    lz = standardize(load)
    sz = standardize(suppression)
    score_map = {
        "pressure": pz,
        "load": lz,
        "suppression": sz,
        "pressure+load": pz + lz,
        "pressure+2load": pz + 2 * lz,
        "2pressure+load": 2 * pz + lz,
        "pressure+load+supp": pz + lz + 0.5 * sz,
    }
    score = score_map[best["score_name"]]
    event = make_future_event(df["surplus"].values, int(best["threshold"]), int(best["horizon"]))

    artifact = df[["date", "fm", "surplus"]].copy()
    artifact["sp"] = best_sp
    artifact["pressure"] = pressure
    artifact["suppression"] = suppression
    artifact["load"] = load
    artifact["score"] = score
    artifact["future_event"] = event
    artifact_path = ROOT / "analysis" / "BN_dam_hazard_daily.csv"
    artifact.to_csv(artifact_path, index=False)

    print("=" * 96)
    print("DAM / FAILURE HAZARD SEARCH")
    print("=" * 96)
    print("\nBest model:")
    print(best.to_string())

    search_path = ROOT / "analysis" / "BN_dam_hazard_search.csv"
    out.to_csv(search_path, index=False)

    print("\nTop 30 models:")
    print(out.head(30).to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    stricter = out[(out["threshold"] >= 750) & (out["event_rate"] >= 0.05) & (out["event_rate"] <= 0.40)].copy()
    if not stricter.empty:
        print("\nBest stricter-event models (threshold >= 750, event rate 5-40%):")
        print(stricter.head(20).to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print("\nTop by SP family:")
    for family in out["sp_family"].unique():
        sub = out[out["sp_family"] == family].head(1)
        print(sub.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print(f"\nArtifact: {search_path}")
    print(f"\nArtifact: {artifact_path}")


if __name__ == "__main__":
    main()
