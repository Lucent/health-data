#!/usr/bin/env python3
"""BK. Infer a latent depletable control stock from observed intake and set-point pressure.

Idea:
  observed surplus = baseline + lipostat_pressure - control_exertion + noise

where:
  lipostat_pressure = k * (SP - FM)

Instead of treating "willpower" as a fixed trait, model it as a stock:
  control_stock[t+1] = control_stock[t] + recovery - depletion

The stock is not directly observed. We infer the minimum control exertion needed to
explain why observed intake stayed below baseline + pressure, then ask whether a
simple depletable stock can sustain that pattern.

This is not yet a full forward behavioral model. It is a first consistency check:
does the subject history look like a finite control reservoir being spent down during
long restriction and partially recovering between bouts?
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PRESSURE_PER_LB = 55.0


def load_data():
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
    df["fm"] = fm
    df["surplus"] = df["calories"] - df["tdee"]
    return df


def sp_ema(fm, hl):
    alpha = 1 - np.exp(-np.log(2) / hl)
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        sp[i] = sp[i - 1] + alpha * (fm[i] - sp[i - 1])
    return sp


def infer_control_stock(required, cmax, recovery, depletion, rest_gain):
    """
    required: daily control expenditure required to explain observed intake
    cmax: stock capacity (cal/day-equivalent)
    recovery: fraction of remaining capacity recovered per day
    depletion: fraction of exerted control removed from stock
    rest_gain: small extra recovery on low-demand days
    """
    n = len(required)
    stock = np.empty(n)
    exerted = np.empty(n)
    shortfall = np.empty(n)
    stock[0] = cmax
    exerted[0] = min(stock[0], required[0])
    shortfall[0] = required[0] - exerted[0]

    for i in range(1, n):
        prev_req = required[i - 1]
        recover = recovery * (cmax - stock[i - 1])
        if prev_req < 100:
            recover += rest_gain
        stock_now = np.clip(stock[i - 1] + recover - depletion * exerted[i - 1], 0, cmax)
        stock[i] = stock_now
        exerted[i] = min(stock_now, required[i])
        shortfall[i] = required[i] - exerted[i]

    return stock, exerted, shortfall


def summarize_episode(df, start, end):
    sub = df[(df["date"] >= start) & (df["date"] <= end)]
    if sub.empty:
        return None
    return {
        "start": start,
        "end": end,
        "mean_surplus": sub["surplus"].mean(),
        "mean_pressure": sub["pressure"].mean(),
        "mean_required": sub["required_control"].mean(),
        "mean_stock": sub["control_stock"].mean(),
        "mean_shortfall": sub["control_shortfall"].mean(),
        "fm_start": sub["fm"].iloc[0],
        "fm_end": sub["fm"].iloc[-1],
    }


def main():
    df = load_data()
    mask = (df["date"] >= "2011-04-21") & (df["date"] < "2017-01-01") & (df["effective_level"] == 0)
    sub = df.loc[mask].copy().reset_index(drop=True)

    rows = []
    for hl in [45, 60, 80, 100]:
        sp = sp_ema(sub["fm"].values, hl)
        pressure = PRESSURE_PER_LB * (sp - sub["fm"].values)

        # Baseline is the intake surplus offset that minimizes mean required control.
        for baseline in [-400, -300, -200, -100, 0]:
            required = np.maximum(0.0, baseline + pressure - sub["surplus"].values)
            for cmax in [500, 900, 1200]:
                for recovery in [0.005, 0.01, 0.02]:
                    for depletion in [0.4, 0.7, 1.0]:
                        for rest_gain in [0, 10, 20]:
                            stock, exerted, shortfall = infer_control_stock(
                                required, cmax, recovery, depletion, rest_gain
                            )
                            impossible_days = (shortfall > 100).mean()
                            mean_shortfall = shortfall.mean()
                            # We want a coherent stock that explains needed control with few impossible days.
                            # Also prefer meaningful use of the stock rather than trivial zero-demand settings.
                            mean_required = required.mean()
                            score = mean_required - 3 * mean_shortfall - 800 * impossible_days
                            rows.append(
                                {
                                    "hl": hl,
                                    "baseline": baseline,
                                    "cmax": cmax,
                                    "recovery": recovery,
                                    "depletion": depletion,
                                    "rest_gain": rest_gain,
                                    "mean_required": mean_required,
                                    "mean_shortfall": mean_shortfall,
                                    "impossible_days": impossible_days,
                                    "score": score,
                                }
                            )

    out = pd.DataFrame(rows).sort_values(["score", "mean_shortfall"], ascending=[False, True]).reset_index(drop=True)
    best = out.iloc[0]

    sp = sp_ema(sub["fm"].values, int(best["hl"]))
    pressure = PRESSURE_PER_LB * (sp - sub["fm"].values)
    required = np.maximum(0.0, best["baseline"] + pressure - sub["surplus"].values)
    stock, exerted, shortfall = infer_control_stock(
        required,
        best["cmax"],
        best["recovery"],
        best["depletion"],
        best["rest_gain"],
    )

    sub["sp"] = sp
    sub["pressure"] = pressure
    sub["required_control"] = required
    sub["control_stock"] = stock
    sub["control_exerted"] = exerted
    sub["control_shortfall"] = shortfall
    sub["surplus_30"] = sub["surplus"].rolling(30, min_periods=30).mean()
    sub["required_30"] = sub["required_control"].rolling(30, min_periods=30).mean()
    sub["stock_30"] = sub["control_stock"].rolling(30, min_periods=30).mean()
    sub["shortfall_30"] = sub["control_shortfall"].rolling(30, min_periods=30).mean()

    print("=" * 88)
    print("CONTROL STOCK MODEL")
    print("=" * 88)
    print("\nBest parameters:")
    print(best.to_string())

    print("\nTop 20 parameter sets:")
    print(
        out[
            [
                "hl",
                "baseline",
                "cmax",
                "recovery",
                "depletion",
                "rest_gain",
                "mean_required",
                "mean_shortfall",
                "impossible_days",
                "score",
            ]
        ]
        .head(20)
        .to_string(index=False, float_format=lambda x: f"{x:.3f}")
    )

    print("\nSelected episodes:")
    for start, end in [
        ("2011-05-01", "2012-12-31"),
        ("2013-01-01", "2013-12-31"),
        ("2014-01-01", "2014-12-31"),
        ("2015-01-01", "2015-12-31"),
        ("2016-01-01", "2016-12-31"),
    ]:
        s = summarize_episode(sub, start, end)
        if s is None:
            continue
        print(
            f"  {start}..{end}  surplus={s['mean_surplus']:+.0f}  pressure={s['mean_pressure']:+.0f}"
            f"  required={s['mean_required']:.0f}  stock={s['mean_stock']:.0f}"
            f"  shortfall={s['mean_shortfall']:.0f}  FM {s['fm_start']:.1f}->{s['fm_end']:.1f}"
        )

    # Save a compact CSV for inspection.
    out_csv = ROOT / "analysis" / "BK_control_stock_daily.csv"
    sub[
        [
            "date",
            "fm",
            "sp",
            "surplus",
            "pressure",
            "required_control",
            "control_stock",
            "control_exerted",
            "control_shortfall",
            "surplus_30",
            "required_30",
            "stock_30",
            "shortfall_30",
        ]
    ].to_csv(out_csv, index=False)
    print(f"\nArtifact: {out_csv}")


if __name__ == "__main__":
    main()
