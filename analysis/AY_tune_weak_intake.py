#!/usr/bin/env python3
"""AY. Tune weak-intake interpolation strength.

Searches a small grid over intake/TDEE prior weights and scores each run on:
  1. anchor fidelity (weight + FM/FFM scan fit)
  2. implied-TDEE plausibility (spread and tails)

Lower score is better.
"""

from pathlib import Path
import itertools
import subprocess
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
AX = ROOT / "analysis" / "AX_weak_intake_bodyfat_interpolation.py"
PY = ROOT / ".venv" / "bin" / "python"


def main():
    weight = pd.read_csv(ROOT / "analysis" / "P1_smoothed_weight.csv", parse_dates=["date"])[
        ["date", "smoothed_weight_lbs"]
    ]
    comp = pd.read_csv(ROOT / "composition" / "composition.csv", parse_dates=["date"])
    comp = comp[["date", "fat_mass_lbs", "lean_mass_lbs", "weight_lbs"]].dropna(
        subset=["fat_mass_lbs", "lean_mass_lbs"]
    )
    comp = comp.sort_values("date").drop_duplicates("date", keep="last")

    results = []
    out = ROOT / "analysis" / "_ax_tmp.csv"

    grid = itertools.product(
        [0.0, 0.0002, 0.0005, 0.0010, 0.0015],
        [0.0, 0.0002, 0.0005, 0.0008],
        [0.0, 0.0005, 0.0010, 0.0020],
        [0.0, 0.0020, 0.0050, 0.0100],
    )

    for we, wt, tf, ts in grid:
        cmd = [
            str(PY), str(AX),
            "--w-energy", str(we),
            "--w-tdee-prior", str(wt),
            "--l2-tdee-first", str(tf),
            "--l2-tdee-second", str(ts),
            "--out-file", str(out),
        ]
        cp = subprocess.run(cmd, capture_output=True, text=True)
        if cp.returncode != 0:
            results.append({
                "w_energy": we, "w_tdee_prior": wt, "l2_tdee_first": tf, "l2_tdee_second": ts,
                "score": np.inf, "status": "fail",
            })
            continue

        df = pd.read_csv(out)
        df["date"] = pd.to_datetime(df["date"])
        wsub = df.merge(weight, on="date", how="inner", suffixes=("", "_obs"))
        asub = df.merge(comp, on="date", how="inner")

        weight_rmse = float(np.sqrt(np.mean((wsub["weight_pred_lbs"] - wsub["smoothed_weight_lbs"]) ** 2)))
        fm_rmse = float(np.sqrt(np.mean((asub["fm_lbs_weak"] - asub["fat_mass_lbs"]) ** 2)))
        ffm_rmse = float(np.sqrt(np.mean((asub["ffm_lbs_weak"] - asub["lean_mass_lbs"]) ** 2)))

        implied = df["calories"].values[:-1] - 3500 * np.diff(df["fm_lbs_weak"].values)
        tdee_sd = float(np.nanstd(implied))
        tdee_p1 = float(np.nanpercentile(implied, 1))
        tdee_p99 = float(np.nanpercentile(implied, 99))

        # Score prioritizes anchor fidelity, then plausibility.
        score = (
            2.0 * fm_rmse +
            1.0 * ffm_rmse +
            1.5 * weight_rmse +
            max(0.0, (tdee_sd - 180.0) / 30.0) +
            max(0.0, (1850.0 - tdee_p1) / 100.0) +
            max(0.0, (tdee_p99 - 2450.0) / 100.0)
        )

        results.append({
            "w_energy": we,
            "w_tdee_prior": wt,
            "l2_tdee_first": tf,
            "l2_tdee_second": ts,
            "weight_rmse": weight_rmse,
            "fm_rmse": fm_rmse,
            "ffm_rmse": ffm_rmse,
            "tdee_sd": tdee_sd,
            "tdee_p1": tdee_p1,
            "tdee_p99": tdee_p99,
            "score": score,
            "status": "ok",
        })

    res = pd.DataFrame(results).sort_values(["score", "fm_rmse", "weight_rmse"])
    print(res.head(20).to_string(index=False, float_format=lambda x: f"{x:0.4f}"))


if __name__ == "__main__":
    main()
