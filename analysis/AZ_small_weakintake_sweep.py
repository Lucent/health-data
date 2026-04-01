#!/usr/bin/env python3
"""AZ. Small curated sweep for weak-intake interpolation.

Faster companion to AY for interactive tuning.
"""

from pathlib import Path
import subprocess
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
AX = ROOT / "analysis" / "AX_weak_intake_bodyfat_interpolation.py"
PY = ROOT / ".venv" / "bin" / "python"

CONFIGS = [
    (0.0,    0.0,    0.0,    0.0),
    (0.0002, 0.0002, 0.0005, 0.0020),
    (0.0003, 0.0002, 0.0005, 0.0010),
    (0.0003, 0.0002, 0.0010, 0.0020),
    (0.0005, 0.0002, 0.0005, 0.0020),
    (0.0005, 0.0005, 0.0005, 0.0020),
    (0.0007, 0.0003, 0.0010, 0.0020),
    (0.0007, 0.0005, 0.0010, 0.0050),
    (0.0010, 0.0005, 0.0010, 0.0050),
    (0.0015, 0.0008, 0.0020, 0.0100),
]


def main():
    weight = pd.read_csv(ROOT / "analysis" / "P1_smoothed_weight.csv", parse_dates=["date"])[
        ["date", "smoothed_weight_lbs"]
    ]
    comp = pd.read_csv(ROOT / "composition" / "composition.csv", parse_dates=["date"])
    comp = comp[["date", "fat_mass_lbs", "lean_mass_lbs", "weight_lbs"]].dropna(
        subset=["fat_mass_lbs", "lean_mass_lbs"]
    )
    comp = comp.sort_values("date").drop_duplicates("date", keep="last")

    rows = []
    out = ROOT / "analysis" / "_ax_tmp_small.csv"
    for we, wt, tf, ts in CONFIGS:
        cp = subprocess.run(
            [
                str(PY), str(AX),
                "--w-energy", str(we),
                "--w-tdee-prior", str(wt),
                "--l2-tdee-first", str(tf),
                "--l2-tdee-second", str(ts),
                "--out-file", str(out),
            ],
            capture_output=True,
            text=True,
        )
        if cp.returncode != 0:
            rows.append({
                "w_energy": we, "w_tdee_prior": wt, "l2_tdee_first": tf, "l2_tdee_second": ts,
                "status": "fail",
            })
            continue

        df = pd.read_csv(out, parse_dates=["date"])
        wsub = df.merge(weight, on="date", how="inner", suffixes=("", "_obs"))
        asub = df.merge(comp, on="date", how="inner")
        implied = df["calories"].values[:-1] - 3500 * np.diff(df["fm_lbs_weak"].values)

        rows.append({
            "w_energy": we,
            "w_tdee_prior": wt,
            "l2_tdee_first": tf,
            "l2_tdee_second": ts,
            "weight_rmse": float(np.sqrt(np.mean((wsub["weight_pred_lbs"] - wsub["smoothed_weight_lbs"]) ** 2))),
            "fm_rmse": float(np.sqrt(np.mean((asub["fm_lbs_weak"] - asub["fat_mass_lbs"]) ** 2))),
            "ffm_rmse": float(np.sqrt(np.mean((asub["ffm_lbs_weak"] - asub["lean_mass_lbs"]) ** 2))),
            "tdee_sd": float(np.std(implied)),
            "tdee_p1": float(np.percentile(implied, 1)),
            "tdee_p99": float(np.percentile(implied, 99)),
            "status": "ok",
        })

    res = pd.DataFrame(rows)
    print(res.to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
