#!/usr/bin/env python3
"""Legacy compatibility wrapper for the set-point model comparison.

Superseded by BC_set_point_unified.py. Kept only as a thin wrapper so older
references still resolve while all shared logic lives in one place.
"""

from set_point_lib import ROOT, build_daily_base, evaluate_primary_signal, load_state_daily
import pandas as pd


def main():
    rows = []

    weight = pd.read_csv(ROOT / "analysis" / "P1_smoothed_weight.csv", parse_dates=["date"])
    obs = weight[["date", "smoothed_weight_lbs"]].dropna().rename(columns={"smoothed_weight_lbs": "state"})
    obs = obs.merge(build_daily_base(), on="date", how="left")
    rows.append(evaluate_primary_signal(obs, "state", "Corrected Weight Only", irregular=True))

    rows.append(
        evaluate_primary_signal(
            load_state_daily(ROOT / "analysis" / "AV_anchor_daily_composition.csv", "fm_lbs_anchor"),
            "state",
            "Pure Anchor FM",
        )
    )
    rows.append(
        evaluate_primary_signal(
            load_state_daily(ROOT / "analysis" / "AX_weakintake_daily_composition.csv", "fm_lbs_weak"),
            "state",
            "Weak-Intake Anchor FM",
        )
    )
    rows.append(
        evaluate_primary_signal(
            load_state_daily(ROOT / "analysis" / "P4_kalman_daily.csv", "fat_mass_lbs"),
            "state",
            "Full Kalman FM",
        )
    )

    out = pd.DataFrame(rows)
    print("=" * 70)
    print("SET-POINT MODEL COMPARISON")
    print("=" * 70)
    print(out.to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
