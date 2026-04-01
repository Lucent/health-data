#!/usr/bin/env python3
"""Legacy compatibility wrapper for fixed-outcome half-life comparison.

Superseded by BC_set_point_unified.py. Kept as a wrapper to the shared library.
"""

from set_point_lib import (
    ROOT,
    add_rolling_outcomes,
    add_window_features,
    build_daily_base,
    evaluate_fixed_outcome_half_life,
    load_state_daily,
)
import pandas as pd


def main():
    base = build_daily_base()
    rows = []

    weight = pd.read_csv(ROOT / "analysis" / "P1_smoothed_weight.csv", parse_dates=["date"])
    obs = weight[["date", "smoothed_weight_lbs"]].dropna().rename(columns={"smoothed_weight_lbs": "state"})
    obs = add_rolling_outcomes(add_window_features(obs.merge(base, on="date", how="left"), 30))
    for outcome in ["next_30d_cal", "delta_30d_cal"]:
        res = evaluate_fixed_outcome_half_life(obs, "state", outcome, irregular=True)
        rows.append({"model": "Corrected Weight Only", "outcome": outcome, **res})

    for path, col, name in [
        (ROOT / "analysis" / "AV_anchor_daily_composition.csv", "fm_lbs_anchor", "Pure Anchor FM"),
        (ROOT / "analysis" / "AX_weakintake_daily_composition.csv", "fm_lbs_weak", "Weak-Intake Anchor FM"),
        (ROOT / "analysis" / "P4_kalman_daily.csv", "fat_mass_lbs", "Full Kalman FM"),
    ]:
        daily = add_rolling_outcomes(add_window_features(load_state_daily(path, col), 30))
        for outcome in ["next_30d_cal", "delta_30d_cal", "mean_90d_cal"]:
            res = evaluate_fixed_outcome_half_life(daily, "state", outcome, irregular=False)
            rows.append({"model": name, "outcome": outcome, **res})

    out = pd.DataFrame(rows)
    print("=" * 70)
    print("HALF-LIFE WITH FIXED OUTCOMES")
    print("=" * 70)
    print(out.to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
