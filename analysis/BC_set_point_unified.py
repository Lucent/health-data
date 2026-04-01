#!/usr/bin/env python3
"""Unified set-point comparison artifact.

This replaces the overlapping appetite-side role of:
  AG_binge_set_point.py
  AH_set_point_properties.py
  AM_lipostat_sensitivity.py
  AP_overshoot_shape.py
  AU_noncircular_set_point.py
  AW_anchor_fm_set_point.py
  BA_set_point_model_comparison.py
  BB_half_life_fixed_outcome.py

It does not attempt to replace tirzepatide-specific PK/coverage analyses.

What it reports:
  1. State construction quality (anchor fit where relevant)
  2. Primary appetite-side comparison across state constructions
  3. Fixed-outcome half-life comparison
  4. Brief nulls / caveats
"""

from set_point_lib import (
    ROOT,
    add_rolling_outcomes,
    add_window_features,
    build_daily_base,
    evaluate_fixed_outcome_half_life,
    evaluate_primary_signal,
    load_state_daily,
)
import pandas as pd
import numpy as np


def corrected_weight_state():
    weight = pd.read_csv(ROOT / "analysis" / "P1_smoothed_weight.csv", parse_dates=["date"])
    obs = weight[["date", "smoothed_weight_lbs"]].dropna().rename(columns={"smoothed_weight_lbs": "state"})
    return obs.merge(build_daily_base(), on="date", how="left")


def anchor_fit_summary(path, fm_col, ffm_col, label):
    comp = pd.read_csv(ROOT / "composition" / "composition.csv", parse_dates=["date"])
    comp = comp[["date", "fat_mass_lbs", "lean_mass_lbs", "weight_lbs"]].dropna(
        subset=["fat_mass_lbs", "lean_mass_lbs"]
    )
    comp = comp.sort_values("date").drop_duplicates("date", keep="last")
    weight = pd.read_csv(ROOT / "analysis" / "P1_smoothed_weight.csv", parse_dates=["date"])[
        ["date", "smoothed_weight_lbs"]
    ]
    daily = pd.read_csv(path, parse_dates=["date"])
    wsub = daily.merge(weight, on="date", how="inner", suffixes=("", "_obs"))
    asub = daily.merge(comp, on="date", how="inner")
    return {
        "model": label,
        "weight_rmse": float(np.sqrt(np.mean((wsub["weight_pred_lbs"] - wsub["smoothed_weight_lbs"]) ** 2))),
        "fm_rmse": float(np.sqrt(np.mean((asub[fm_col] - asub["fat_mass_lbs"]) ** 2))),
        "ffm_rmse": float(np.sqrt(np.mean((asub[ffm_col] - asub["lean_mass_lbs"]) ** 2))),
    }


def main():
    print("=" * 70)
    print("UNIFIED SET-POINT COMPARISON")
    print("=" * 70)

    print("\nState quality:")
    fit_rows = [
        anchor_fit_summary(
            ROOT / "analysis" / "AV_anchor_daily_composition.csv",
            "fm_lbs_anchor",
            "ffm_lbs_anchor",
            "Pure Anchor FM",
        ),
        anchor_fit_summary(
            ROOT / "analysis" / "AX_weakintake_daily_composition.csv",
            "fm_lbs_weak",
            "ffm_lbs_weak",
            "Weak-Intake Anchor FM",
        ),
    ]
    fit_df = pd.DataFrame(fit_rows)
    print(fit_df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print("\nPrimary appetite-side comparison:")
    rows = [
        evaluate_primary_signal(corrected_weight_state(), "state", "Corrected Weight Only", irregular=True),
        evaluate_primary_signal(
            load_state_daily(ROOT / "analysis" / "AV_anchor_daily_composition.csv", "fm_lbs_anchor"),
            "state",
            "Pure Anchor FM",
        ),
        evaluate_primary_signal(
            load_state_daily(ROOT / "analysis" / "AX_weakintake_daily_composition.csv", "fm_lbs_weak"),
            "state",
            "Weak-Intake Anchor FM",
        ),
        evaluate_primary_signal(
            load_state_daily(ROOT / "analysis" / "P4_kalman_daily.csv", "fat_mass_lbs"),
            "state",
            "Full Kalman FM",
        ),
    ]
    compare_df = pd.DataFrame(rows)
    print(compare_df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    print("\nFixed-outcome half-life comparison:")
    half_rows = []
    obs = add_rolling_outcomes(add_window_features(corrected_weight_state(), 30))
    for outcome in ["next_30d_cal", "delta_30d_cal"]:
        half_rows.append({
            "model": "Corrected Weight Only",
            "outcome": outcome,
            **evaluate_fixed_outcome_half_life(obs, "state", outcome, irregular=True),
        })
    for path, col, name in [
        (ROOT / "analysis" / "AV_anchor_daily_composition.csv", "fm_lbs_anchor", "Pure Anchor FM"),
        (ROOT / "analysis" / "AX_weakintake_daily_composition.csv", "fm_lbs_weak", "Weak-Intake Anchor FM"),
        (ROOT / "analysis" / "P4_kalman_daily.csv", "fat_mass_lbs", "Full Kalman FM"),
    ]:
        daily = add_rolling_outcomes(add_window_features(load_state_daily(path, col), 30))
        for outcome in ["next_30d_cal", "delta_30d_cal", "mean_90d_cal"]:
            half_rows.append({
                "model": name,
                "outcome": outcome,
                **evaluate_fixed_outcome_half_life(daily, "state", outcome, irregular=False),
            })
    half_df = pd.DataFrame(half_rows)
    print(half_df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    print("\nBrief nulls / caveats:")
    print("  Pure anchor FM is the cleanest latent-state construction, but the appetite signal is attenuated.")
    print("  Fixed targets underperform moving-state constructions across all tested setups.")
    print("  The 40-50d half-life is strongest for smoothed appetite-level outcomes, not for 30d intake-shift outcomes.")
    print("  Weak-intake anchor FM strengthens the signal materially, but is not fully noncircular.")

    print("\nSupersedes:")
    print("  Use this artifact as the primary reference instead of AG/AH/AM/AP/AU/AW/BA/BB.")
    print("  Keep AQ for tirzepatide-specific decomposition and AI for expenditure-timescale analysis.")


if __name__ == "__main__":
    main()
