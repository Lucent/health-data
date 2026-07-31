"""Audit same-device composition, protein, and strength during tirzepatide.

This is descriptive. It uses observed H30 scans, raw intake, and logged
strength-session dates. It does not use interpolated or Kalman states.

Run:

    uv run --with numpy --with pandas python analysis/EE_tirzepatide_composition_support.py

Outputs:
    analysis/EE_tirzepatide_composition_summary.csv
    analysis/EE_tirzepatide_support_periods.csv
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
DRUG_START = pd.Timestamp("2024-09-17")
INTAKE_END = pd.Timestamp("2026-02-28")
H30_DEVICE = "H30"
WINDOW_DAYS = 90


def metric_summary(scans: pd.DataFrame) -> pd.DataFrame:
    first = scans.iloc[0]
    last = scans.iloc[-1]
    rows = []
    for metric in [
        "weight_lbs",
        "fat_mass_lbs",
        "lean_mass_lbs",
        "smm_lbs",
        "visceral_fat_level",
        "waist_hip_ratio",
    ]:
        if metric not in scans or scans[metric].notna().sum() == 0:
            continue
        rows.append(
            {
                "device": H30_DEVICE,
                "metric": metric,
                "first_date": first["date"].strftime("%Y-%m-%d"),
                "last_date": last["date"].strftime("%Y-%m-%d"),
                "scans": int(scans[metric].notna().sum()),
                "first_value": first[metric],
                "last_value": last[metric],
                "change": last[metric] - first[metric],
                "minimum": scans[metric].min(),
                "maximum": scans[metric].max(),
            }
        )
    return pd.DataFrame(rows)


def support_period_summary(
    intake: pd.DataFrame,
    strength: pd.DataFrame,
    label: str,
    start: pd.Timestamp,
    end: pd.Timestamp,
    weight_reference_lbs: float,
    lean_reference_lbs: float,
) -> dict[str, object]:
    diet = intake[intake["date"].between(start, end)].copy()
    sessions = strength[strength["date"].between(start, end)].copy()
    calendar_days = (end - start).days + 1
    mean_protein = diet["protein_g"].mean()
    return {
        "period": label,
        "start": start.strftime("%Y-%m-%d"),
        "end": end.strftime("%Y-%m-%d"),
        "calendar_days": calendar_days,
        "logged_intake_days": len(diet),
        "mean_calories": diet["calories"].mean(),
        "mean_protein_g": mean_protein,
        "median_protein_g": diet["protein_g"].median(),
        "protein_below_60g_days": int((diet["protein_g"] < 60).sum()),
        "protein_below_80g_days": int((diet["protein_g"] < 80).sum()),
        "protein_at_least_80g_days": int((diet["protein_g"] >= 80).sum()),
        "mean_protein_g_per_kg_total_weight_at_end": (
            mean_protein / (weight_reference_lbs * 0.453592)
        ),
        "mean_protein_g_per_kg_lean_reference": (
            mean_protein / (lean_reference_lbs * 0.453592)
        ),
        "strength_sessions": len(sessions),
        "strength_sessions_per_week": len(sessions) / calendar_days * 7,
        "strength_duration_field_minutes": sessions["duration_min"].sum(),
    }


def main() -> None:
    composition = pd.read_csv(
        ROOT / "composition" / "composition.csv",
        parse_dates=["date"],
    )
    intake = pd.read_csv(
        ROOT / "intake" / "intake_daily.csv",
        parse_dates=["date"],
    )
    strength = pd.read_csv(
        ROOT / "workout" / "strength.csv",
        parse_dates=["date"],
    )
    h30 = composition[composition["device"] == H30_DEVICE].sort_values("date")
    metrics = metric_summary(h30)

    first = h30.iloc[0]
    last = h30.iloc[-1]
    total_loss = first["weight_lbs"] - last["weight_lbs"]
    fat_loss = first["fat_mass_lbs"] - last["fat_mass_lbs"]
    measured_fat_share = fat_loss / total_loss if total_loss > 0 else np.nan
    lean_reference_lbs = float(last["lean_mass_lbs"])
    weight_reference_lbs = float(last["weight_lbs"])

    periods = [
        (
            "full logged tirzepatide era",
            DRUG_START,
            INTAKE_END,
        ),
        (
            "first 90 logged tirzepatide days",
            DRUG_START,
            DRUG_START + pd.Timedelta(days=WINDOW_DAYS - 1),
        ),
        (
            "last 90 logged tirzepatide days",
            INTAKE_END - pd.Timedelta(days=WINDOW_DAYS - 1),
            INTAKE_END,
        ),
        (
            "same-device H30 interval",
            first["date"],
            last["date"],
        ),
    ]
    support = pd.DataFrame(
        support_period_summary(
            intake,
            strength,
            label,
            start,
            end,
            weight_reference_lbs,
            lean_reference_lbs,
        )
        for label, start, end in periods
    )

    metrics.to_csv(
        ROOT / "analysis" / "EE_tirzepatide_composition_summary.csv",
        index=False,
    )
    support.to_csv(
        ROOT / "analysis" / "EE_tirzepatide_support_periods.csv",
        index=False,
    )

    print("TIRZEPATIDE COMPOSITION / SUPPORT AUDIT")
    print("=" * 112)
    print(
        f"H30 scans: {len(h30)}, {first['date'].date()} to {last['date'].date()}. "
        f"The first H30 scan is {(first['date'] - DRUG_START).days} days after "
        "tirzepatide began."
    )
    print("\nSame-device observed changes")
    print(metrics.round(3).to_string(index=False))
    print(
        f"\nMeasured fat loss / measured total loss: {measured_fat_share:.1%}. "
        "Values above 100% can occur when BIA lean/water estimates rise."
    )
    print("\nDiet and logged training support")
    print(
        support[
            [
                "period",
                "calendar_days",
                "mean_calories",
                "mean_protein_g",
                "protein_below_60g_days",
                "protein_at_least_80g_days",
                "mean_protein_g_per_kg_lean_reference",
                "strength_sessions",
                "strength_sessions_per_week",
            ]
        ].round(3).to_string(index=False)
    )
    print(
        "\nThe strength file records session dates only and assigns every session "
        "30 minutes; it cannot establish progressive overload or training quality."
    )


if __name__ == "__main__":
    main()
