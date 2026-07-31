"""Describe observed tirzepatide-era outcomes and identification limits.

This script uses the injection log, raw daily intake, and observed scale
weights. It does not use modeled blood levels, effective levels, P2/P3/P4,
set-point states, or an assumed calorie undercount multiplier.

Run:

    uv run --with numpy --with pandas python analysis/ED_tirzepatide_direct_audit.py

Outputs:
    analysis/ED_tirzepatide_phase_summary.csv
    analysis/ED_tirzepatide_window_summary.csv
    analysis/ED_tirzepatide_cycle_summary.csv
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
HIGH_INTAKE_KCAL = 3000
LOW_INTAKE_KCAL = 1500
PRE_WINDOW_DAYS = 365
EARLY_LATE_WINDOW_DAYS = 90


def load_injections() -> pd.DataFrame:
    medicine = pd.read_csv(ROOT / "drugs" / "medicine.csv", parse_dates=["date"])
    injections = medicine[
        medicine["drug"].str.contains(
            "tirzepatide|zepbound",
            case=False,
            na=False,
        )
    ].copy()
    injections["dose_mg"] = pd.to_numeric(injections["dose_mg"], errors="coerce")
    injections["clicks"] = pd.to_numeric(injections["clicks"], errors="coerce")
    return injections.sort_values("date").reset_index(drop=True)


def assign_injection_context(
    intake: pd.DataFrame,
    injections: pd.DataFrame,
) -> pd.DataFrame:
    context = injections[["date", "dose_mg", "clicks"]].rename(
        columns={"date": "injection_date"}
    )
    daily = pd.merge_asof(
        intake.sort_values("date"),
        context.sort_values("injection_date"),
        left_on="date",
        right_on="injection_date",
        direction="backward",
    )
    daily["on_drug"] = daily["injection_date"].notna()
    daily["days_since_injection"] = (
        daily["date"] - daily["injection_date"]
    ).dt.days
    daily["weekday"] = daily["date"].dt.dayofweek
    daily["high_intake"] = daily["calories"] >= HIGH_INTAKE_KCAL
    daily["low_intake"] = daily["calories"] <= LOW_INTAKE_KCAL
    return daily


def observed_weight_summary(
    weights: pd.DataFrame,
    start: pd.Timestamp,
    end: pd.Timestamp,
) -> dict[str, object]:
    frame = weights[weights["date"].between(start, end)].sort_values("date")
    if frame.empty:
        return {
            "weight_observations": 0,
            "first_weight_date": "",
            "first_weight_lbs": np.nan,
            "last_weight_date": "",
            "last_weight_lbs": np.nan,
            "observed_weight_change_lbs": np.nan,
        }
    first = frame.iloc[0]
    last = frame.iloc[-1]
    return {
        "weight_observations": len(frame),
        "first_weight_date": first["date"].strftime("%Y-%m-%d"),
        "first_weight_lbs": float(first["weight_lbs"]),
        "last_weight_date": last["date"].strftime("%Y-%m-%d"),
        "last_weight_lbs": float(last["weight_lbs"]),
        "observed_weight_change_lbs": float(
            last["weight_lbs"] - first["weight_lbs"]
        ),
    }


def summarize_frame(
    label: str,
    frame: pd.DataFrame,
    weights: pd.DataFrame,
) -> dict[str, object]:
    start = frame["date"].min()
    end = frame["date"].max()
    row = {
        "period": label,
        "start": start.strftime("%Y-%m-%d"),
        "end": end.strftime("%Y-%m-%d"),
        "days": len(frame),
        "mean_calories": frame["calories"].mean(),
        "median_calories": frame["calories"].median(),
        "calorie_sd": frame["calories"].std(),
        "high_intake_rate": frame["high_intake"].mean(),
        "low_intake_rate": frame["low_intake"].mean(),
        "mean_protein_g": frame["protein_g"].mean(),
        "mean_fiber_g": frame["fiber_g"].mean(),
    }
    row.update(observed_weight_summary(weights, start, end))
    return row


def dose_phase_summary(
    daily: pd.DataFrame,
    weights: pd.DataFrame,
) -> pd.DataFrame:
    on_drug = daily[daily["on_drug"]].copy()
    rows = []
    for dose, frame in on_drug.groupby("dose_mg", sort=True):
        row = summarize_frame(f"{dose:g} mg assigned phase", frame, weights)
        row["dose_mg"] = dose
        rows.append(row)
    return pd.DataFrame(rows)


def comparison_window_summary(
    daily: pd.DataFrame,
    weights: pd.DataFrame,
    drug_start: pd.Timestamp,
) -> pd.DataFrame:
    intake_end = daily["date"].max()
    windows = [
        (
            "pre-drug 365 days",
            drug_start - pd.Timedelta(days=PRE_WINDOW_DAYS),
            drug_start - pd.Timedelta(days=1),
        ),
        (
            "pre-drug 90 days",
            drug_start - pd.Timedelta(days=EARLY_LATE_WINDOW_DAYS),
            drug_start - pd.Timedelta(days=1),
        ),
        (
            "first 90 on-drug days",
            drug_start,
            drug_start + pd.Timedelta(days=EARLY_LATE_WINDOW_DAYS - 1),
        ),
        (
            "last 90 logged on-drug days",
            intake_end - pd.Timedelta(days=EARLY_LATE_WINDOW_DAYS - 1),
            intake_end,
        ),
    ]
    rows = []
    for label, start, end in windows:
        frame = daily[daily["date"].between(start, end)].copy()
        rows.append(summarize_frame(label, frame, weights))
    return pd.DataFrame(rows)


def cycle_summary(
    daily: pd.DataFrame,
    drug_start: pd.Timestamp,
) -> tuple[pd.DataFrame, float, pd.Series]:
    on_drug = daily[
        daily["on_drug"] & daily["days_since_injection"].between(0, 6)
    ].copy()
    pre = daily[
        daily["date"].between(
            drug_start - pd.Timedelta(days=PRE_WINDOW_DAYS),
            drug_start - pd.Timedelta(days=1),
        )
    ].copy()
    pre_weekday = pre.groupby("weekday")["calories"].mean()

    rows = []
    for cycle_day, frame in on_drug.groupby("days_since_injection", sort=True):
        weekday_counts = frame["weekday"].value_counts()
        dominant_weekday = int(weekday_counts.index[0])
        rows.append(
            {
                "cycle_day": int(cycle_day),
                "days": len(frame),
                "dominant_weekday": dominant_weekday,
                "dominant_weekday_share": weekday_counts.iloc[0] / len(frame),
                "mean_calories": frame["calories"].mean(),
                "median_calories": frame["calories"].median(),
                "high_intake_rate": frame["high_intake"].mean(),
                "pre365_same_weekday_mean_calories": pre_weekday.loc[
                    dominant_weekday
                ],
                "change_vs_pre_same_weekday": (
                    frame["calories"].mean() - pre_weekday.loc[dominant_weekday]
                ),
            }
        )

    cross = pd.crosstab(
        on_drug["days_since_injection"],
        on_drug["weekday"],
    )
    mapping_purity = float(cross.max(axis=1).sum() / cross.to_numpy().sum())
    return pd.DataFrame(rows), mapping_purity, on_drug["weekday"].value_counts()


def print_report(
    injections: pd.DataFrame,
    phases: pd.DataFrame,
    windows: pd.DataFrame,
    cycles: pd.DataFrame,
    mapping_purity: float,
) -> None:
    drug_start = injections["date"].min()
    q3d_start = injections.loc[injections["clicks"].notna(), "date"].min()
    weekly = injections[
        injections["clicks"].isna()
        & (injections["date"] < q3d_start if pd.notna(q3d_start) else True)
    ]
    injection_weekdays = weekly["date"].dt.day_name().value_counts()

    print("DIRECT TIRZEPATIDE AUDIT")
    print("=" * 112)
    print(
        f"Injections logged: {len(injections)} | first: {drug_start.date()} | "
        f"raw intake ends: {windows['end'].max()}"
    )
    print("\nObserved comparison windows")
    print(
        windows[
            [
                "period",
                "days",
                "mean_calories",
                "median_calories",
                "high_intake_rate",
                "mean_protein_g",
                "first_weight_lbs",
                "last_weight_lbs",
                "observed_weight_change_lbs",
            ]
        ].round(2).to_string(index=False)
    )

    print("\nAssigned dose phases: descriptive only; dose rises with calendar time")
    print(
        phases[
            [
                "period",
                "days",
                "mean_calories",
                "high_intake_rate",
                "first_weight_lbs",
                "last_weight_lbs",
                "observed_weight_change_lbs",
            ]
        ].round(2).to_string(index=False)
    )

    print("\nRaw weekly-cycle table")
    print(cycles.round(3).to_string(index=False))
    print(
        f"\nCycle-day to weekday mapping purity: {mapping_purity:.1%}. "
        "A value near 100% means cycle day and weekday are nearly the same variable."
    )
    print("Weekly injection weekdays:")
    print(injection_weekdays.to_string())
    if pd.notna(q3d_start):
        print(
            f"\nClick/q3d entries begin {q3d_start.date()}. "
            "They occur after the raw intake endpoint and therefore have no "
            "canonical intake outcome yet."
        )


def main() -> None:
    injections = load_injections()
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    weights = pd.read_csv(ROOT / "weight" / "weight.csv", parse_dates=["date"])
    daily = assign_injection_context(intake, injections)
    drug_start = injections["date"].min()

    phases = dose_phase_summary(daily, weights)
    windows = comparison_window_summary(daily, weights, drug_start)
    cycles, mapping_purity, _ = cycle_summary(daily, drug_start)

    phases.to_csv(
        ROOT / "analysis" / "ED_tirzepatide_phase_summary.csv",
        index=False,
    )
    windows.to_csv(
        ROOT / "analysis" / "ED_tirzepatide_window_summary.csv",
        index=False,
    )
    cycles.to_csv(
        ROOT / "analysis" / "ED_tirzepatide_cycle_summary.csv",
        index=False,
    )
    print_report(injections, phases, windows, cycles, mapping_purity)


if __name__ == "__main__":
    main()
