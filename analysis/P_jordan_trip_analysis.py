"""Summarize the two month-long Jordan trips for unexpected TDEE or weight swing."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.append(str(ROOT))

import pandas as pd

import O_diet_epoch_analysis as diet_epoch_analysis
PRE_WINDOW_DAYS = 28
POST_WINDOW_DAYS = 28


@dataclass
class Trip:
    label: str
    start: pd.Timestamp
    end: pd.Timestamp


TRIPS = [
    Trip("Jordan 2015", pd.Timestamp("2015-09-07"), pd.Timestamp("2015-10-07")),
    Trip("Jordan 2019", pd.Timestamp("2019-06-02"), pd.Timestamp("2019-06-30")),
]


def load_daily_with_weight() -> pd.DataFrame:
    daily = diet_epoch_analysis.load_daily()
    weight = pd.read_csv(ROOT / "analysis" / "P1_smoothed_weight.csv", parse_dates=["date"])
    daily = daily.merge(weight[["date", "smoothed_weight_lbs"]], on="date", how="left")
    return daily


def summarize_segment(df: pd.DataFrame) -> dict[str, float]:
    summary = diet_epoch_analysis.summarize_slice(df)
    summary["mean_tdee"] = round(df["tdee"].mean(), 1)
    if df["smoothed_weight_lbs"].notna().any():
        start_weight = df["smoothed_weight_lbs"].ffill().iloc[0]
        end_weight = df["smoothed_weight_lbs"].ffill().iloc[-1]
        summary.update(
            {
                "weight_start_lbs": round(start_weight, 2),
                "weight_end_lbs": round(end_weight, 2),
                "weight_delta_lbs": round(end_weight - start_weight, 2),
            }
        )
    else:
        summary.update(
            {
                "weight_start_lbs": None,
                "weight_end_lbs": None,
                "weight_delta_lbs": None,
            }
        )
    return summary


def build_windows(daily: pd.DataFrame, trip: Trip) -> dict[str, pd.DataFrame]:
    pre = daily[
        (daily["date"] >= trip.start - pd.Timedelta(days=PRE_WINDOW_DAYS))
        & (daily["date"] < trip.start)
        & (daily["effective_level"] == 0)
    ]
    during = daily[
        (daily["date"] >= trip.start)
        & (daily["date"] <= trip.end)
        & (daily["effective_level"] == 0)
    ]
    post = daily[
        (daily["date"] > trip.end)
        & (daily["date"] <= trip.end + pd.Timedelta(days=POST_WINDOW_DAYS))
        & (daily["effective_level"] == 0)
    ]
    return {"pre": pre, "during": during, "post": post}


def summarize_trips(daily: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for trip in TRIPS:
        windows = build_windows(daily, trip)
        for phase, df in windows.items():
            if df.empty:
                continue
            row = {"trip": trip.label, "phase": phase}
            row.update(summarize_segment(df))
            rows.append(row)
    return pd.DataFrame(rows)


def summarize_deltas(summary_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for trip in TRIPS:
        trip_rows = summary_df[summary_df["trip"] == trip.label]
        phases = {phase: row for phase, row in zip(trip_rows["phase"], trip_rows.to_dict("records"))}
        if not {"pre", "during", "post"} <= set(phases):
            continue
        trip_summary = phases["during"]
        pre_summary = phases["pre"]
        post_summary = phases["post"]
        rows.append(
            {
                "trip": trip.label,
                "during_minus_pre_tdee": round(
                    trip_summary["mean_tdee"] - pre_summary["mean_tdee"], 1
                ),
                "during_minus_pre_ratio": round(
                    trip_summary["mean_tdee_rmr_ratio"] - pre_summary["mean_tdee_rmr_ratio"], 4
                ),
                "weight_delta_during": trip_summary["weight_delta_lbs"],
                "fat_delta_during": trip_summary["fat_delta_lbs"],
                "post_minus_during_ratio": round(
                    post_summary["mean_tdee_rmr_ratio"] - trip_summary["mean_tdee_rmr_ratio"], 4
                ),
            }
        )
    return pd.DataFrame(rows)


def save_outputs(summary: pd.DataFrame, deltas: pd.DataFrame) -> None:
    summary.to_csv(ROOT / "analysis" / "P_jordan_trip_summary.csv", index=False)
    deltas.to_csv(ROOT / "analysis" / "P_jordan_trip_delta.csv", index=False)


def print_report(summary: pd.DataFrame, deltas: pd.DataFrame) -> None:
    print("\n=== Jordan Trip Summary ===")
    print(summary[["trip", "phase", "days", "mean_calories", "mean_tdee", "mean_tdee_rmr_ratio", "weight_delta_lbs", "fat_delta_lbs"]].to_string(index=False))
    print("\nDelta highlights:")
    print(deltas.to_string(index=False))


def main() -> None:
    daily = load_daily_with_weight()
    summary = summarize_trips(daily)
    deltas = summarize_deltas(summary)
    save_outputs(summary, deltas)
    print_report(summary, deltas)


if __name__ == "__main__":
    main()
