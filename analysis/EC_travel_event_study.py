"""Paired event study of travel and intake using raw observations only.

Trips come from travel/trips.md. The unavailable underlying Google Timeline
export is not reconstructed here. Each trip is paired with the 14 calendar
days before and after it. P2/P3/P4 and modeled drug exposure are not used.

Run:

    uv run --with numpy --with pandas python analysis/EC_travel_event_study.py

Outputs:
    analysis/EC_travel_trip_results.csv
    analysis/EC_travel_summary.csv
"""

from __future__ import annotations

from pathlib import Path
import math

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
COMPARISON_DAYS = 14
MIN_TRIP_DAYS = 2
MIN_BASELINE_DAYS = 10
HIGH_INTAKE_KCAL = 3000
TIRZEPATIDE_START = pd.Timestamp("2024-09-17")


def parse_trips() -> pd.DataFrame:
    rows = []
    for line in (ROOT / "travel" / "trips.md").read_text().splitlines():
        parts = [part.strip() for part in line.strip().strip("|").split("|")]
        if len(parts) < 6:
            continue
        try:
            start = pd.Timestamp(parts[0])
            end = pd.Timestamp(parts[1])
        except (ValueError, TypeError):
            continue
        rows.append(
            {
                "start": start,
                "end": end,
                "listed_days": int(parts[2]),
                "location": parts[3],
                "distance": parts[4],
                "notes": parts[5],
            }
        )
    return pd.DataFrame(rows).sort_values("start").reset_index(drop=True)


def nearest_weight_change(
    weights: pd.DataFrame,
    start: pd.Timestamp,
    end: pd.Timestamp,
) -> tuple[float, str, str]:
    before = weights[
        weights["date"].between(
            start - pd.Timedelta(days=COMPARISON_DAYS),
            start - pd.Timedelta(days=1),
        )
    ].sort_values("date")
    after = weights[
        weights["date"].between(
            end + pd.Timedelta(days=1),
            end + pd.Timedelta(days=COMPARISON_DAYS),
        )
    ].sort_values("date")
    if before.empty or after.empty:
        return np.nan, "", ""
    before_row = before.iloc[-1]
    after_row = after.iloc[0]
    return (
        float(after_row["weight_lbs"] - before_row["weight_lbs"]),
        before_row["date"].strftime("%Y-%m-%d"),
        after_row["date"].strftime("%Y-%m-%d"),
    )


def build_trip_results() -> pd.DataFrame:
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    weights = pd.read_csv(ROOT / "weight" / "weight.csv", parse_dates=["date"])
    trips = parse_trips()
    intake_start = intake["date"].min()
    intake_end = intake["date"].max()

    rows = []
    for trip in trips.itertuples():
        trip_days = (trip.end - trip.start).days + 1
        if (
            trip_days < MIN_TRIP_DAYS
            or trip.start < intake_start
            or trip.end > intake_end
        ):
            continue

        pre = intake[
            intake["date"].between(
                trip.start - pd.Timedelta(days=COMPARISON_DAYS),
                trip.start - pd.Timedelta(days=1),
            )
        ]
        during = intake[intake["date"].between(trip.start, trip.end)]
        post = intake[
            intake["date"].between(
                trip.end + pd.Timedelta(days=1),
                trip.end + pd.Timedelta(days=COMPARISON_DAYS),
            )
        ]
        if len(pre) < MIN_BASELINE_DAYS or len(post) < MIN_BASELINE_DAYS:
            continue

        weight_delta, weight_before_date, weight_after_date = nearest_weight_change(
            weights,
            trip.start,
            trip.end,
        )
        rows.append(
            {
                "start": trip.start,
                "end": trip.end,
                "trip_days": trip_days,
                "location": trip.location,
                "notes": trip.notes,
                "drug_era": (
                    "tirzepatide"
                    if trip.start >= TIRZEPATIDE_START
                    else "pre_tirzepatide"
                ),
                "pre_days": len(pre),
                "during_days": len(during),
                "post_days": len(post),
                "pre_mean_calories": pre["calories"].mean(),
                "during_mean_calories": during["calories"].mean(),
                "post_mean_calories": post["calories"].mean(),
                "during_minus_pre": (
                    during["calories"].mean() - pre["calories"].mean()
                ),
                "post_minus_pre": (
                    post["calories"].mean() - pre["calories"].mean()
                ),
                "pre_high_intake_rate": (
                    pre["calories"] >= HIGH_INTAKE_KCAL
                ).mean(),
                "during_high_intake_rate": (
                    during["calories"] >= HIGH_INTAKE_KCAL
                ).mean(),
                "post_high_intake_rate": (
                    post["calories"] >= HIGH_INTAKE_KCAL
                ).mean(),
                "observed_weight_change_lbs": weight_delta,
                "weight_before_date": weight_before_date,
                "weight_after_date": weight_after_date,
            }
        )
    return pd.DataFrame(rows)


def exact_sign_p(values: pd.Series) -> float:
    nonzero = values.dropna()
    nonzero = nonzero[nonzero != 0]
    n = len(nonzero)
    if n == 0:
        return np.nan
    positives = int((nonzero > 0).sum())
    lower_tail = min(positives, n - positives)
    probability = sum(math.comb(n, count) for count in range(lower_tail + 1)) / (
        2**n
    )
    return min(1.0, 2 * probability)


def summarize_group(label: str, frame: pd.DataFrame) -> dict[str, object]:
    delta = frame["during_minus_pre"]
    post_delta = frame["post_minus_pre"]
    return {
        "group": label,
        "trips": len(frame),
        "travel_days": int(frame["during_days"].sum()),
        "mean_during_minus_pre": delta.mean(),
        "median_during_minus_pre": delta.median(),
        "trips_higher_than_pre": int((delta > 0).sum()),
        "trips_lower_than_pre": int((delta < 0).sum()),
        "during_sign_test_p": exact_sign_p(delta),
        "mean_post_minus_pre": post_delta.mean(),
        "median_post_minus_pre": post_delta.median(),
        "pre_high_intake_rate_unweighted": frame["pre_high_intake_rate"].mean(),
        "during_high_intake_rate_unweighted": frame[
            "during_high_intake_rate"
        ].mean(),
        "post_high_intake_rate_unweighted": frame["post_high_intake_rate"].mean(),
        "trips_with_weight_anchors": int(
            frame["observed_weight_change_lbs"].notna().sum()
        ),
        "median_observed_weight_change_lbs": frame[
            "observed_weight_change_lbs"
        ].median(),
    }


def build_summary(trips: pd.DataFrame) -> pd.DataFrame:
    groups = [
        ("all eligible trips", trips),
        (
            "pre-tirzepatide",
            trips[trips["drug_era"] == "pre_tirzepatide"],
        ),
        (
            "pre-tirzepatide excluding Jordan",
            trips[
                (trips["drug_era"] == "pre_tirzepatide")
                & ~trips["location"].str.contains("Jordan", case=False, na=False)
            ],
        ),
        (
            "short trips 2-6 days",
            trips[trips["trip_days"].between(2, 6)],
        ),
        (
            "long trips 7+ days",
            trips[trips["trip_days"] >= 7],
        ),
        (
            "modern MHTML pre-tirzepatide",
            trips[
                (trips["start"] >= pd.Timestamp("2018-03-01"))
                & (trips["drug_era"] == "pre_tirzepatide")
            ],
        ),
        (
            "tirzepatide-era",
            trips[trips["drug_era"] == "tirzepatide"],
        ),
    ]
    return pd.DataFrame(
        summarize_group(label, frame)
        for label, frame in groups
        if not frame.empty
    )


def print_report(trips: pd.DataFrame, summary: pd.DataFrame) -> None:
    print("TRAVEL EVENT STUDY")
    print("=" * 112)
    print(
        "Each trip is paired with the 14 days before and after it. "
        "Only raw intake and observed scale weights are used."
    )
    print("\nPaired summaries")
    print(
        summary[
            [
                "group",
                "trips",
                "travel_days",
                "mean_during_minus_pre",
                "median_during_minus_pre",
                "trips_higher_than_pre",
                "trips_lower_than_pre",
                "during_sign_test_p",
                "mean_post_minus_pre",
            ]
        ].round(2).to_string(index=False)
    )

    print("\nTrips sorted by intake change versus their own prior 14 days")
    display = trips[
        [
            "start",
            "location",
            "trip_days",
            "drug_era",
            "pre_mean_calories",
            "during_mean_calories",
            "during_minus_pre",
            "post_minus_pre",
            "observed_weight_change_lbs",
        ]
    ].sort_values("during_minus_pre", ascending=False)
    numeric_columns = display.select_dtypes(include=[np.number]).columns
    display = display.copy()
    display[numeric_columns] = display[numeric_columns].round(1)
    print(
        display.to_string(index=False)
    )


def main() -> None:
    trips = build_trip_results()
    summary = build_summary(trips)
    trips.to_csv(ROOT / "analysis" / "EC_travel_trip_results.csv", index=False)
    summary.to_csv(ROOT / "analysis" / "EC_travel_summary.csv", index=False)
    print_report(trips, summary)


if __name__ == "__main__":
    main()
