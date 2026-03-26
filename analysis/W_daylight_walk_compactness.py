#!/usr/bin/env python3
"""Compare daylight walk days to step-matched smeared-step days.

Samsung's minute-level pedometer export in the current backup only covers a
recent slice of 2026. That is enough to test the structural question:
are midday/afternoon daylight walks different from similar numbers of steps
spread thinly across a longer span?

Outputs:
  analysis/daylight_walk_compactness_daily.csv
  analysis/daylight_walk_compactness_matches.csv
  analysis/daylight_walk_compactness_summary.csv
"""

from __future__ import annotations

import csv
from collections import defaultdict
from datetime import datetime
from pathlib import Path
import subprocess
import tempfile

import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
ARCHIVE = Path("/mnt/c/Users/Lucent/OneDrive/Documents/Backup/Samsung/samsunghealth_michael_20260324183778.7z")
STEP_COUNT_MEMBER = "samsunghealth_michael_20260324183778/com.samsung.shealth.tracker.pedometer_step_count.20260324183778.csv"

DAILY_PATH = ROOT / "analysis" / "W_daylight_walk_compactness_daily.csv"
MATCHES_PATH = ROOT / "analysis" / "W_daylight_walk_compactness_matches.csv"
SUMMARY_PATH = ROOT / "analysis" / "W_daylight_walk_compactness_summary.csv"


# Walk window parameters (shared with V_exercise_walk_analysis.py)
HOUR_START = 12   # noon
HOUR_END = 19     # 7pm
MIN_YEAR = 2026   # minute-level pedometer data only available for recent slice


def load_minute_steps() -> pd.DataFrame:
    per_day: dict[str, list[tuple[datetime, float]]] = defaultdict(list)
    with tempfile.TemporaryDirectory() as tmpdir:
        subprocess.run(
            ["7z", "x", f"-o{tmpdir}", str(ARCHIVE), STEP_COUNT_MEMBER],
            capture_output=True,
            text=True,
            check=True,
        )
        path = Path(tmpdir) / STEP_COUNT_MEMBER
        with open(path, encoding="utf-8-sig") as fh:
            next(fh)
            reader = csv.DictReader(fh)
            for row in reader:
                start = row.get("com.samsung.health.step_count.start_time")
                if not start:
                    continue
                dt = datetime.strptime(start[:19], "%Y-%m-%d %H:%M:%S")
                if dt.year != MIN_YEAR:
                    continue
                hour = dt.hour + dt.minute / 60
                if not (HOUR_START <= hour < HOUR_END):
                    continue
                count = float(row.get("com.samsung.health.step_count.count") or 0)
                if count <= 0:
                    continue
                per_day[dt.strftime("%Y-%m-%d")].append((dt, count))

    rows = []
    for date, points in sorted(per_day.items()):
        points = sorted(points)
        series = pd.Series(
            [count for _, count in points],
            index=pd.DatetimeIndex([dt for dt, _ in points]),
        ).resample("1min").sum()
        total = float(series.sum())
        active_minutes = int((series > 0).sum())
        span_minutes = (series.index.max() - series.index.min()).total_seconds() / 60 + 1
        rows.append(
            {
                "date": date,
                "window_steps": round(total, 2),
                "active_minutes": active_minutes,
                "span_minutes": round(span_minutes, 2),
                "steps_per_active_min": round(total / active_minutes, 4),
                "steps_per_span_min": round(total / span_minutes, 4),
                "top30_share": round(series.rolling(30, min_periods=1).sum().max() / total, 4),
                "top60_share": round(series.rolling(60, min_periods=1).sum().max() / total, 4),
            }
        )
    return pd.DataFrame(rows)


def classify_walk_days() -> pd.DataFrame:
    exercises = pd.read_csv(ROOT / "steps-sleep" / "exercises_samsung.csv", parse_dates=["start_time", "end_time"])
    walk = exercises[
        (exercises["exercise_type"].astype(str) == "1001")
        & (exercises["start_time"].dt.year == MIN_YEAR)
    ].copy()
    walk["duration_min"] = pd.to_numeric(walk["duration_min"], errors="coerce")
    walk["date"] = walk["start_time"].dt.strftime("%Y-%m-%d")
    walk["hour"] = walk["start_time"].dt.hour + walk["start_time"].dt.minute / 60
    walk = walk[(walk["hour"] >= HOUR_START) & (walk["hour"] < HOUR_END)]

    rows = []
    for date, group in walk.sort_values(["date", "start_time"]).groupby("date"):
        group = group.reset_index(drop=True)
        regime = None
        for i in range(len(group) - 1):
            break_min = (group.loc[i + 1, "start_time"] - group.loc[i, "end_time"]).total_seconds() / 60
            if (
                20 <= group.loc[i, "duration_min"] <= 45
                and 20 <= group.loc[i + 1, "duration_min"] <= 45
                and 5 <= break_min <= 40
            ):
                regime = "paired_daylight_walk"
                break
        if regime is None and ((group["duration_min"] >= 20) & (group["duration_min"] <= 45)).any():
            regime = "single_daylight_walk"
        if regime is not None:
            rows.append({"date": date, "walk_regime": regime})
    return pd.DataFrame(rows).drop_duplicates()


def build_matches(daily: pd.DataFrame) -> pd.DataFrame:
    singles = daily[daily["walk_regime"] == "single_daylight_walk"].copy()
    smeared_pool = daily[daily["walk_regime"].isna()].copy()

    matches = []
    for _, row in singles.iterrows():
        candidates = smeared_pool[
            smeared_pool["window_steps"].between(row["window_steps"] * 0.7, row["window_steps"] * 1.3)
        ].copy()
        candidates = candidates[
            (candidates["active_minutes"] >= row["active_minutes"] + 5)
            | (candidates["span_minutes"] >= row["span_minutes"] + 45)
            | (candidates["steps_per_span_min"] <= row["steps_per_span_min"] * 0.67)
        ]
        if candidates.empty:
            continue
        score = (
            (candidates["window_steps"] - row["window_steps"]).abs() / max(row["window_steps"], 1)
            + (candidates["steps_per_span_min"] - row["steps_per_span_min"]).abs() / max(row["steps_per_span_min"], 1e-6)
        )
        match = candidates.loc[score.idxmin()]
        matches.append(
            {
                "walk_date": row["date"],
                "smeared_date": match["date"],
                "walk_steps": row["window_steps"],
                "smeared_steps": match["window_steps"],
                "step_diff": round(abs(row["window_steps"] - match["window_steps"]), 2),
                "walk_active_minutes": row["active_minutes"],
                "smeared_active_minutes": match["active_minutes"],
                "walk_span_minutes": row["span_minutes"],
                "smeared_span_minutes": match["span_minutes"],
                "walk_steps_per_span_min": row["steps_per_span_min"],
                "smeared_steps_per_span_min": match["steps_per_span_min"],
                "walk_top60_share": row["top60_share"],
                "smeared_top60_share": match["top60_share"],
            }
        )
    return pd.DataFrame(matches)


def summarize(daily: pd.DataFrame, matches: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for regime in ["single_daylight_walk", "paired_daylight_walk"]:
        group = daily[daily["walk_regime"] == regime]
        rows.append(
            {
                "group": regime,
                "days": int(len(group)),
                "mean_window_steps": round(group["window_steps"].mean(), 2) if len(group) else "",
                "mean_active_minutes": round(group["active_minutes"].mean(), 2) if len(group) else "",
                "mean_span_minutes": round(group["span_minutes"].mean(), 2) if len(group) else "",
                "mean_steps_per_span_min": round(group["steps_per_span_min"].mean(), 4) if len(group) else "",
                "mean_top60_share": round(group["top60_share"].mean(), 4) if len(group) else "",
            }
        )
    if not matches.empty:
        rows.append(
            {
                "group": "single_walk_minus_smeared_match",
                "days": int(len(matches)),
                "mean_window_steps": round((matches["walk_steps"] - matches["smeared_steps"]).mean(), 2),
                "mean_active_minutes": round((matches["walk_active_minutes"] - matches["smeared_active_minutes"]).mean(), 2),
                "mean_span_minutes": round((matches["walk_span_minutes"] - matches["smeared_span_minutes"]).mean(), 2),
                "mean_steps_per_span_min": round(
                    (matches["walk_steps_per_span_min"] - matches["smeared_steps_per_span_min"]).mean(), 4
                ),
                "mean_top60_share": round((matches["walk_top60_share"] - matches["smeared_top60_share"]).mean(), 4),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    daily = load_minute_steps().merge(classify_walk_days(), on="date", how="left")
    matches = build_matches(daily)
    summary = summarize(daily, matches)

    daily.to_csv(DAILY_PATH, index=False)
    matches.to_csv(MATCHES_PATH, index=False)
    summary.to_csv(SUMMARY_PATH, index=False)

    single = summary[summary["group"] == "single_daylight_walk"].iloc[0]
    print(
        f"Minute-step coverage window: {daily['date'].min()} to {daily['date'].max()} "
        f"({len(daily)} noon-to-7pm days)."
    )
    print(
        f"Single daylight walks: {int(single['days'])} days, "
        f"mean {single['mean_window_steps']} noon-to-7pm steps, "
        f"mean span {single['mean_span_minutes']} min."
    )
    if not matches.empty:
        diff = summary[summary["group"] == "single_walk_minus_smeared_match"].iloc[0]
        print(
            f"Single walk vs smeared-step matches: active minutes {diff['mean_active_minutes']:+.2f}, "
            f"span {diff['mean_span_minutes']:+.2f}, "
            f"steps/span/min {diff['mean_steps_per_span_min']:+.4f}, "
            f"top60 share {diff['mean_top60_share']:+.4f}."
        )


if __name__ == "__main__":
    main()
