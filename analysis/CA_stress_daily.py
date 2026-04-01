#!/usr/bin/env python3
"""Test whether Samsung stress scores track intake, binges, or sleep timing.

FINDINGS
- Samsung stress coverage is sparse: 962 session rows collapse to only 208
  daily stress days, almost all in 2021-2022.
- The strongest signal is same-day binge risk, but it is small and based on
  only 16 binge days total.
- Stress is more correlated with later bedtime than with sleep duration.
- This is suggestive, not a strong enough series to anchor major conclusions.

Artifacts
- analysis/CA_stress_daily.csv
- analysis/CA_stress_summary.csv
"""

from __future__ import annotations

import csv
import statistics as stats
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
SAMSUNG_ARCHIVE = Path(
    "/mnt/c/Users/Lucent/OneDrive/Documents/Backup/Samsung/samsunghealth_michael_20260324183778.7z"
)
TMP_ROOT = Path("/tmp")
BINGE_THRESHOLD = 2800


def extract_stress_csv() -> Path:
    import subprocess
    import tempfile

    tmpdir = Path(tempfile.mkdtemp(prefix="ca_stress_", dir=str(TMP_ROOT)))
    subprocess.run(
        [
            "7z",
            "e",
            "-y",
            f"-o{tmpdir}",
            str(SAMSUNG_ARCHIVE),
            "samsunghealth_michael_20260324183778/com.samsung.shealth.stress.20260324183778.csv",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return tmpdir / "com.samsung.shealth.stress.20260324183778.csv"


def bedtime_numeric(value: str | float) -> float:
    if pd.isna(value):
        return np.nan
    hour, minute = map(int, str(value).split(":"))
    out = hour + minute / 60.0
    if out < 12:
        out += 24.0
    return out


def partial_r(frame: pd.DataFrame, x: str, y: str, controls: list[str]) -> tuple[int, float]:
    cols = [x, y] + controls
    d = frame[cols].dropna()
    if len(d) < 10:
        return len(d), np.nan
    X = d[controls].to_numpy(dtype=float)
    X = np.column_stack([np.ones(len(X)), X])
    xv = d[x].to_numpy(dtype=float)
    yv = d[y].to_numpy(dtype=float)
    bx = np.linalg.lstsq(X, xv, rcond=None)[0]
    by = np.linalg.lstsq(X, yv, rcond=None)[0]
    rx = xv - X @ bx
    ry = yv - X @ by
    return len(d), float(np.corrcoef(rx, ry)[0, 1])


def main() -> None:
    stress_path = extract_stress_csv()
    with stress_path.open(newline="", encoding="utf-8-sig") as fh:
        next(fh)
        reader = csv.DictReader(fh)
        rows = []
        for row in reader:
            start_time = row.get("start_time")
            score = row.get("score")
            if not start_time or not score:
                continue
            rows.append(
                {
                    "date": start_time[:10],
                    "score": float(score),
                    "stress_min": float(row["min"]) if row.get("min") else np.nan,
                    "stress_max_raw": float(row["max"]) if row.get("max") else np.nan,
                    "tag_id": row.get("tag_id", ""),
                }
            )

    daily = (
        pd.DataFrame(rows)
        .groupby("date")
        .agg(
            stress_n=("score", "size"),
            stress_mean=("score", "mean"),
            stress_max=("score", "max"),
            stress_10011_n=("tag_id", lambda s: int((s == "10011").sum())),
            stress_10011_mean=("score", lambda s: float(s.mean())),
        )
        .reset_index()
    )
    daily["date"] = pd.to_datetime(daily["date"])

    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    sleep = pd.read_csv(ROOT / "steps-sleep" / "sleep.csv", parse_dates=["date"])
    steps = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])

    sleep["bedtime"] = sleep["sleep_start"].map(bedtime_numeric)

    daily = daily.merge(intake[["date", "calories"]], on="date", how="left")
    daily = daily.merge(sleep[["date", "sleep_hours", "bedtime"]], on="date", how="left")
    daily = daily.merge(steps[["date", "steps"]], on="date", how="left")
    daily = daily.sort_values("date").reset_index(drop=True)
    daily["year"] = daily["date"].dt.year.astype(float)
    daily["weekend"] = (daily["date"].dt.dayofweek >= 5).astype(float)
    daily["binge"] = (daily["calories"] > BINGE_THRESHOLD).astype(float)
    daily["cal_next"] = daily["calories"].shift(-1)
    daily["binge_next"] = daily["binge"].shift(-1)

    median_stress = float(daily["stress_mean"].median())
    high = daily[daily["stress_mean"] >= median_stress]
    low = daily[daily["stress_mean"] < median_stress]

    summary_rows = [
        {
            "metric": "coverage_days",
            "n": len(daily),
            "value": f"{len(daily):.0f}",
            "notes": f"{daily['date'].min().date()} to {daily['date'].max().date()}",
        },
        {
            "metric": "stress_rows",
            "n": len(rows),
            "value": f"{len(rows):.0f}",
            "notes": "",
        },
        {
            "metric": "median_measures_per_day",
            "n": len(daily),
            "value": f"{stats.median(daily['stress_n']):.1f}",
            "notes": "",
        },
        {
            "metric": "same_day_calories_r",
            "n": partial_r(daily, "stress_mean", "calories", [])[0],
            "value": f"{partial_r(daily, 'stress_mean', 'calories', [])[1]:.3f}",
            "notes": "raw Pearson r",
        },
        {
            "metric": "same_day_calories_partial_r",
            "n": partial_r(daily, "stress_mean", "calories", ["steps", "weekend", "year"])[0],
            "value": f"{partial_r(daily, 'stress_mean', 'calories', ['steps', 'weekend', 'year'])[1]:.3f}",
            "notes": "controls: steps, weekend, year",
        },
        {
            "metric": "same_day_binge_partial_r",
            "n": partial_r(daily, "stress_mean", "binge", ["steps", "weekend", "year"])[0],
            "value": f"{partial_r(daily, 'stress_mean', 'binge', ['steps', 'weekend', 'year'])[1]:.3f}",
            "notes": "controls: steps, weekend, year",
        },
        {
            "metric": "next_day_calories_partial_r",
            "n": partial_r(daily, "stress_mean", "cal_next", ["steps", "weekend", "year"])[0],
            "value": f"{partial_r(daily, 'stress_mean', 'cal_next', ['steps', 'weekend', 'year'])[1]:.3f}",
            "notes": "controls: steps, weekend, year",
        },
        {
            "metric": "later_bedtime_partial_r",
            "n": partial_r(daily, "stress_mean", "bedtime", ["steps", "weekend", "year"])[0],
            "value": f"{partial_r(daily, 'stress_mean', 'bedtime', ['steps', 'weekend', 'year'])[1]:.3f}",
            "notes": "controls: steps, weekend, year",
        },
        {
            "metric": "sleep_hours_partial_r",
            "n": partial_r(daily, "stress_mean", "sleep_hours", ["steps", "weekend", "year"])[0],
            "value": f"{partial_r(daily, 'stress_mean', 'sleep_hours', ['steps', 'weekend', 'year'])[1]:.3f}",
            "notes": "controls: steps, weekend, year",
        },
        {
            "metric": "high_vs_low_binge_rate_pct",
            "n": len(daily),
            "value": f"{100 * high['binge'].mean():.1f} vs {100 * low['binge'].mean():.1f}",
            "notes": f"median stress split at {median_stress:.1f}",
        },
        {
            "metric": "high_vs_low_calories",
            "n": len(daily),
            "value": f"{high['calories'].mean():.1f} vs {low['calories'].mean():.1f}",
            "notes": f"delta {high['calories'].mean() - low['calories'].mean():+.1f} cal",
        },
    ]

    daily_out = daily[
        [
            "date",
            "year",
            "stress_n",
            "stress_mean",
            "stress_max",
            "stress_10011_n",
            "stress_10011_mean",
            "calories",
            "binge",
            "cal_next",
            "binge_next",
            "sleep_hours",
            "bedtime",
            "steps",
        ]
    ].copy()
    daily_out["date"] = daily_out["date"].dt.strftime("%Y-%m-%d")

    daily_out.to_csv(ROOT / "analysis" / "CA_stress_daily.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(ROOT / "analysis" / "CA_stress_summary.csv", index=False)

    print(f"stress rows={len(rows)}")
    print(f"daily days={len(daily)}")
    print(f"median stress split={median_stress:.1f}")
    print(f"high-stress binge rate={100 * high['binge'].mean():.1f}%")
    print(f"low-stress binge rate={100 * low['binge'].mean():.1f}%")


if __name__ == "__main__":
    main()
