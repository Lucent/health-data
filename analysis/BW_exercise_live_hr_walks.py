#!/usr/bin/env python3
"""Analyze pace-adjusted walking heart rate from Samsung exercise live-data.

FINDINGS
- The per-exercise heart-rate chart lives in
  `jsons/com.samsung.shealth.exercise/*/*.com.samsung.health.exercise.live_data.json`,
  not in the generic heart-rate tracker tables.
- Using workout live-data yields a much larger and cleaner sample than the
  hourly tracker fallback for walking fitness comparisons.
- The early 2015-2017 live-HR era is on a different scale from 2021+ and
  should not be pooled without caution.

Artifacts
- analysis/BW_walk_live_hr_sessions.csv
- analysis/BW_walk_live_hr_year_summary.csv
"""

from __future__ import annotations

import csv
import json
import os
import statistics as stats
from collections import defaultdict
from datetime import date, datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SAMSUNG_ARCHIVE = Path(
    "/mnt/c/Users/Lucent/OneDrive/Documents/Backup/Samsung/samsunghealth_michael_20260324183778.7z"
)
TMP_ROOT = Path("/tmp")


def parse_ts(value: str) -> datetime:
    return datetime.strptime(value, "%Y-%m-%d %H:%M:%S.%f")


def extract_needed_files() -> Path:
    import subprocess
    import tempfile

    tmpdir = Path(tempfile.mkdtemp(prefix="bw_walk_hr_", dir=str(TMP_ROOT)))
    subprocess.run(
        [
            "7z",
            "x",
            "-y",
            f"-o{tmpdir}",
            str(SAMSUNG_ARCHIVE),
            "samsunghealth_michael_20260324183778/com.samsung.shealth.exercise.20260324183778.csv",
            "samsunghealth_michael_20260324183778/jsons/com.samsung.shealth.exercise/*/*.com.samsung.health.exercise.live_data.json",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return tmpdir / "samsunghealth_michael_20260324183778"


def load_live_map(base: Path) -> dict[str, Path]:
    live_root = base / "jsons" / "com.samsung.shealth.exercise"
    out: dict[str, Path] = {}
    for root, _, files in os.walk(live_root):
        for fn in files:
            suffix = ".com.samsung.health.exercise.live_data.json"
            if fn.endswith(suffix):
                out[fn[: -len(suffix)]] = Path(root) / fn
    return out


def load_walk_rows(base: Path) -> list[dict[str, float | int | str | datetime]]:
    live_map = load_live_map(base)
    rows: list[dict[str, float | int | str | datetime]] = []
    exercise_path = base / "com.samsung.shealth.exercise.20260324183778.csv"
    with exercise_path.open(newline="", encoding="utf-8-sig") as fh:
        next(fh)
        reader = csv.DictReader(fh)
        for row in reader:
            if row["com.samsung.health.exercise.exercise_type"] != "1001":
                continue
            uuid = row["com.samsung.health.exercise.datauuid"]
            live_path = live_map.get(uuid)
            if live_path is None:
                continue

            arr = json.loads(live_path.read_text(encoding="utf-8"))
            hrs = [
                float(obj["heart_rate"])
                for obj in arr
                if "heart_rate" in obj and 35 <= float(obj["heart_rate"]) <= 220
            ]
            if not hrs:
                continue

            dist_m = float(row.get("total_distance") or 0)
            if dist_m <= 0:
                dist_samples = [float(obj["distance"]) for obj in arr if "distance" in obj]
                if dist_samples:
                    dist_m = sum(dist_samples)
            dur_ms = float(row["com.samsung.health.exercise.duration"] or 0)
            dur_min = dur_ms / 60000.0
            if dur_min <= 0 or dist_m <= 0:
                continue

            kph = (dist_m / 1000.0) / (dur_min / 60.0)
            if not (2.5 <= kph <= 7.5) or dur_min < 15 or len(hrs) < 10:
                continue

            rows.append(
                {
                    "date": parse_ts(row["com.samsung.health.exercise.start_time"]).date(),
                    "start_time": parse_ts(row["com.samsung.health.exercise.start_time"]),
                    "duration_min": dur_min,
                    "distance_m": dist_m,
                    "kph": kph,
                    "hr_mean": stats.fmean(hrs),
                    "hr_median": stats.median(hrs),
                    "hr_max": max(hrs),
                    "hr_samples": len(hrs),
                    "datauuid": uuid,
                }
            )
    rows.sort(key=lambda r: r["start_time"])  # type: ignore[index]
    return rows


def fit_line(xs: list[float], ys: list[float]) -> tuple[float, float]:
    mx, my = stats.mean(xs), stats.mean(ys)
    varx = sum((x - mx) ** 2 for x in xs)
    if varx == 0:
        return my, 0.0
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / varx
    intercept = my - slope * mx
    return intercept, slope


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    base = extract_needed_files()
    rows = load_walk_rows(base)
    if not rows:
        raise SystemExit("No usable walking live-HR rows found.")

    xs = [float(r["kph"]) for r in rows]
    ys = [float(r["hr_mean"]) for r in rows]
    intercept, slope = fit_line(xs, ys)
    for row in rows:
        row["hr_resid_all"] = float(row["hr_mean"]) - (intercept + slope * float(row["kph"]))

    modern = [r for r in rows if int(r["date"].year) >= 2021]  # type: ignore[attr-defined]
    m_intercept, m_slope = fit_line([float(r["kph"]) for r in modern], [float(r["hr_mean"]) for r in modern])
    for row in rows:
        row["hr_resid_2021p"] = float(row["hr_mean"]) - (m_intercept + m_slope * float(row["kph"]))

    session_rows = []
    by_year: defaultdict[int, list[dict]] = defaultdict(list)
    for row in rows:
        year = int(row["date"].year)  # type: ignore[attr-defined]
        by_year[year].append(row)
        session_rows.append(
            {
                "date": row["date"].isoformat(),
                "start_time": row["start_time"].isoformat(sep=" "),
                "duration_min": f"{float(row['duration_min']):.2f}",
                "distance_m": f"{float(row['distance_m']):.1f}",
                "kph": f"{float(row['kph']):.3f}",
                "hr_mean": f"{float(row['hr_mean']):.2f}",
                "hr_median": f"{float(row['hr_median']):.2f}",
                "hr_max": f"{float(row['hr_max']):.0f}",
                "hr_samples": int(row["hr_samples"]),
                "hr_resid_all": f"{float(row['hr_resid_all']):.2f}",
                "hr_resid_2021p": f"{float(row['hr_resid_2021p']):.2f}",
                "datauuid": row["datauuid"],
            }
        )

    summary_rows = []
    for year in sorted(by_year):
        ys = by_year[year]
        summary_rows.append(
            {
                "year": year,
                "n": len(ys),
                "kph_median": f"{stats.median(float(r['kph']) for r in ys):.3f}",
                "hr_mean_median": f"{stats.median(float(r['hr_mean']) for r in ys):.2f}",
                "hr_resid_all_median": f"{stats.median(float(r['hr_resid_all']) for r in ys):.2f}",
                "hr_resid_2021p_median": f"{stats.median(float(r['hr_resid_2021p']) for r in ys):.2f}",
            }
        )

    write_csv(
        ROOT / "analysis" / "BW_walk_live_hr_sessions.csv",
        session_rows,
        [
            "date",
            "start_time",
            "duration_min",
            "distance_m",
            "kph",
            "hr_mean",
            "hr_median",
            "hr_max",
            "hr_samples",
            "hr_resid_all",
            "hr_resid_2021p",
            "datauuid",
        ],
    )
    write_csv(
        ROOT / "analysis" / "BW_walk_live_hr_year_summary.csv",
        summary_rows,
        ["year", "n", "kph_median", "hr_mean_median", "hr_resid_all_median", "hr_resid_2021p_median"],
    )

    print("=== BW. Exercise Live-Data Walk HR ===")
    print(f"Usable walking sessions: {len(rows)}")
    print(f"Date range: {rows[0]['date']} to {rows[-1]['date']}")
    print(f"All-era fit: HR = {intercept:.2f} + {slope:.2f} * kph")
    print(f"2021+ fit:   HR = {m_intercept:.2f} + {m_slope:.2f} * kph")
    print("\nYear summary:")
    for row in summary_rows:
        print(
            f"  {row['year']}: n={row['n']}, kph_med={row['kph_median']}, "
            f"hr_med={row['hr_mean_median']}, resid_2021p_med={row['hr_resid_2021p_median']}"
        )


if __name__ == "__main__":
    main()
