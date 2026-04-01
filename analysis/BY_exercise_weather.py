#!/usr/bin/env python3
"""Summarize Samsung's per-exercise weather snapshot table.

FINDINGS
- Samsung stores one weather snapshot per exercise with `exercise_id`,
  start-time weather, latitude, longitude, temperature, humidity, phrase,
  wind, and provider.
- This is not hourly historical weather through the workout. It is a single
  start-time snapshot.
- Coverage is uneven by year and disappears in the newest archive after 2023.
- At least one humidity value is invalid (`-1`), so this table needs light
  quality control before climate analysis.

Artifacts
- analysis/BY_exercise_weather_joined.csv
- analysis/BY_exercise_weather_summary.csv
"""

from __future__ import annotations

import csv
import statistics as stats
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SAMSUNG_ARCHIVE = Path(
    "/mnt/c/Users/Lucent/OneDrive/Documents/Backup/Samsung/samsunghealth_michael_20260324183778.7z"
)
TMP_ROOT = Path("/tmp")

EXERCISE_TYPE_MAP = {
    "1001": "walking",
    "1002": "running",
    "11007": "bike",
    "13001": "hiking",
    "15003": "indoor_bike",
    "9001": "pilates",
    "9002": "yoga",
    "0": "other",
}


def parse_ts(value: str) -> datetime:
    return datetime.strptime(value, "%Y-%m-%d %H:%M:%S.%f")


def extract_needed_files() -> Path:
    import subprocess
    import tempfile

    tmpdir = Path(tempfile.mkdtemp(prefix="by_ex_weather_", dir=str(TMP_ROOT)))
    subprocess.run(
        [
            "7z",
            "e",
            "-y",
            f"-o{tmpdir}",
            str(SAMSUNG_ARCHIVE),
            "samsunghealth_michael_20260324183778/com.samsung.shealth.exercise.20260324183778.csv",
            "samsunghealth_michael_20260324183778/com.samsung.shealth.exercise.weather.20260324183778.csv",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return tmpdir


def infer_label(exercise_type: str) -> str:
    return EXERCISE_TYPE_MAP.get(str(exercise_type), f"code_{exercise_type}")


def load_exercises(base: Path) -> dict[str, dict[str, object]]:
    path = base / "com.samsung.shealth.exercise.20260324183778.csv"
    out: dict[str, dict[str, object]] = {}
    with path.open(newline="", encoding="utf-8-sig") as fh:
        next(fh)
        reader = csv.DictReader(fh)
        for row in reader:
            raw_start = row.get("com.samsung.health.exercise.start_time")
            uuid = row.get("com.samsung.health.exercise.datauuid")
            if not raw_start or not uuid:
                continue
            start_time = parse_ts(raw_start)
            duration_min = float(row["com.samsung.health.exercise.duration"] or 0) / 60000.0
            distance_km = float(row.get("com.samsung.health.exercise.distance") or 0) / 1000.0
            exercise_type = row["com.samsung.health.exercise.exercise_type"]
            out[uuid] = {
                "exercise_id": uuid,
                "date": start_time.date().isoformat(),
                "start_time": start_time.isoformat(sep=" "),
                "year": start_time.year,
                "exercise_type_code": exercise_type,
                "exercise_type": infer_label(exercise_type),
                "duration_min": duration_min,
                "distance_km": distance_km,
            }
    return out


def load_weather(base: Path) -> list[dict[str, str]]:
    path = base / "com.samsung.shealth.exercise.weather.20260324183778.csv"
    with path.open(newline="", encoding="utf-8-sig") as fh:
        next(fh)
        reader = csv.DictReader(fh)
        return list(reader)


def fmt(value: float | None, digits: int = 2) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}f}"


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    base = extract_needed_files()
    exercises = load_exercises(base)
    weather_rows = load_weather(base)

    joined_rows: list[dict[str, object]] = []
    by_year_total = Counter(int(v["year"]) for v in exercises.values())
    by_year_match = Counter()
    by_type_total = Counter(str(v["exercise_type"]) for v in exercises.values())
    by_type_match = Counter()
    by_provider = Counter()
    temps: list[float] = []
    humidity_all: list[float] = []
    humidity_valid: list[float] = []
    invalid_humidity = 0

    for row in weather_rows:
        exercise_id = row["exercise_id"]
        exercise = exercises.get(exercise_id)
        if exercise is None:
            continue
        year = int(exercise["year"])
        label = str(exercise["exercise_type"])
        by_year_match[year] += 1
        by_type_match[label] += 1
        provider = (row.get("content_provider") or "").strip()
        by_provider[provider] += 1

        temp_c = float(row["temperature"]) if row.get("temperature") else None
        humidity = float(row["humidity"]) if row.get("humidity") else None
        if temp_c is not None:
            temps.append(temp_c)
        if humidity is not None:
            humidity_all.append(humidity)
            if 0 <= humidity <= 100:
                humidity_valid.append(humidity)
            else:
                invalid_humidity += 1

        joined_rows.append(
            {
                "date": exercise["date"],
                "start_time": exercise["start_time"],
                "year": year,
                "exercise_id": exercise_id,
                "exercise_type": label,
                "exercise_type_code": exercise["exercise_type_code"],
                "duration_min": fmt(float(exercise["duration_min"])),
                "distance_km": fmt(float(exercise["distance_km"]), 3),
                "latitude": row.get("latitude", ""),
                "longitude": row.get("longitude", ""),
                "temperature_c": fmt(temp_c),
                "humidity": fmt(humidity),
                "humidity_valid": "" if humidity is None else int(0 <= humidity <= 100),
                "phrase": row.get("phrase", ""),
                "wind_speed": row.get("wind_speed", ""),
                "wind_direction": row.get("wind_direction", ""),
                "uv_index": row.get("uv_index", ""),
                "provider": provider,
            }
        )

    joined_rows.sort(key=lambda r: (str(r["start_time"]), str(r["exercise_id"])))

    summary_rows: list[dict[str, object]] = []
    summary_rows.append(
        {
            "dimension": "overall",
            "key": "all",
            "total_exercises": len(exercises),
            "weather_rows": len(joined_rows),
            "coverage_pct": fmt(100.0 * len(joined_rows) / len(exercises)),
            "temperature_median_c": fmt(stats.median(temps) if temps else None),
            "temperature_min_c": fmt(min(temps) if temps else None),
            "temperature_max_c": fmt(max(temps) if temps else None),
            "humidity_median": fmt(stats.median(humidity_valid) if humidity_valid else None),
            "humidity_min": fmt(min(humidity_valid) if humidity_valid else None),
            "humidity_max": fmt(max(humidity_valid) if humidity_valid else None),
            "invalid_humidity_rows": invalid_humidity,
            "notes": "single start-time weather snapshot per exercise",
        }
    )

    for year in sorted(by_year_total):
        matches = by_year_match[year]
        total = by_year_total[year]
        summary_rows.append(
            {
                "dimension": "year",
                "key": str(year),
                "total_exercises": total,
                "weather_rows": matches,
                "coverage_pct": fmt(100.0 * matches / total if total else None),
                "temperature_median_c": "",
                "temperature_min_c": "",
                "temperature_max_c": "",
                "humidity_median": "",
                "humidity_min": "",
                "humidity_max": "",
                "invalid_humidity_rows": "",
                "notes": "",
            }
        )

    for label in sorted(by_type_total):
        matches = by_type_match[label]
        total = by_type_total[label]
        summary_rows.append(
            {
                "dimension": "exercise_type",
                "key": label,
                "total_exercises": total,
                "weather_rows": matches,
                "coverage_pct": fmt(100.0 * matches / total if total else None),
                "temperature_median_c": "",
                "temperature_min_c": "",
                "temperature_max_c": "",
                "humidity_median": "",
                "humidity_min": "",
                "humidity_max": "",
                "invalid_humidity_rows": "",
                "notes": "",
            }
        )

    for provider, count in sorted(by_provider.items()):
        summary_rows.append(
            {
                "dimension": "provider",
                "key": provider,
                "total_exercises": "",
                "weather_rows": count,
                "coverage_pct": "",
                "temperature_median_c": "",
                "temperature_min_c": "",
                "temperature_max_c": "",
                "humidity_median": "",
                "humidity_min": "",
                "humidity_max": "",
                "invalid_humidity_rows": "",
                "notes": "",
            }
        )

    write_csv(
        ROOT / "analysis" / "BY_exercise_weather_joined.csv",
        joined_rows,
        [
            "date",
            "start_time",
            "year",
            "exercise_id",
            "exercise_type",
            "exercise_type_code",
            "duration_min",
            "distance_km",
            "latitude",
            "longitude",
            "temperature_c",
            "humidity",
            "humidity_valid",
            "phrase",
            "wind_speed",
            "wind_direction",
            "uv_index",
            "provider",
        ],
    )
    write_csv(
        ROOT / "analysis" / "BY_exercise_weather_summary.csv",
        summary_rows,
        [
            "dimension",
            "key",
            "total_exercises",
            "weather_rows",
            "coverage_pct",
            "temperature_median_c",
            "temperature_min_c",
            "temperature_max_c",
            "humidity_median",
            "humidity_min",
            "humidity_max",
            "invalid_humidity_rows",
            "notes",
        ],
    )

    print(f"exercises={len(exercises)}")
    print(f"weather_rows={len(weather_rows)}")
    print(f"matched_weather_rows={len(joined_rows)}")
    print(f"invalid_humidity_rows={invalid_humidity}")


if __name__ == "__main__":
    main()
