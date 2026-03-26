#!/usr/bin/env python3
"""Extract daily steps, sleep, and exercise sessions from Samsung Health export.

Finds the most recent export files matching the expected naming pattern,
extracts steps, sleep, and exercise sessions, and writes CSVs.

Steps: uses step_daily_trend with source_type=-2 (deduplicated watch+phone).
       Takes the latest update per day.

Sleep: uses sleep start/end times (stored in UTC), converts to local time
       using the per-record time_offset. Assigns sleep to the date of wake-up
       (end_time in local), since a 2am-10am sleep belongs to that calendar day.

Exercises: uses the session-level exercise export with local start/end,
           duration, exercise type, step count, distance, and calories.
           Also writes a daily coverage table comparing summed exercise-session
           steps against the full daily step total.

Outputs:
  steps.csv  — date, steps, distance, speed
  sleep.csv  — date, sleep_start, sleep_end, sleep_hours, time_offset
  exercises.csv — session-level exercise bouts with local start/end times
                  plus inferred exercise labels
  exercise_daily.csv — daily exercise-session sums vs total steps
"""

import csv
import fnmatch
import subprocess
import sys
import tempfile
from datetime import datetime, timezone, timedelta
from pathlib import Path

EXPORT_DIR = Path(__file__).parent
SAMSUNG_BACKUP = Path("/mnt/c/Users/Lucent/OneDrive/Documents/Backup/Samsung")

# Files we need from the archive (glob-style, matched against member names)
NEEDED_FILES = [
    "com.samsung.shealth.step_daily_trend.*.csv",
    "com.samsung.shealth.sleep.*.csv",
    "com.samsung.health.sleep.*.csv",
    "com.samsung.shealth.exercise.*.csv",
]


EXERCISE_TYPE_MAP = {
    "1001": ("walking", "high"),
    "1002": ("running", "high"),
    "11007": ("bike", "high"),
    "15003": ("indoor_bike", "high"),
}


def find_newest_7z():
    """Find the newest .7z file in the Samsung backup directory."""
    archives = sorted(SAMSUNG_BACKUP.glob("*.7z"))
    if not archives:
        print(f"  No .7z files in {SAMSUNG_BACKUP}", file=sys.stderr)
        sys.exit(1)
    return archives[-1]


def extract_from_7z(archive, tmpdir):
    """Extract only the needed CSVs from a 7z archive into tmpdir.

    Returns the tmpdir Path (files land in a subdirectory named after the archive).
    """
    # List archive contents to find exact paths for our needed files
    result = subprocess.run(
        ["7z", "l", "-ba", str(archive)],
        capture_output=True, text=True, check=True,
    )
    to_extract = []
    for line in result.stdout.splitlines():
        # -ba format: "attr size compressed name" — name is last field
        parts = line.split()
        if len(parts) < 4:
            continue
        name = parts[-1]
        basename = name.split("/")[-1] if "/" in name else name
        for pattern in NEEDED_FILES:
            if fnmatch.fnmatch(basename, pattern):
                to_extract.append(name)
                break

    if not to_extract:
        print("  No matching files found in archive", file=sys.stderr)
        sys.exit(1)

    print(f"  Extracting from {archive.name}:", file=sys.stderr)
    for f in to_extract:
        print(f"    {f}", file=sys.stderr)

    subprocess.run(
        ["7z", "x", f"-o{tmpdir}", str(archive)] + to_extract,
        capture_output=True, text=True, check=True,
    )
    return Path(tmpdir)


def find_latest_file(search_dir, pattern):
    """Find the most recent Samsung Health export file matching pattern."""
    matches = sorted(search_dir.rglob(pattern))
    if not matches:
        print(f"  No file matching {pattern}", file=sys.stderr)
        return None
    return matches[-1]


def find_latest_exercise_file(search_dir):
    """Find the main exercise session export, not exercise.* side tables."""
    matches = sorted(search_dir.rglob("com.samsung.shealth.exercise.*.csv"))
    matches = [m for m in matches if m.name.count(".") == 5]
    if not matches:
        print("  No main exercise export file found", file=sys.stderr)
        return None
    return matches[-1]


def parse_offset(offset_str):
    """Parse 'UTC-0500' or 'UTC+0300' to a timedelta."""
    if not offset_str or not offset_str.startswith("UTC"):
        return timedelta(hours=0)
    sign = 1 if "+" in offset_str else -1
    num = offset_str.replace("UTC", "").replace("+", "").replace("-", "")
    hours = int(num[:2]) if len(num) >= 2 else 0
    minutes = int(num[2:4]) if len(num) >= 4 else 0
    return timedelta(hours=sign * hours, minutes=sign * minutes)


def parse_local_datetime(ts, offset_str):
    """Parse a Samsung timestamp and convert it to local time using offset."""
    if not ts:
        return None
    for fmt in ("%Y-%m-%d %H:%M:%S.%f", "%Y-%m-%d %H:%M:%S"):
        try:
            return datetime.strptime(ts[: len(datetime.now().strftime(fmt))], fmt) + parse_offset(offset_str)
        except ValueError:
            continue
    return None


def infer_exercise_label(exercise_type):
    """Return an inferred label and confidence for a Samsung exercise code."""
    label, confidence = EXERCISE_TYPE_MAP.get(str(exercise_type), ("unknown", "low"))
    return label, confidence


def extract_steps(search_dir):
    """Extract daily step counts from step_daily_trend CSV."""
    f = find_latest_file(search_dir, "com.samsung.shealth.step_daily_trend.*.csv")
    if not f:
        return []

    with open(f) as fh:
        next(fh)  # skip metadata row
        reader = csv.DictReader(fh)
        rows = [r for r in reader if r.get("source_type") == "-2"]

    print(f"  Steps: {len(rows)} source_type=-2 rows from {f.name}", file=sys.stderr)

    # Group by date, keep latest update per day
    by_date = {}
    for r in rows:
        try:
            ms = int(r["day_time"])
        except (ValueError, KeyError):
            continue
        date = datetime.fromtimestamp(ms / 1000, tz=timezone.utc).strftime("%Y-%m-%d")
        update = r.get("update_time", "")

        if date not in by_date or update > by_date[date]["update_time"]:
            by_date[date] = {
                "date": date,
                "steps": r.get("count", "0"),
                "distance": r.get("distance", ""),
                "speed": r.get("speed", ""),
                "calorie": r.get("calorie", ""),
                "update_time": update,
            }

    results = []
    for date in sorted(by_date):
        d = by_date[date]
        # Convert active time: speed * distance / count gives rough check,
        # but we don't have active_time in this file. Just output what we have.
        results.append({
            "date": d["date"],
            "steps": d["steps"],
            "distance": d["distance"],
            "speed": d["speed"],
        })

    return results


def extract_sleep(search_dir):
    """Extract sleep periods from sleep CSV."""
    f = find_latest_file(search_dir, "com.samsung.health.sleep.*.csv") or \
        find_latest_file(search_dir, "com.samsung.shealth.sleep.*.csv")
    if not f:
        return []

    with open(f, encoding="utf-8-sig") as fh:
        next(fh)  # skip metadata row
        reader = csv.DictReader(fh)
        rows = list(reader)

    print(f"  Sleep: {len(rows)} rows from {f.name}", file=sys.stderr)

    results = []
    for r in rows:
        start_utc = r.get("com.samsung.health.sleep.start_time", "")
        end_utc = r.get("com.samsung.health.sleep.end_time", "")
        offset_str = r.get("com.samsung.health.sleep.time_offset", "")

        if not start_utc or not end_utc:
            continue

        try:
            start_dt = datetime.strptime(start_utc[:23], "%Y-%m-%d %H:%M:%S.%f")
        except ValueError:
            try:
                start_dt = datetime.strptime(start_utc[:19], "%Y-%m-%d %H:%M:%S")
            except ValueError:
                continue

        try:
            end_dt = datetime.strptime(end_utc[:23], "%Y-%m-%d %H:%M:%S.%f")
        except ValueError:
            try:
                end_dt = datetime.strptime(end_utc[:19], "%Y-%m-%d %H:%M:%S")
            except ValueError:
                continue

        offset = parse_offset(offset_str)
        start_local = start_dt + offset
        end_local = end_dt + offset
        hours = (end_dt - start_dt).total_seconds() / 3600

        # Assign to the date of wake-up (end_time local)
        date = end_local.strftime("%Y-%m-%d")

        results.append({
            "date": date,
            "sleep_start": start_local.strftime("%H:%M"),
            "sleep_end": end_local.strftime("%H:%M"),
            "sleep_hours": f"{hours:.2f}",
            "time_offset": offset_str,
        })

    # If multiple sleep entries per day, keep the longest (primary sleep)
    by_date = {}
    for r in results:
        d = r["date"]
        if d not in by_date or float(r["sleep_hours"]) > float(by_date[d]["sleep_hours"]):
            by_date[d] = r

    return [by_date[d] for d in sorted(by_date)]


def extract_exercises(search_dir):
    """Extract exercise sessions from Samsung Health export."""
    f = find_latest_exercise_file(search_dir)
    if not f:
        return []

    with open(f, encoding="utf-8-sig") as fh:
        next(fh)  # skip metadata row
        reader = csv.DictReader(fh)
        rows = list(reader)

    print(f"  Exercises: {len(rows)} rows from {f.name}", file=sys.stderr)

    sessions = {}
    for r in rows:
        datauuid = r.get("com.samsung.health.exercise.datauuid", "")
        start_local = parse_local_datetime(
            r.get("com.samsung.health.exercise.start_time", ""),
            r.get("com.samsung.health.exercise.time_offset", ""),
        )
        end_local = parse_local_datetime(
            r.get("com.samsung.health.exercise.end_time", ""),
            r.get("com.samsung.health.exercise.time_offset", ""),
        )
        if not datauuid or start_local is None or end_local is None:
            continue

        update_time = r.get("com.samsung.health.exercise.update_time", "")
        row = {
            "date": start_local.strftime("%Y-%m-%d"),
            "start_time": start_local.strftime("%Y-%m-%d %H:%M:%S"),
            "end_time": end_local.strftime("%Y-%m-%d %H:%M:%S"),
            "duration_min": round(
                float(r.get("com.samsung.health.exercise.duration", 0) or 0) / 60000, 2
            ),
            "exercise_type": r.get("com.samsung.health.exercise.exercise_type", ""),
            "title": r.get("title", ""),
            "source_type": r.get("source_type", ""),
            "count": r.get("com.samsung.health.exercise.count", ""),
            "distance": r.get("com.samsung.health.exercise.distance", ""),
            "calorie": r.get("com.samsung.health.exercise.calorie", ""),
            "time_offset": r.get("com.samsung.health.exercise.time_offset", ""),
            "pkg_name": r.get("com.samsung.health.exercise.pkg_name", ""),
            "datauuid": datauuid,
            "_update_time": update_time,
        }
        row["exercise_label"], row["exercise_label_confidence"] = infer_exercise_label(
            row["exercise_type"]
        )

        if datauuid not in sessions or update_time > sessions[datauuid]["_update_time"]:
            sessions[datauuid] = row

    results = []
    for datauuid in sorted(sessions, key=lambda k: sessions[k]["start_time"]):
        row = sessions[datauuid].copy()
        row.pop("_update_time", None)
        results.append(row)
    return results


def summarize_exercises_daily(exercises, steps):
    """Aggregate exercise sessions by day and compare with total step counts."""
    step_map = {row["date"]: int(float(row["steps"])) for row in steps}
    by_date = {}
    for row in exercises:
        d = row["date"]
        by_date.setdefault(
            d,
            {
                "date": d,
                "exercise_sessions": 0,
                "exercise_duration_min": 0.0,
                "exercise_steps": 0.0,
                "exercise_distance": 0.0,
                "exercise_calorie": 0.0,
            },
        )
        by_date[d]["exercise_sessions"] += 1
        by_date[d]["exercise_duration_min"] += float(row["duration_min"] or 0)
        by_date[d]["exercise_steps"] += float(row["count"] or 0)
        by_date[d]["exercise_distance"] += float(row["distance"] or 0)
        by_date[d]["exercise_calorie"] += float(row["calorie"] or 0)

    results = []
    for d in sorted(set(step_map) | set(by_date)):
        summary = by_date.get(
            d,
            {
                "date": d,
                "exercise_sessions": 0,
                "exercise_duration_min": 0.0,
                "exercise_steps": 0.0,
                "exercise_distance": 0.0,
                "exercise_calorie": 0.0,
            },
        )
        total_steps = step_map.get(d, 0)
        summary["total_steps"] = total_steps
        summary["exercise_step_fraction"] = (
            round(summary["exercise_steps"] / total_steps, 4) if total_steps else ""
        )
        summary["exercise_duration_min"] = round(summary["exercise_duration_min"], 2)
        summary["exercise_steps"] = round(summary["exercise_steps"], 2)
        summary["exercise_distance"] = round(summary["exercise_distance"], 2)
        summary["exercise_calorie"] = round(summary["exercise_calorie"], 2)
        results.append(summary)
    return results


def main():
    archive = find_newest_7z()
    print(f"Extracting Samsung Health data from {archive.name}...", file=sys.stderr)

    with tempfile.TemporaryDirectory() as tmpdir:
        extract_from_7z(archive, tmpdir)
        steps = extract_steps(Path(tmpdir))
        sleep = extract_sleep(Path(tmpdir))
        exercises = extract_exercises(Path(tmpdir))

    exercise_daily = summarize_exercises_daily(exercises, steps)

    # Write steps
    steps_path = EXPORT_DIR / "steps.csv"
    with open(steps_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["date", "steps", "distance", "speed"])
        writer.writeheader()
        writer.writerows(steps)

    # Write sleep
    sleep_path = EXPORT_DIR / "sleep.csv"
    with open(sleep_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["date", "sleep_start", "sleep_end",
                                                "sleep_hours", "time_offset"])
        writer.writeheader()
        writer.writerows(sleep)

    # Write exercise sessions
    exercise_path = EXPORT_DIR / "exercises.csv"
    with open(exercise_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "date",
                "start_time",
                "end_time",
                "duration_min",
                "exercise_type",
                "exercise_label",
                "exercise_label_confidence",
                "title",
                "source_type",
                "count",
                "distance",
                "calorie",
                "time_offset",
                "pkg_name",
                "datauuid",
            ],
        )
        writer.writeheader()
        writer.writerows(exercises)

    # Write daily exercise coverage summary
    exercise_daily_path = EXPORT_DIR / "exercise_daily.csv"
    with open(exercise_daily_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "date",
                "exercise_sessions",
                "exercise_duration_min",
                "exercise_steps",
                "exercise_distance",
                "exercise_calorie",
                "total_steps",
                "exercise_step_fraction",
            ],
        )
        writer.writeheader()
        writer.writerows(exercise_daily)

    print(f"\nSteps: {len(steps)} days ({steps[0]['date']} to {steps[-1]['date']})",
          file=sys.stderr)
    print(f"Sleep: {len(sleep)} days ({sleep[0]['date']} to {sleep[-1]['date']})",
          file=sys.stderr)
    if exercises:
        fractions = [row["exercise_step_fraction"] for row in exercise_daily if row["exercise_step_fraction"] != ""]
        mean_fraction = sum(fractions) / len(fractions) if fractions else 0
        print(
            f"Exercises: {len(exercises)} sessions ({exercises[0]['date']} to {exercises[-1]['date']})",
            file=sys.stderr,
        )
        print(f"Exercise-step coverage mean: {mean_fraction:.3f}", file=sys.stderr)
    print(f"Written: {steps_path}", file=sys.stderr)
    print(f"Written: {sleep_path}", file=sys.stderr)
    print(f"Written: {exercise_path}", file=sys.stderr)
    print(f"Written: {exercise_daily_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
