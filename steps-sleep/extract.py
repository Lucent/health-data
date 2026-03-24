#!/usr/bin/env python3
"""Extract daily steps and sleep from Samsung Health export.

Finds the most recent export files matching the expected naming pattern,
extracts steps and sleep, and writes daily CSVs.

Steps: uses step_daily_trend with source_type=-2 (deduplicated watch+phone).
       Takes the latest update per day.

Sleep: uses sleep start/end times (stored in UTC), converts to local time
       using the per-record time_offset. Assigns sleep to the date of wake-up
       (end_time in local), since a 2am-10am sleep belongs to that calendar day.

Outputs:
  steps.csv  — date, steps, distance, active_minutes, speed
  sleep.csv  — date, sleep_start, sleep_end, sleep_hours, time_offset
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
]


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


def parse_offset(offset_str):
    """Parse 'UTC-0500' or 'UTC+0300' to a timedelta."""
    if not offset_str or not offset_str.startswith("UTC"):
        return timedelta(hours=0)
    sign = 1 if "+" in offset_str else -1
    num = offset_str.replace("UTC", "").replace("+", "").replace("-", "")
    hours = int(num[:2]) if len(num) >= 2 else 0
    minutes = int(num[2:4]) if len(num) >= 4 else 0
    return timedelta(hours=sign * hours, minutes=sign * minutes)


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


def main():
    archive = find_newest_7z()
    print(f"Extracting Samsung Health data from {archive.name}...", file=sys.stderr)

    with tempfile.TemporaryDirectory() as tmpdir:
        extract_from_7z(archive, tmpdir)
        steps = extract_steps(Path(tmpdir))
        sleep = extract_sleep(Path(tmpdir))

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

    print(f"\nSteps: {len(steps)} days ({steps[0]['date']} to {steps[-1]['date']})",
          file=sys.stderr)
    print(f"Sleep: {len(sleep)} days ({sleep[0]['date']} to {sleep[-1]['date']})",
          file=sys.stderr)
    print(f"Written: {steps_path}", file=sys.stderr)
    print(f"Written: {sleep_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
