#!/usr/bin/env python3
"""Merge step counts from all sources into canonical steps.csv.

Sources (Samsung takes priority on date overlap):
  steps_samsung.csv   — Samsung Health daily steps (2014-04+)
  exercises_mfp.csv   — MFP exercises converted to step estimates
  steps_calendar.csv  — hospital shifts from calendar not in MFP

MFP → steps rules:
  Hospital              → 5500 steps (4-hour shift walking hospital floors)
  Vision, P3, L5        → treadmill at ~3.0 mph, steps from duration
  Walking, * mph, *     → steps from mph × duration (2000 steps/mile)
  Running, * mph, *     → steps from mph × duration (1600 steps/mile)
  Bicycling, *          → no steps

Produces:
  steps.csv — canonical merged (Samsung + backfill)
"""

import csv
import re
import sys
from pathlib import Path

DIR = Path(__file__).parent

HOSPITAL_STEPS = 5500
WALKING_STEPS_PER_MILE = 2000
RUNNING_STEPS_PER_MILE = 1600
TREADMILL_MPH = 3.0

MPH_RE = re.compile(r"(\d+\.?\d*)\s*mph")


def estimate_steps(name, minutes):
    """Estimate step count from exercise name and duration."""
    minutes = int(minutes)

    if name == "Hospital":
        return HOSPITAL_STEPS

    if name.startswith("Vision"):
        miles = TREADMILL_MPH * (minutes / 60)
        return round(miles * WALKING_STEPS_PER_MILE)

    mph_match = MPH_RE.search(name)
    if not mph_match:
        return None
    mph = float(mph_match.group(1))
    miles = mph * (minutes / 60)

    if "running" in name.lower() or "jogging" in name.lower():
        return round(miles * RUNNING_STEPS_PER_MILE)
    if "walking" in name.lower():
        return round(miles * WALKING_STEPS_PER_MILE)

    return None


def main():
    # Samsung steps
    samsung_rows = []
    samsung_dates = set()
    with open(DIR / "steps_samsung.csv") as f:
        for r in csv.DictReader(f):
            samsung_rows.append(r)
            samsung_dates.add(r["date"])

    # MFP exercises → step estimates (skip dates Samsung covers)
    with open(DIR / "mfp_exercises.csv") as f:
        mfp = list(csv.DictReader(f))

    backfill = {}  # date → steps
    for r in mfp:
        d = r["date"]
        if d in samsung_dates:
            continue
        steps = estimate_steps(r["name"], r["minutes"])
        if steps is not None:
            backfill[d] = backfill.get(d, 0) + steps

    # Calendar hospital dates (skip dates already covered)
    with open(DIR / "steps_calendar.csv") as f:
        for r in csv.DictReader(f):
            d = r["date"]
            if d not in samsung_dates and d not in backfill:
                backfill[d] = int(r["steps"])

    # Merge: Samsung + backfill
    merged = []
    for r in samsung_rows:
        merged.append({"date": r["date"], "steps": r["steps"],
                        "distance": r["distance"], "speed": r["speed"]})
    for d in sorted(backfill):
        merged.append({"date": d, "steps": backfill[d], "distance": "", "speed": ""})
    merged.sort(key=lambda r: r["date"])

    steps_path = DIR / "steps.csv"
    with open(steps_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["date", "steps", "distance", "speed"])
        writer.writeheader()
        writer.writerows(merged)

    print(f"Samsung: {len(samsung_rows)} days", file=sys.stderr)
    print(f"Backfill: {len(backfill)} days (MFP + calendar)", file=sys.stderr)
    print(f"Canonical steps: {len(merged)} days ({merged[0]['date']} to {merged[-1]['date']})", file=sys.stderr)
    print(f"Written: {steps_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
