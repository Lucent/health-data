#!/usr/bin/env python3
"""Merge cardio exercise sessions into canonical exercises.csv.

Sources:
  exercises_samsung.csv — Samsung Health auto-detected sessions (2014-04+)
  exercises_mfp.csv     — pre-Samsung running entries only (dates < Samsung start)

Samsung provides: date, type, duration, distance, calorie, start/end times.
MFP pre-Samsung provides: date, name (→ type), duration only.

Produces:
  exercises.csv — canonical merged cardio exercises
"""

import csv
import re
import sys
from pathlib import Path

DIR = Path(__file__).parent

MFP_TYPE_MAP = {
    "running": "running",
    "jogging": "running",
    "walking": "walking",
}


def classify_mfp(name):
    """Map MFP exercise name to a canonical type, or None to skip."""
    lower = name.lower()
    for keyword, etype in MFP_TYPE_MAP.items():
        if keyword in lower:
            return etype
    return None


def main():
    # Samsung exercises
    samsung = []
    with open(DIR / "exercises_samsung.csv") as f:
        for r in csv.DictReader(f):
            samsung.append(r)

    samsung_start = samsung[0]["date"]

    # MFP exercises — only pre-Samsung cardio
    mfp_rows = []
    with open(DIR / "mfp_exercises.csv") as f:
        for r in csv.DictReader(f):
            if r["date"] >= samsung_start:
                continue
            etype = classify_mfp(r["name"])
            if etype is None:
                continue
            mfp_rows.append({
                "date": r["date"],
                "type": etype,
                "duration_min": r["minutes"],
                "distance": "",
                "calorie": "",
                "source": "mfp",
            })

    # Normalize Samsung rows
    canonical = []
    for r in mfp_rows:
        canonical.append(r)

    for r in samsung:
        canonical.append({
            "date": r["date"],
            "type": r["exercise_label"],
            "duration_min": r["duration_min"],
            "distance": r["distance"],
            "calorie": r["calorie"],
            "source": "samsung",
        })

    canonical.sort(key=lambda r: (r["date"], r["source"]))

    out_path = DIR / "exercises.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "date", "type", "duration_min", "distance", "calorie", "source",
        ])
        writer.writeheader()
        writer.writerows(canonical)

    print(f"Samsung: {len(samsung)} sessions", file=sys.stderr)
    print(f"MFP pre-Samsung: {len(mfp_rows)} sessions", file=sys.stderr)
    print(f"Canonical exercises: {len(canonical)} ({canonical[0]['date']} to {canonical[-1]['date']})", file=sys.stderr)
    print(f"Written: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
