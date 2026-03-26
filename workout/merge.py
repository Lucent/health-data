#!/usr/bin/env python3
"""Merge strength training dates from all sources into canonical strength.csv.

Sources (priority on date overlap: PDF > Chloe > MFP):
  workout/*.pdf            — ActivTrax YMCA sessions (2018-01 to 2025-10)
  workout/Chloe Workout.xlsx — personal trainer sessions (2016-09 to 2017-01)
  exercises_mfp.csv        — "Circuit training, general" entries

All sessions normalized to 30 minutes.

MFP Circuit dates are excluded if:
  - The date already has a PDF or Chloe session
  - The next day has a PDF (off-by-one from after-midnight MFP logging)

Produces:
  workout/strength.csv — date, duration_min, source
"""

import csv
import os
import sys
from datetime import datetime, timedelta
from pathlib import Path

import openpyxl

WORKOUT_DIR = Path(__file__).parent
REPO_ROOT = WORKOUT_DIR.parent

DURATION_MIN = 30


def load_pdf_dates():
    """Extract dates from PDF filenames (YYYY-MM-DD.pdf)."""
    dates = set()
    for f in WORKOUT_DIR.iterdir():
        if f.suffix == ".pdf" and f.stem[:4].isdigit() and len(f.stem) == 10:
            dates.add(f.stem)
    return dates


def load_chloe_dates():
    """Extract dates from Chloe Workout.xlsx."""
    wb = openpyxl.load_workbook(WORKOUT_DIR / "Chloe Workout.xlsx")
    ws = wb.active
    dates = set()
    for row in ws.iter_rows(values_only=True):
        if row[0] and hasattr(row[0], "strftime"):
            dates.add(row[0].strftime("%Y-%m-%d"))
    return dates


def load_mfp_circuit_dates():
    """Extract Circuit training dates from MFP exercises."""
    dates = set()
    with open(REPO_ROOT / "steps-sleep" / "mfp_exercises.csv") as f:
        for r in csv.DictReader(f):
            if "circuit" in r["name"].lower():
                dates.add(r["date"])
    return dates


def main():
    pdf_dates = load_pdf_dates()
    chloe_dates = load_chloe_dates()
    mfp_dates = load_mfp_circuit_dates()

    # Filter MFP: remove dates covered by PDF or Chloe, and off-by-one dupes
    mfp_filtered = set()
    for d in mfp_dates:
        if d in pdf_dates or d in chloe_dates:
            continue
        # Off-by-one: if next day has a PDF, this MFP entry is the same session
        next_day = (datetime.strptime(d, "%Y-%m-%d") + timedelta(days=1)).strftime("%Y-%m-%d")
        if next_day in pdf_dates:
            continue
        mfp_filtered.add(d)

    # Build canonical list with source priority
    rows = []
    all_dates = pdf_dates | chloe_dates | mfp_filtered
    for d in sorted(all_dates):
        if d in pdf_dates:
            source = "pdf"
        elif d in chloe_dates:
            source = "chloe"
        else:
            source = "mfp"
        rows.append({"date": d, "duration_min": DURATION_MIN, "source": source})

    out_path = WORKOUT_DIR / "strength.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["date", "duration_min", "source"])
        writer.writeheader()
        writer.writerows(rows)

    from collections import Counter
    counts = Counter(r["source"] for r in rows)
    print(f"PDF: {counts['pdf']} dates", file=sys.stderr)
    print(f"Chloe: {counts['chloe']} dates", file=sys.stderr)
    print(f"MFP Circuit: {counts['mfp']} dates (of {len(mfp_dates)} total, {len(mfp_dates) - len(mfp_filtered)} excluded)", file=sys.stderr)
    print(f"Canonical strength: {len(rows)} dates ({rows[0]['date']} to {rows[-1]['date']})", file=sys.stderr)
    print(f"Written: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
