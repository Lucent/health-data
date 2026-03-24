#!/usr/bin/env python3
"""Compare OXPS and HTML/MHTML extraction results for overlapping months.

For each month that has both OXPS and HTML/MHTML data, runs both extractors
and compares:
  1. Daily TOTAL rows (should match exactly)
  2. Individual food items (OXPS may be a subset due to print truncation)

Outputs a summary report to stdout.
"""

import csv
import io
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

INTAKE_DIR = Path(__file__).parent


def run_extractor(script, filepath):
    """Run an extractor script and return parsed rows."""
    result = subprocess.run(
        [sys.executable, str(script), str(filepath)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"  ERROR: {script.name} on {filepath}: {result.stderr}", file=sys.stderr)
        return []

    reader = csv.DictReader(io.StringIO(result.stdout))
    return list(reader)


def rows_by_date(rows):
    """Group rows by date."""
    by_date = defaultdict(list)
    for row in rows:
        by_date[row["date"]].append(row)
    return by_date


def compare_totals(oxps_rows, html_rows):
    """Compare daily TOTAL rows. Returns (matches, mismatches, missing)."""
    oxps_totals = {r["date"]: r for r in oxps_rows if r["meal"] == "TOTAL"}
    html_totals = {r["date"]: r for r in html_rows if r["meal"] == "TOTAL"}

    matches = 0
    mismatches = []
    missing_from_oxps = []
    missing_from_html = []

    all_dates = sorted(set(oxps_totals.keys()) | set(html_totals.keys()))
    for date in all_dates:
        if date not in oxps_totals:
            missing_from_oxps.append(date)
            continue
        if date not in html_totals:
            missing_from_html.append(date)
            continue

        o = oxps_totals[date]
        h = html_totals[date]

        # Compare calorie totals (the most reliable field)
        match = True
        diffs = []
        for col in ["calories", "carbs_g", "fat_g", "protein_g",
                     "cholest_mg", "sodium_mg", "sugars_g", "fiber_g"]:
            ov = o.get(col, "").strip()
            hv = h.get(col, "").strip()
            # Treat 0 and empty as equivalent for comparison
            # (OXPS renders -- as 0, HTML as empty)
            if ov == "0" and hv == "":
                continue
            if ov == hv:
                continue
            # Tolerate ±1 rounding (internal MFP rounding from exports at different times)
            try:
                if abs(int(ov) - int(hv)) <= 1:
                    continue
            except (ValueError, TypeError):
                pass
            match = False
            diffs.append(f"{col}: OXPS={ov!r} HTML={hv!r}")

        if match:
            matches += 1
        else:
            mismatches.append((date, diffs))

    return matches, mismatches, missing_from_oxps, missing_from_html


def normalize_food(name):
    """Normalize food name for comparison: collapse whitespace, strip."""
    return re.sub(r"\s+", " ", name).strip()


def compare_items(oxps_rows, html_rows):
    """Compare food item counts. Returns (oxps_count, html_count, oxps_only, html_only)."""
    oxps_items = [(r["date"], r["meal"], normalize_food(r["food"])) for r in oxps_rows if r["meal"] != "TOTAL"]
    html_items = [(r["date"], r["meal"], normalize_food(r["food"])) for r in html_rows if r["meal"] != "TOTAL"]

    oxps_set = set(oxps_items)
    html_set = set(html_items)

    return len(oxps_items), len(html_items), oxps_set - html_set, html_set - oxps_set


def find_overlap_months():
    """Find months that have both OXPS and (HTML or MHTML) files."""
    overlaps = []

    for year_dir in sorted(INTAKE_DIR.iterdir()):
        if not year_dir.is_dir() or year_dir.name == "quarter-depends":
            continue

        oxps_files = sorted(year_dir.glob("*.oxps"))
        html_files = sorted(year_dir.glob("*.html")) + sorted(year_dir.glob("*.mhtml"))
        html_files = [f for f in html_files if "new layout" not in f.name]

        if not oxps_files or not html_files:
            continue

        # For quarterly HTML files, we need to match by date range
        # For monthly files, match by month
        for oxps in oxps_files:
            # Skip annual summary files like "Food 2012.oxps"
            if not any(c.isdigit() for c in oxps.stem.split("-")[-1] if c.isdigit()):
                # Check if it's a monthly file
                pass
            overlaps.append((oxps, html_files))

    return overlaps


def main():
    oxps_script = INTAKE_DIR / "extract_oxps.py"
    html_script = INTAKE_DIR / "extract_html.py"

    total_match = 0
    total_mismatch = 0
    total_oxps_items = 0
    total_html_items = 0

    # Find all years with OXPS files
    for year_dir in sorted(INTAKE_DIR.iterdir()):
        if not year_dir.is_dir() or year_dir.name == "quarter-depends":
            continue

        oxps_files = sorted(year_dir.glob("*.oxps"))
        html_files = sorted(year_dir.glob("*.html")) + sorted(year_dir.glob("*.mhtml"))
        html_files = [f for f in html_files if "new layout" not in f.name]

        if not oxps_files or not html_files:
            continue

        year = year_dir.name
        print(f"\n{'='*60}")
        print(f"YEAR: {year}")
        print(f"{'='*60}")

        # Extract all OXPS data for this year
        oxps_rows = []
        for f in oxps_files:
            # Skip annual summaries
            if f.stem.endswith(year):
                continue
            oxps_rows.extend(run_extractor(oxps_script, f))

        # Extract all HTML/MHTML data for this year
        html_rows = []
        for f in html_files:
            html_rows.extend(run_extractor(html_script, f))

        # Filter HTML rows to only dates present in OXPS
        oxps_dates = set(r["date"] for r in oxps_rows)
        html_rows_overlap = [r for r in html_rows if r["date"] in oxps_dates]

        print(f"  OXPS dates: {len(oxps_dates)}")
        print(f"  HTML dates (overlap): {len(set(r['date'] for r in html_rows_overlap))}")

        # Compare totals
        matches, mismatches, missing_oxps, missing_html = compare_totals(oxps_rows, html_rows_overlap)
        total_match += matches
        total_mismatch += len(mismatches)

        print(f"  Daily totals: {matches} match, {len(mismatches)} mismatch")
        if mismatches:
            for date, diffs in mismatches[:5]:
                print(f"    {date}: {'; '.join(diffs)}")
            if len(mismatches) > 5:
                print(f"    ... and {len(mismatches) - 5} more")

        # Compare items
        oxps_count, html_count, oxps_only, html_only = compare_items(oxps_rows, html_rows_overlap)
        total_oxps_items += oxps_count
        total_html_items += html_count

        print(f"  Food items: OXPS={oxps_count}, HTML={html_count} (HTML has {html_count - oxps_count} more)")
        if oxps_only:
            print(f"  Items only in OXPS: {len(oxps_only)}")
            for item in sorted(oxps_only)[:3]:
                print(f"    {item}")
        if missing_oxps:
            print(f"  Dates missing from OXPS: {missing_oxps[:5]}")
        if missing_html:
            print(f"  Dates missing from HTML: {missing_html[:5]}")

    print(f"\n{'='*60}")
    print(f"OVERALL SUMMARY")
    print(f"{'='*60}")
    print(f"  Daily totals: {total_match} match, {total_mismatch} mismatch")
    print(f"  Food items: OXPS={total_oxps_items}, HTML={total_html_items}")
    print(f"  OXPS captures {total_oxps_items/max(total_html_items,1)*100:.1f}% of HTML food items")


if __name__ == "__main__":
    main()
