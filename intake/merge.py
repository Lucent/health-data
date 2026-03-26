#!/usr/bin/env python3
"""Merge all intake data into unified CSVs.

Priority waterfall for overlapping dates: HTML > MHTML > OXPS
(HTML has the most complete data; OXPS is fallback for 2013-2017).

Produces:
  intake/intake_foods.csv        — every food item with all 8 nutrients
  intake/intake_daily.csv        — daily totals (from TOTAL rows)
  steps-sleep/mfp_exercises.csv  — MFP exercise entries (date, name, calories, minutes)
"""

import csv
import subprocess
import sys
import io
from collections import defaultdict
from pathlib import Path

INTAKE_DIR = Path(__file__).parent
REPO_ROOT = INTAKE_DIR.parent

NUTRIENT_COLS = ["calories", "carbs_g", "fat_g", "protein_g",
                 "cholest_mg", "sodium_mg", "sugars_g", "fiber_g"]

EXERCISE_COLS = ["date", "name", "calories", "minutes"]

FOOD_HEADER = ["date", "meal", "food"] + NUTRIENT_COLS + ["source"]
DAILY_HEADER = ["date"] + NUTRIENT_COLS
EXERCISE_HEADER = EXERCISE_COLS + ["source"]


def run_extractor(script, files):
    """Run an extractor and return (food_rows, exercise_rows).

    Extractors output food CSV, then '---EXERCISES---' sentinel, then exercise CSV.
    """
    if not files:
        return [], []
    args = [sys.executable, str(script)] + [str(f) for f in files]
    r = subprocess.run(args, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  ERROR running {script.name}: {r.stderr[:200]}", file=sys.stderr)
        return [], []

    parts = r.stdout.split("---EXERCISES---\n", 1)
    food_rows = list(csv.DictReader(io.StringIO(parts[0])))
    exercise_rows = []
    if len(parts) > 1:
        exercise_rows = list(csv.DictReader(io.StringIO(parts[1])))
    return food_rows, exercise_rows


def collect_all_files():
    """Collect all source files organized by year and type."""
    years = {}
    for year_dir in sorted(INTAKE_DIR.iterdir()):
        if not year_dir.is_dir() or year_dir.name in ("quarter-depends",):
            continue
        try:
            year = int(year_dir.name)
        except ValueError:
            continue

        html_files = sorted(year_dir.glob("*.html"))
        monthly_mhtml = sorted(f for f in year_dir.glob("*.mhtml")
                               if "new layout" not in f.name and "Q" not in f.name)
        quarterly_mhtml = sorted(f for f in year_dir.glob("*.mhtml")
                                 if "Q" in f.name and "new layout" not in f.name)
        # Only use quarterly MHTMLs if no monthly MHTMLs exist for that quarter
        mhtml_files = list(monthly_mhtml)
        if quarterly_mhtml and not monthly_mhtml:
            mhtml_files.extend(quarterly_mhtml)
        oxps_files = sorted(f for f in year_dir.glob("*.oxps")
                            if not f.stem.endswith(str(year)))  # skip annual summaries

        years[year] = {
            "html": html_files,
            "mhtml": mhtml_files,
            "oxps": oxps_files,
        }
    return years


def main():
    oxps_script = INTAKE_DIR / "extract_oxps.py"
    html_script = INTAKE_DIR / "extract_html.py"

    years = collect_all_files()
    print(f"Found {len(years)} years: {min(years)}-{max(years)}", file=sys.stderr)

    # Collect all data, tagged by source
    all_foods = []  # (date, meal, food, nutrients..., source)
    all_totals = {}  # date -> {nutrient cols}
    all_exercises = []  # (date, name, calories, minutes, source)

    # Track which dates are covered by which source (for priority)
    dates_by_source = {"html": set(), "mhtml": set(), "oxps": set()}
    exercise_dates_by_source = {"html": set(), "mhtml": set(), "oxps": set()}

    for year in sorted(years):
        files = years[year]
        print(f"  {year}...", end="", file=sys.stderr, flush=True)

        # Extract from each source
        for source, script, file_list in [
            ("html", html_script, files["html"]),
            ("mhtml", html_script, files["mhtml"]),
            ("oxps", oxps_script, files["oxps"]),
        ]:
            if not file_list:
                continue

            food_rows, exercise_rows = run_extractor(script, file_list)
            for r in food_rows:
                date = r["date"]
                dates_by_source[source].add(date)

                if r["meal"] == "TOTAL":
                    # Store total, will pick best source later
                    key = (date, source)
                    all_totals.setdefault(date, {})[source] = r
                else:
                    all_foods.append({
                        "date": date,
                        "meal": r["meal"],
                        "food": r["food"],
                        **{col: r.get(col, "") for col in NUTRIENT_COLS},
                        "source": source,
                    })

            for r in exercise_rows:
                date = r["date"]
                exercise_dates_by_source[source].add(date)
                all_exercises.append({
                    "date": date,
                    "name": r["name"],
                    "calories": r["calories"],
                    "minutes": r["minutes"],
                    "source": source,
                })

        print(" done", file=sys.stderr)

    # Apply priority waterfall: for each date, keep only the highest-priority source
    # Priority: mhtml > html > oxps
    # MHTML takes priority over HTML because MHTMLs are more recent exports
    # with corrected food entries. Quarterly HTMLs (2011-2012, 2019 Q1-Q3)
    # are older snapshots. But 2011-2012 HTMLs are the ONLY HTML source for
    # those dates (no MHTML exists), so HTML is still used when MHTML is absent.
    date_source = {}
    for date in sorted(set(r["date"] for r in all_foods)):
        if date in dates_by_source["mhtml"]:
            date_source[date] = "mhtml"
        elif date in dates_by_source["html"]:
            date_source[date] = "html"
        elif date in dates_by_source["oxps"]:
            date_source[date] = "oxps"

    # Filter foods to best source per date
    filtered_foods = [r for r in all_foods if r["source"] == date_source.get(r["date"])]
    filtered_foods.sort(key=lambda r: r["date"])

    # Build daily totals from best source
    daily_rows = []
    for date in sorted(date_source):
        source = date_source[date]
        if date in all_totals and source in all_totals[date]:
            t = all_totals[date][source]
            daily_rows.append({
                "date": date,
                **{col: t.get(col, "") for col in NUTRIENT_COLS},
            })

    # Write foods CSV
    foods_path = INTAKE_DIR / "intake_foods.csv"
    with open(foods_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FOOD_HEADER)
        writer.writeheader()
        writer.writerows(filtered_foods)

    # Write daily CSV
    daily_path = INTAKE_DIR / "intake_daily.csv"
    with open(daily_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=DAILY_HEADER)
        writer.writeheader()
        writer.writerows(daily_rows)

    # Apply same priority waterfall to exercises
    exercise_date_source = {}
    exercise_dates_all = set(r["date"] for r in all_exercises)
    for date in sorted(exercise_dates_all):
        if date in exercise_dates_by_source["mhtml"]:
            exercise_date_source[date] = "mhtml"
        elif date in exercise_dates_by_source["html"]:
            exercise_date_source[date] = "html"
        elif date in exercise_dates_by_source["oxps"]:
            exercise_date_source[date] = "oxps"

    filtered_exercises = [r for r in all_exercises
                          if r["source"] == exercise_date_source.get(r["date"])]
    filtered_exercises.sort(key=lambda r: r["date"])

    # Write exercises CSV
    exercises_path = REPO_ROOT / "steps-sleep" / "mfp_exercises.csv"
    with open(exercises_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=EXERCISE_HEADER)
        writer.writeheader()
        writer.writerows(filtered_exercises)

    # Summary
    all_dates = sorted(date_source.keys())
    source_counts = defaultdict(int)
    for d, s in date_source.items():
        source_counts[s] += 1

    ex_source_counts = defaultdict(int)
    for d, s in exercise_date_source.items():
        ex_source_counts[s] += 1

    print(f"\n=== MERGE COMPLETE ===", file=sys.stderr)
    print(f"Date range: {all_dates[0]} to {all_dates[-1]}", file=sys.stderr)
    print(f"Total dates: {len(all_dates)}", file=sys.stderr)
    print(f"  from html:  {source_counts['html']}", file=sys.stderr)
    print(f"  from mhtml: {source_counts['mhtml']}", file=sys.stderr)
    print(f"  from oxps:  {source_counts['oxps']}", file=sys.stderr)
    print(f"Food items: {len(filtered_foods)}", file=sys.stderr)
    print(f"Daily totals: {len(daily_rows)}", file=sys.stderr)
    print(f"Exercise entries: {len(filtered_exercises)} ({len(exercise_date_source)} days)", file=sys.stderr)
    print(f"  from html:  {ex_source_counts['html']}", file=sys.stderr)
    print(f"  from mhtml: {ex_source_counts['mhtml']}", file=sys.stderr)
    print(f"  from oxps:  {ex_source_counts['oxps']}", file=sys.stderr)
    print(f"Written: {foods_path}", file=sys.stderr)
    print(f"Written: {daily_path}", file=sys.stderr)
    print(f"Written: {exercises_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
