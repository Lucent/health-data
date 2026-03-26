#!/usr/bin/env python3
"""Extract food diary data from OXPS files (MyFitnessPal printable diary).

OXPS files are ZIP archives containing XML fixed-page documents.
Food data lives in <Glyphs UnicodeString="..."> attributes in page order.

Output:
  Food CSV to stdout:  date, meal, food, calories, carbs_g, fat_g, protein_g, cholest_mg, sodium_mg, sugars_g, fiber_g
  Exercise CSV to file (EXERCISE_FILE env var): date, name, calories, minutes
"""

import csv
import html as html_mod
import os
import re
import sys
import zipfile
from pathlib import Path

NUTRIENT_COLS = ["calories", "carbs_g", "fat_g", "protein_g",
                 "cholest_mg", "sodium_mg", "sugars_g", "fiber_g"]

MEALS = {"Breakfast", "Lunch", "Dinner", "Snacks", "Supper"}

HEADER_WORDS = {"FOODS", "Foods", "Calories", "Carbs", "Fat", "Protein",
                "Cholest", "Sodium", "Sugars", "Sugar", "Fiber",
                "Minutes", "Sets"}

EXERCISE_HEADERS = {"EXERCISES", "Exercises"}
EXERCISE_CATEGORIES = {"Cardiovascular", "Strength"}

SKIP_PREFIXES = ("Printable Diary for", "Food Diary", "Food Notes",
                 "Exercise Diary", "Exercise notes", "From:", "To:",
                 "Show:", "change report")

# Matches: "January 1, 2013", "May 15, 2017", etc.
DATE_RE = re.compile(
    r"^(January|February|March|April|May|June|July|August|September|"
    r"October|November|December)\s+\d{1,2},\s+\d{4}$"
)

# Matches a nutrient value: bare number (possibly with commas), or number + g/mg
# Also matches "--" placeholder values
NUTRIENT_RE = re.compile(r"^(?:--|-?\d[\d,]*)(?:g|mg)?$")

# ISO date at bottom of page (e.g. "2013-01-01")
ISO_DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")

MONTH_MAP = {
    "January": "01", "February": "02", "March": "03", "April": "04",
    "May": "05", "June": "06", "July": "07", "August": "08",
    "September": "09", "October": "10", "November": "11", "December": "12",
}


def parse_date(s):
    """Convert 'January 1, 2013' to '2013-01-01'."""
    m = re.match(r"(\w+)\s+(\d{1,2}),\s+(\d{4})", s)
    if not m:
        return None
    month, day, year = m.group(1), int(m.group(2)), m.group(3)
    return f"{year}-{MONTH_MAP[month]}-{int(day):02d}"


def parse_nutrient(s):
    """Normalize a nutrient string: strip units, handle '--' as empty."""
    s = s.strip()
    if s == "" or s.startswith("--"):
        return ""
    # Remove units
    s = s.replace("mg", "").replace("g", "")
    # Remove commas
    s = s.replace(",", "")
    return s


def extract_rows_from_page(zf, page_path):
    """Extract table rows from an fpage, using horizontal rules as row boundaries.

    The XPS renders horizontal lines (<Path>) between table rows. We use these
    as cell boundaries: all glyphs between two consecutive horizontal lines
    belong to the same table cell/row.

    Within each cell, glyphs are split into text (food name column, X < boundary)
    and nutrients (nutrient columns, X >= boundary).

    Returns list of rows. Each row is a dict:
      {"texts": [...], "nutrients": [...]}
    """
    with zf.open(page_path) as f:
        content = f.read().decode("utf-8")

    # --- Extract horizontal separator lines ---
    # Table row separators are rendered as thin rectangles or lines spanning
    # the full table width. They come in two styles:
    #   1. Gray (#e6e6e6) full-width lines: "F 1 M 800,Y L 0,Y 0,Y2 800,Y2 Z"
    #   2. Blue (#2b9acb) cell border rectangles: many small rects at the same Y
    # Both use the same raw coordinate space as Glyphs.
    #
    # Strategy: find all Y values where a Path element touches BOTH x≈0 and x≈800
    # on the same horizontal line (same Y within tolerance). These are row borders.
    # Collect all (x, y) from horizontal edges across ALL Path elements.
    # A horizontal edge is two consecutive points with the same Y in a path.
    # Then find Y values that span from x≈0 to x≈800 across all paths combined.
    from collections import defaultdict
    y_x_ranges = defaultdict(set)  # y -> set of x values

    for m in re.finditer(r'<Path[^>]*Data="([^"]*)"', content):
        d = m.group(1)
        coords = [(float(xm.group(1)), float(xm.group(2)))
                  for xm in re.finditer(r'(\d+\.?\d*),(\d+\.?\d*)', d)]
        # For each pair of adjacent coords with same Y (horizontal segment),
        # record the x range at that Y
        for i in range(len(coords) - 1):
            x1, y1 = coords[i]
            x2, y2 = coords[i + 1]
            if abs(y1 - y2) < 2.0:  # horizontal segment
                y_key = round(min(y1, y2))
                y_x_ranges[y_key].add(round(x1))
                y_x_ranges[y_key].add(round(x2))

    hline_ys = set()
    for y, xs in y_x_ranges.items():
        # A row border spans the full table: has points at both x≈0 and x≈800
        has_left = any(x < 5 for x in xs)
        has_right = any(x > 795 for x in xs)
        if has_left and has_right:
            hline_ys.add(y)

    hlines = sorted(hline_ys)

    # --- Extract glyphs ---
    glyphs = []
    for m in re.finditer(r"<Glyphs\s([^>]+)/>", content):
        attrs = m.group(1)
        text_m = re.search(r'UnicodeString="([^"]*)"', attrs)
        x_m = re.search(r'OriginX="([^"]*)"', attrs)
        y_m = re.search(r'OriginY="([^"]*)"', attrs)
        if text_m and x_m and y_m:
            text = html_mod.unescape(text_m.group(1))
            x = float(x_m.group(1))
            y = float(y_m.group(1))
            glyphs.append((y, x, text))

    if not glyphs:
        return []

    glyphs.sort()

    # --- Determine nutrient column X boundary ---
    nutrient_x_boundary = 400.0
    foods_y = None
    for y, x, text in glyphs:
        if text.strip() == "FOODS":
            foods_y = y
            break
    if foods_y is not None:
        for y, x, text in glyphs:
            if text.strip() == "Calories" and abs(y - foods_y) < 1.0:
                nutrient_x_boundary = x - 10
                break

    # --- Group glyphs by horizontal-line cells ---
    # For each glyph, find which pair of horizontal lines it falls between.
    # All glyphs in the same cell are one logical row.
    def cell_index(gy):
        """Return the index of the horizontal line just above this glyph."""
        idx = 0
        for i, hy in enumerate(hlines):
            if hy <= gy + 0.5:
                idx = i
            else:
                break
        return idx

    if not hlines:
        # No horizontal lines — fall back to Y-grouping
        hlines_available = False
    else:
        hlines_available = True

    rows = []

    if hlines_available:
        # Group glyphs by cell index, keeping full (y, x, text) for sorting
        from collections import OrderedDict
        cells = OrderedDict()
        for y, x, text in glyphs:
            ci = cell_index(y)
            if ci not in cells:
                cells[ci] = []
            cells[ci].append((y, x, text))

        for ci in cells:
            # Sort by Y then X so multi-line food names read top-to-bottom
            cell_glyphs = sorted(cells[ci])
            texts = [t for gy, gx, t in cell_glyphs if gx < nutrient_x_boundary]
            nutrients = [t for gy, gx, t in cell_glyphs if gx >= nutrient_x_boundary]
            rows.append({"texts": texts, "nutrients": nutrients})
    else:
        # Fallback: group by Y
        current_y = glyphs[0][0]
        current_row_glyphs = []
        for y, x, text in glyphs:
            if abs(y - current_y) > 1.0:
                if current_row_glyphs:
                    texts = [t for xv, t in sorted(current_row_glyphs) if xv < nutrient_x_boundary]
                    nutrients = [t for xv, t in sorted(current_row_glyphs) if xv >= nutrient_x_boundary]
                    rows.append({"texts": texts, "nutrients": nutrients})
                current_row_glyphs = [(x, text)]
                current_y = y
            else:
                current_row_glyphs.append((x, text))
        if current_row_glyphs:
            texts = [t for xv, t in sorted(current_row_glyphs) if xv < nutrient_x_boundary]
            nutrients = [t for xv, t in sorted(current_row_glyphs) if xv >= nutrient_x_boundary]
            rows.append({"texts": texts, "nutrients": nutrients})

    return rows


def extract_oxps(filepath):
    """Extract all food and exercise entries from an OXPS file.

    Uses Y-coordinate grouping to reconstruct table rows from the XPS layout,
    eliminating page-boundary issues entirely.

    Returns (food_results, exercise_results):
      food_results:     list of dicts with date, meal, food, + nutrient columns
      exercise_results: list of dicts with date, name, calories, minutes
    """
    food_results = []
    exercise_results = []

    with zipfile.ZipFile(filepath, "r") as zf:
        pages = sorted(
            [n for n in zf.namelist() if n.endswith(".fpage")],
            key=lambda p: int(re.search(r"/(\d+)\.fpage$", p).group(1))
        )

        all_rows = []
        for page in pages:
            all_rows.extend(extract_rows_from_page(zf, page))

    # Now we have a list of visual rows (each row = list of strings on same Y).
    # Classify each row and build food entries.
    current_date = None
    current_meal = None
    in_exercise_zone = False  # True after EXERCISES header, until next date
    dates_with_total = set()
    pending_name_parts = []  # multi-line food name accumulator
    orphaned_nutrients = None  # nutrients whose food name is on next row(s)

    for row in all_rows:
        texts = [s.strip() for s in row["texts"] if s.strip()]
        nutrients = [s.strip() for s in row["nutrients"] if s.strip()]
        first_text = texts[0] if texts else ""
        all_content = " ".join(texts + nutrients)

        # Skip empty rows
        if not all_content:
            continue

        # Date header — check FIRST, before skipping metadata,
        # because the first page's header cell may contain both
        # metadata AND the date (e.g. "April 20, 2011" mixed with
        # "2011-04-01", "Printable Diary...", etc.)
        date_candidate = None
        for s in texts + nutrients:
            if DATE_RE.match(s):
                date_candidate = s
                break
        if date_candidate:
            # Flush pending food data before switching dates
            if orphaned_nutrients is not None and pending_name_parts and not in_exercise_zone:
                food_name = " ".join(pending_name_parts)
                if food_name and current_meal:
                    entry = {"date": current_date, "meal": current_meal, "food": food_name}
                    for j, col in enumerate(NUTRIENT_COLS):
                        entry[col] = parse_nutrient(orphaned_nutrients[j])
                    food_results.append(entry)
            current_date = parse_date(date_candidate)
            current_meal = None
            in_exercise_zone = False
            pending_name_parts = []
            orphaned_nutrients = None
            continue

        # Skip footer/metadata — check all strings since some end up
        # in nutrient columns due to X positioning
        all_strings_flat = texts + nutrients
        if any(ISO_DATE_RE.match(s) for s in all_strings_flat):
            continue
        if any(s.startswith(p) for s in all_strings_flat for p in SKIP_PREFIXES):
            continue
        if any(s.startswith("\ue001") or s == "change report" for s in all_strings_flat):
            continue

        # Column header row (food headers)
        if first_text in HEADER_WORDS:
            continue

        # Exercise header — enter exercise zone
        if first_text in EXERCISE_HEADERS:
            # Flush any pending food data
            if orphaned_nutrients is not None and pending_name_parts:
                food_name = " ".join(pending_name_parts)
                if food_name and current_meal:
                    entry = {"date": current_date, "meal": current_meal, "food": food_name}
                    for j, col in enumerate(NUTRIENT_COLS):
                        entry[col] = parse_nutrient(orphaned_nutrients[j])
                    food_results.append(entry)
                orphaned_nutrients = None
            pending_name_parts = []
            in_exercise_zone = True
            continue

        # --- EXERCISE ZONE ---
        if in_exercise_zone:
            # Exercise column headers (Calories, Minutes, Sets, Reps, Weight)
            if first_text in ("Calories", "Minutes", "Sets", "Reps", "Weight"):
                continue

            # Exercise TOTALS — end of exercise section
            if any(t in ("TOTAL:", "TOTALS:") for t in texts):
                continue

            # Exercise category header (Cardiovascular, Strength)
            if first_text in EXERCISE_CATEGORIES and not nutrients:
                continue

            # Filter footer junk from nutrients
            valid_nutrients = []
            for s in nutrients:
                if not NUTRIENT_RE.match(s):
                    continue
                stripped = s.replace(",", "").replace("g", "").replace("m", "")
                if stripped.isdigit() and len(stripped) == 4 and "," not in s and not s.endswith("g"):
                    continue
                if "/" in s:
                    continue
                valid_nutrients.append(s)

            # Exercise entry: name + at least 2 values (calories, minutes)
            if texts and len(valid_nutrients) >= 2:
                exercise_results.append({
                    "date": current_date,
                    "name": " ".join(texts),
                    "calories": valid_nutrients[0],
                    "minutes": valid_nutrients[1],
                })
            continue

        # --- FOOD ZONE (unchanged logic below) ---

        # TOTAL row — check any text, not just first (TOTAL can merge with footer URL)
        if any(t in ("TOTAL:", "TOTALS:") for t in texts):
            # Flush any orphaned nutrients before processing TOTAL
            if orphaned_nutrients is not None and pending_name_parts:
                food_name = " ".join(pending_name_parts)
                if food_name and current_meal:
                    entry = {"date": current_date, "meal": current_meal, "food": food_name}
                    for j, col in enumerate(NUTRIENT_COLS):
                        entry[col] = parse_nutrient(orphaned_nutrients[j])
                    food_results.append(entry)
                orphaned_nutrients = None
            pending_name_parts = []
            # Filter footer junk from nutrients. When TOTAL merges with the
            # page footer, date fragments like "9/", "1/", "2016" land in
            # the nutrient column. Real calorie totals always have commas
            # (e.g. "1,825") or units ("282g"). Bare numbers without commas
            # or units that are 4 digits are year fragments, not nutrients.
            valid_nutrients = []
            for s in nutrients:
                if not NUTRIENT_RE.match(s):
                    continue
                # Bare number (no comma, no units) that's exactly 4 digits = year
                stripped = s.replace(",", "").replace("g", "").replace("m", "")
                if stripped.isdigit() and len(stripped) == 4 and "," not in s and not s.endswith("g"):
                    continue
                # Date fragments with slashes
                if "/" in s:
                    continue
                valid_nutrients.append(s)
            if len(valid_nutrients) == 8 and current_date not in dates_with_total:
                entry = {"date": current_date, "meal": "TOTAL", "food": "TOTAL"}
                for j, col in enumerate(NUTRIENT_COLS):
                    entry[col] = parse_nutrient(valid_nutrients[j])
                food_results.append(entry)
                dates_with_total.add(current_date)
            continue

        # Filter footer junk from texts (URLs that merge with food names)
        texts = [t for t in texts if not t.startswith("http")]
        first_text = texts[0] if texts else ""

        # Filter footer junk from nutrients. Bare 4-digit numbers without
        # commas or units are years (e.g. "2016"). Slashes are date parts.
        valid_nutrients = []
        for s in nutrients:
            if not NUTRIENT_RE.match(s):
                continue
            stripped = s.replace(",", "").replace("g", "").replace("m", "")
            if stripped.isdigit() and len(stripped) == 4 and "," not in s and not s.endswith("g"):
                continue
            valid_nutrients.append(s)

        # Meal header (text column only, no real nutrients)
        if first_text in MEALS and not valid_nutrients:
            # Flush any orphaned nutrients with accumulated name
            if orphaned_nutrients is not None and pending_name_parts:
                food_name = " ".join(pending_name_parts)
                if food_name and current_meal:
                    entry = {"date": current_date, "meal": current_meal, "food": food_name}
                    for j, col in enumerate(NUTRIENT_COLS):
                        entry[col] = parse_nutrient(orphaned_nutrients[j])
                    food_results.append(entry)
                orphaned_nutrients = None
                pending_name_parts = []
            current_meal = first_text
            pending_name_parts = []
            continue

        # --- Row has text + 8 nutrient columns: a complete food entry ---
        if len(valid_nutrients) == 8 and texts:
            # First, flush any orphaned nutrients with the name accumulated so far
            if orphaned_nutrients is not None:
                food_name = " ".join(pending_name_parts)
                if food_name and current_meal:
                    entry = {"date": current_date, "meal": current_meal, "food": food_name}
                    for j, col in enumerate(NUTRIENT_COLS):
                        entry[col] = parse_nutrient(orphaned_nutrients[j])
                    food_results.append(entry)
                orphaned_nutrients = None
                pending_name_parts = []

            # Now emit this row's food entry
            food_name_parts = pending_name_parts + texts
            if food_name_parts and food_name_parts[0] in MEALS:
                current_meal = food_name_parts[0]
                food_name_parts = food_name_parts[1:]
            food_name = " ".join(food_name_parts)
            pending_name_parts = []

            if food_name and current_meal:
                entry = {"date": current_date, "meal": current_meal, "food": food_name}
                for j, col in enumerate(NUTRIENT_COLS):
                    entry[col] = parse_nutrient(valid_nutrients[j])
                food_results.append(entry)

        # --- Row has ONLY nutrient columns: orphaned nutrients ---
        # The food name is split: some parts came before (in pending_name_parts),
        # and more parts may follow on subsequent text-only rows.
        elif len(valid_nutrients) == 8 and not texts:
            # If there were ALREADY orphaned nutrients, flush them first
            if orphaned_nutrients is not None:
                food_name = " ".join(pending_name_parts)
                if food_name and current_meal:
                    entry = {"date": current_date, "meal": current_meal, "food": food_name}
                    for j, col in enumerate(NUTRIENT_COLS):
                        entry[col] = parse_nutrient(orphaned_nutrients[j])
                    food_results.append(entry)
                pending_name_parts = []
            # Stash these nutrients — keep pending_name_parts as the name
            # collected so far, and continue accumulating text rows.
            orphaned_nutrients = list(valid_nutrients)

        # --- Row has only text: food name line or continuation ---
        elif texts and not nutrients:
            pending_name_parts.extend(texts)

        # --- Partial nutrients (not exercise anymore — unknown) — skip ---
        else:
            pending_name_parts = []

    return food_results, exercise_results


EXERCISE_COLS = ["date", "name", "calories", "minutes"]


def main():
    if len(sys.argv) < 2:
        print("Usage: extract_oxps.py <file_or_dir> [file_or_dir ...]", file=sys.stderr)
        print("  Extracts food + exercise data from OXPS files to CSV on stdout.", file=sys.stderr)
        sys.exit(1)

    all_food_rows = []
    all_exercise_rows = []
    paths = []
    for arg in sys.argv[1:]:
        p = Path(arg)
        if p.is_dir():
            paths.extend(sorted(p.rglob("*.oxps")))
        else:
            paths.append(p)

    for filepath in sorted(paths):
        try:
            food_rows, exercise_rows = extract_oxps(filepath)
            all_food_rows.extend(food_rows)
            all_exercise_rows.extend(exercise_rows)
        except Exception as e:
            print(f"ERROR processing {filepath}: {e}", file=sys.stderr)

    writer = csv.DictWriter(sys.stdout, fieldnames=["date", "meal", "food"] + NUTRIENT_COLS)
    writer.writeheader()
    for row in all_food_rows:
        writer.writerow(row)

    print("---EXERCISES---")

    writer = csv.DictWriter(sys.stdout, fieldnames=EXERCISE_COLS)
    writer.writeheader()
    for row in all_exercise_rows:
        writer.writerow(row)


if __name__ == "__main__":
    main()
