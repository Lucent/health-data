#!/usr/bin/env python3
"""Extract food diary data from HTML and MHTML files (MyFitnessPal printable diary).

Handles three layout eras:
  - Old layout (2011 Q2 through 2022-09): <h2 id="date">, <td> cells
  - New MUI layout (2022-10 through 2026+): <p> date tags, <th scope="row"> cells
  - MHTML files: quoted-printable decoded, then parsed same as HTML

Output: food CSV then '---EXERCISES---' sentinel then exercise CSV, both to stdout.
  Food:     date, meal, food, calories, carbs_g, fat_g, protein_g, cholest_mg, sodium_mg, sugars_g, fiber_g
  Exercise: date, name, calories, minutes
"""

import csv
import html as html_mod
import io
import quopri
import re
import sys
from pathlib import Path

NUTRIENT_COLS = ["calories", "carbs_g", "fat_g", "protein_g",
                 "cholest_mg", "sodium_mg", "sugars_g", "fiber_g"]

MEALS = {"Breakfast", "Lunch", "Dinner", "Snacks", "Supper"}

EXERCISE_CATEGORIES = {"Cardiovascular", "Strength"}

MONTH_MAP = {
    "January": "01", "February": "02", "March": "03", "April": "04",
    "May": "05", "June": "06", "July": "07", "August": "08",
    "September": "09", "October": "10", "November": "11", "December": "12",
    "Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04",
    "Jun": "06", "Jul": "07", "Aug": "08",
    "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12",
}

# "January 1, 2013" or "Jan 1, 2026"
DATE_RE = re.compile(
    r"^(January|February|March|April|May|June|July|August|September|"
    r"October|November|December|Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)"
    r"\s+(\d{1,2}),\s+(\d{4})$"
)


def parse_date(s):
    """Convert 'January 1, 2013' or 'Jan 1, 2026' to '2013-01-01'."""
    m = DATE_RE.match(s.strip())
    if not m:
        return None
    month_name, day, year = m.group(1), int(m.group(2)), m.group(3)
    return f"{year}-{MONTH_MAP[month_name]}-{int(day):02d}"


def parse_nutrient(s):
    """Normalize a nutrient string."""
    s = s.strip()
    if s == "--" or s.startswith("--") or s == "":
        return ""
    s = s.replace("mg", "").replace("g", "").replace(",", "")
    return s


def decode_mhtml(filepath):
    """Extract the HTML part from an MHTML file, decoding quoted-printable."""
    with open(filepath, "rb") as f:
        raw = f.read()

    # Find boundary
    text_raw = raw.decode("utf-8", errors="replace")

    # Find the first text/html part
    parts = []
    boundary_match = re.search(r'boundary="([^"]+)"', text_raw)
    if not boundary_match:
        # Not actually MHTML, just decode the whole thing
        inp = io.BytesIO(raw)
        out = io.BytesIO()
        quopri.decode(inp, out)
        return out.getvalue().decode("utf-8", errors="replace")

    boundary = boundary_match.group(1)
    sections = raw.split(boundary.encode())

    for section in sections:
        section_str = section.decode("utf-8", errors="replace")
        if "text/html" in section_str[:500]:
            # Find where headers end and content begins
            header_end = section_str.find("\r\n\r\n")
            if header_end == -1:
                header_end = section_str.find("\n\n")
            if header_end == -1:
                continue

            headers = section_str[:header_end]
            content_bytes = section[header_end + 4:] if b"\r\n\r\n" in section[:header_end + 10] else section[header_end + 2:]

            if "quoted-printable" in headers.lower():
                inp = io.BytesIO(content_bytes)
                out = io.BytesIO()
                quopri.decode(inp, out)
                return out.getvalue().decode("utf-8", errors="replace")
            elif "base64" in headers.lower():
                import base64
                return base64.b64decode(content_bytes).decode("utf-8", errors="replace")
            else:
                return content_bytes.decode("utf-8", errors="replace")

    # Fallback: decode whole file
    inp = io.BytesIO(raw)
    out = io.BytesIO()
    quopri.decode(inp, out)
    return out.getvalue().decode("utf-8", errors="replace")


def detect_layout(html):
    """Detect whether this is old layout or new MUI layout."""
    if 'id="date"' in html:
        return "old"
    if "MuiTable" in html or "MuiTypography" in html:
        return "new"
    # Fallback: check for old-style table
    if 'class="table0"' in html or 'id="food"' in html:
        return "old"
    return "new"


def extract_old_layout(html):
    """Parse old MyFitnessPal layout (2011-2022-09).

    Structure:
      <h2 class="main-title-2" id="date">January 1, 2012</h2>
      <table id="food">
        <tr class="title"><td colspan="9">Breakfast</td></tr>
        <tr><td class="first">Food name</td><td>120</td>...<td class="last">2g</td></tr>
        ...
        <tr class="total"><td class="first">TOTAL:</td><td>1234</td>...</tr>
      </table>
    """
    rows = []
    current_date = None
    current_meal = None

    # Split by date headers
    # Match <h2...id="date">DATE</h2>
    date_splits = re.split(r'<h2[^>]*id="date"[^>]*>', html)

    for chunk in date_splits[1:]:  # skip everything before first date
        # Extract date
        date_end = chunk.find("</h2>")
        if date_end == -1:
            continue
        date_str = chunk[:date_end].strip()
        date_str = html_mod.unescape(re.sub(r"<[^>]+>", "", date_str)).strip()
        current_date = parse_date(date_str)
        if not current_date:
            continue

        # Find the food table
        food_table_match = re.search(r'<table[^>]*id="food"[^>]*>(.*?)</table>',
                                      chunk, re.DOTALL)
        if not food_table_match:
            continue
        table_html = food_table_match.group(1)

        # Process rows
        tr_chunks = re.split(r"<tr[^>]*>", table_html)
        for tr in tr_chunks:
            # Strip closing tag
            tr = tr.split("</tr>")[0] if "</tr>" in tr else tr

            # Extract all cell contents (td or th)
            cells = re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", tr, re.DOTALL)
            cells = [html_mod.unescape(re.sub(r"<[^>]+>", "", c)).strip() for c in cells]
            cells = [c for c in cells if c != ""]

            if not cells:
                continue

            # Check for meal header (single cell spanning all columns)
            if len(cells) == 1 and cells[0] in MEALS:
                current_meal = cells[0]
                continue

            # Check for column headers
            if cells[0] in ("Foods", "FOODS"):
                continue

            # Check for TOTAL row
            if cells[0] in ("TOTAL:", "TOTALS:"):
                if len(cells) >= 9:
                    row = {
                        "date": current_date,
                        "meal": "TOTAL",
                        "food": "TOTAL",
                    }
                    for j, col in enumerate(NUTRIENT_COLS):
                        row[col] = parse_nutrient(cells[1 + j])
                    rows.append(row)
                continue

            # Food row: food name + 8 nutrients
            if len(cells) >= 9 and current_meal:
                food_name = cells[0]
                row = {
                    "date": current_date,
                    "meal": current_meal,
                    "food": food_name,
                }
                for j, col in enumerate(NUTRIENT_COLS):
                    row[col] = parse_nutrient(cells[1 + j])
                rows.append(row)

    return rows


def extract_new_layout(html):
    """Parse new MUI MyFitnessPal layout (2022-10 onward).

    Structure:
      <p class="MuiTypography-root ...">Jan 1, 2026</p>
      <div ...>
        <table>
          <tr><th>FOODS</th><th>Calories</th>...</tr>
          <td colspan="9">Breakfast</td>
          <tr><th scope="row">Food name</th><th scope="row">120</th>...</tr>
          ...
          <tr class="...total..."><th>TOTALS</th><th>1234</th>...</tr>
        </table>
      </div>
    """
    rows = []
    current_date = None
    current_meal = None

    # Find all dates. They appear as text in <p> tags, or in aria-labels.
    # Strategy: find date strings and split the document by them.
    date_pattern = (
        r"(?:January|February|March|April|May|June|July|August|September|"
        r"October|November|December|Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)"
        r"\s+\d{1,2},\s+\d{4}"
    )

    # Find dates that appear as content (not in aria-label attributes, which are duplicates)
    # Use a different approach: find <p> tags or text nodes containing dates
    # followed by table content

    # Split by date paragraphs - look for the date followed by table data
    # The key pattern: date text appears, then food table follows
    pieces = re.split(r"(?=>(?:" + date_pattern + r")</)", html)

    # Better approach: find all date occurrences and their positions
    date_positions = []
    for m in re.finditer(r">(" + date_pattern + r")<", html):
        date_str = m.group(1)
        # Skip dates that are in aria-label (they appear as "selected date is ...")
        context_before = html[max(0, m.start() - 100):m.start()]
        if "selected date" in context_before or "aria-label" in context_before:
            continue
        parsed = parse_date(date_str)
        if parsed:
            date_positions.append((m.end(), parsed))

    if not date_positions:
        return rows

    # Process each date section
    for idx, (pos, date_val) in enumerate(date_positions):
        current_date = date_val
        current_meal = None

        # Get HTML chunk until next date.
        # Use the start of the >DATE< match (not -200 fudge) to avoid
        # cutting off the TOTAL row that sits just before the next date.
        if idx + 1 < len(date_positions):
            # Find the start of the date's container element (the > before the date text)
            next_date_end = date_positions[idx + 1][0]
            # Search backwards from the date text for a structural break
            # The date sits inside a <p> or similar — find the opening < before it
            search_start = max(pos, next_date_end - 500)
            # Find last </table> or </div> before the next date
            table_end = html.rfind("</table>", search_start, next_date_end)
            div_end = html.rfind("</div>", search_start, next_date_end)
            end_pos = max(table_end, div_end)
            if end_pos <= pos:
                end_pos = next_date_end  # fallback
        else:
            end_pos = len(html)
        chunk = html[pos:end_pos]

        # Find all table rows in this chunk
        # Extract th and td content from tr elements
        tr_chunks = re.split(r"<tr[^>]*>", chunk)

        for tr in tr_chunks:
            tr = tr.split("</tr>")[0] if "</tr>" in tr else tr

            # Check for meal header:
            # Old MUI: <td colspan="9">Breakfast</td>
            # New MUI (2024+): <td>Breakfast</td><td colspan="9"></td>
            meal_match = re.search(r'colspan="9"[^>]*>(.*?)</td>', tr, re.DOTALL)
            if meal_match:
                meal_text = html_mod.unescape(re.sub(r"<[^>]+>", "", meal_match.group(1))).strip()
                if meal_text in MEALS:
                    current_meal = meal_text
                    continue
                # Empty colspan cell — meal name might be in a prior td
                if not meal_text:
                    all_td = re.findall(r"<td[^>]*>(.*?)</td>", tr, re.DOTALL)
                    for td_content in all_td:
                        td_text = html_mod.unescape(re.sub(r"<[^>]+>", "", td_content)).strip()
                        if td_text in MEALS:
                            current_meal = td_text
                            break
                    continue

            # Extract cells (th with scope="row" or just th/td)
            cells = re.findall(r"<th[^>]*>(.*?)</th>", tr, re.DOTALL)
            if not cells:
                cells = re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", tr, re.DOTALL)
            cells = [html_mod.unescape(re.sub(r"<[^>]+>", "", c)).strip() for c in cells]
            cells = [c for c in cells if c != ""]

            if not cells:
                continue

            # Column headers
            if cells[0] in ("Foods", "FOODS"):
                continue

            # TOTAL row
            if cells[0] in ("TOTAL:", "TOTALS:", "TOTAL", "TOTALS"):
                if len(cells) >= 9:
                    row = {
                        "date": current_date,
                        "meal": "TOTAL",
                        "food": "TOTAL",
                    }
                    for j, col in enumerate(NUTRIENT_COLS):
                        row[col] = parse_nutrient(cells[1 + j])
                    rows.append(row)
                continue

            # Food row
            if len(cells) >= 9 and current_meal:
                food_name = cells[0]
                row = {
                    "date": current_date,
                    "meal": current_meal,
                    "food": food_name,
                }
                for j, col in enumerate(NUTRIENT_COLS):
                    row[col] = parse_nutrient(cells[1 + j])
                rows.append(row)

    return rows


def extract_exercises_old_layout(html):
    """Parse exercise data from old MyFitnessPal layout (2011-2022-09).

    Structure:
      <h2 class="main-title-2" id="date">January 1, 2012</h2>
      ... food table ...
      <table id="excercise">                (note: MFP misspelling)
        <thead><tr><td>Exercises</td><td>Calories</td><td>Minutes</td>...</tr></thead>
        <tbody>
          <tr class="title"><td colspan="6">Cardiovascular</td></tr>
          <tr><td>Hospital</td><td>300</td><td>240</td>...</tr>
        </tbody>
        <tfoot><tr><td>TOTALS:</td>...</tr></tfoot>
      </table>
    """
    rows = []
    current_date = None

    # Split by date headers
    date_splits = re.split(r'<h2[^>]*id="date"[^>]*>', html)

    for chunk in date_splits[1:]:
        date_end = chunk.find("</h2>")
        if date_end == -1:
            continue
        date_str = chunk[:date_end].strip()
        date_str = html_mod.unescape(re.sub(r"<[^>]+>", "", date_str)).strip()
        current_date = parse_date(date_str)
        if not current_date:
            continue

        # Find exercise table(s) — note MFP misspelling "excercise"
        for ex_match in re.finditer(
            r'<table[^>]*id="excercise"[^>]*>(.*?)</table>', chunk, re.DOTALL
        ):
            table_html = ex_match.group(1)
            tr_chunks = re.split(r"<tr[^>]*>", table_html)
            for tr in tr_chunks:
                tr = tr.split("</tr>")[0] if "</tr>" in tr else tr
                cells = re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", tr, re.DOTALL)
                cells = [html_mod.unescape(re.sub(r"<[^>]+>", "", c)).strip()
                         for c in cells]
                cells = [c for c in cells if c not in ("", "&nbsp;", "\xa0")]

                if not cells:
                    continue
                # Skip header row, category headers, TOTALS
                if cells[0] in ("Exercises", "EXERCISES"):
                    continue
                if cells[0] in EXERCISE_CATEGORIES:
                    continue
                if cells[0] in ("TOTAL:", "TOTALS:"):
                    continue
                # Exercise row: name, calories, minutes [, sets, reps, weight]
                if len(cells) >= 3:
                    rows.append({
                        "date": current_date,
                        "name": cells[0],
                        "calories": cells[1],
                        "minutes": cells[2],
                    })

    return rows


def extract_file(filepath):
    """Extract food and exercise data from an HTML or MHTML file.

    Returns (food_rows, exercise_rows).
    """
    filepath = Path(filepath)

    if filepath.suffix.lower() == ".mhtml":
        html = decode_mhtml(filepath)
    else:
        with open(filepath, "r", encoding="utf-8", errors="replace") as f:
            html = f.read()

    layout = detect_layout(html)
    if layout == "old":
        return extract_old_layout(html), extract_exercises_old_layout(html)
    else:
        return extract_new_layout(html), []


EXERCISE_COLS = ["date", "name", "calories", "minutes"]


def main():
    if len(sys.argv) < 2:
        print("Usage: extract_html.py <file_or_dir> [file_or_dir ...]", file=sys.stderr)
        sys.exit(1)

    all_food_rows = []
    all_exercise_rows = []
    paths = []
    for arg in sys.argv[1:]:
        p = Path(arg)
        if p.is_dir():
            for ext in ("*.html", "*.mhtml"):
                paths.extend(sorted(p.rglob(ext)))
        else:
            paths.append(p)

    # Filter out "new layout" duplicate files
    paths = [p for p in paths if "new layout" not in p.name]
    paths = sorted(paths)

    for filepath in paths:
        try:
            food_rows, exercise_rows = extract_file(filepath)
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
