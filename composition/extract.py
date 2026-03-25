"""Extract body composition measurements from XLSX to CSV.

Handles 3 measurement eras:
- BOD POD (air displacement): date, weight, fat_mass, fat_pct only
- InBody Partial (2017-03 to 2017-10): adds water, lean_mass, smm, bmi, muscle segmental
- InBody Full (2017-11 to 2019-11): adds water compartments, visceral fat, fat segmental

Idempotent. Re-run to regenerate CSV from XLSX.

Output: composition/composition.csv
"""

import csv
from pathlib import Path
import openpyxl

ROOT = Path(__file__).resolve().parent
XLSX = ROOT / "Body composition.xlsx"
OUT = ROOT / "composition.csv"

# Column indices in the XLSX (1-based, from openpyxl)
# Row 1 = category headers, Row 2 = column names, Row 3+ = data
COL = {
    "date": 1,
    "intracellular": 2,
    "extracellular": 3,
    "water": 4,
    "lean_mass": 5,
    "fat_mass": 6,
    "weight": 7,
    "smm": 8,
    "bmi": 9,
    "fat_pct": 10,
    "arm_muscle_r": 11,
    "arm_muscle_l": 12,
    "trunk_muscle": 13,
    "leg_muscle_r": 14,
    "leg_muscle_l": 15,
    "ecw_tbw": 16,
    "visceral": 17,
    "arm_fat_r": 18,
    "arm_fat_l": 19,
    "trunk_fat": 20,
    "leg_fat_r": 21,
    "leg_fat_l": 22,
    "score": 23,
}

# Output fields (core + extended)
CORE_FIELDS = [
    "date", "weight_lbs", "fat_mass_lbs", "fat_pct", "lean_mass_lbs",
    "smm_lbs", "bmi", "era",
]
EXTENDED_FIELDS = [
    "intracellular_lbs", "extracellular_lbs", "water_lbs",
    "visceral_fat_level", "ecw_tbw",
    "arm_muscle_r_lbs", "arm_muscle_l_lbs", "trunk_muscle_lbs",
    "leg_muscle_r_lbs", "leg_muscle_l_lbs",
    "arm_fat_r_lbs", "arm_fat_l_lbs", "trunk_fat_lbs",
    "leg_fat_r_lbs", "leg_fat_l_lbs",
]
ALL_FIELDS = CORE_FIELDS + EXTENDED_FIELDS


def cell(row, col_name):
    """Get cell value by column name."""
    return row[COL[col_name] - 1]


def to_float(v):
    """Convert cell to float, returning None for empty/non-numeric."""
    if v is None:
        return None
    try:
        return float(v)
    except (ValueError, TypeError):
        return None


def classify_era(row_data):
    """Classify measurement era based on which columns are populated."""
    has_visceral = to_float(cell(row_data, "visceral")) is not None
    has_lean = to_float(cell(row_data, "lean_mass")) is not None

    if has_visceral:
        return "inbody_full"
    elif has_lean:
        return "inbody_partial"
    else:
        return "bodpod"


def extract():
    wb = openpyxl.load_workbook(XLSX, data_only=True)
    ws = wb.active

    rows = list(ws.iter_rows(min_row=3, values_only=True))  # skip 2 header rows
    print(f"Read {len(rows)} rows from {XLSX.name}")

    records = []
    for row in rows:
        date_val = row[COL["date"] - 1]
        if date_val is None:
            continue

        date_str = date_val.strftime("%Y-%m-%d") if hasattr(date_val, "strftime") else str(date_val)[:10]

        weight = to_float(cell(row, "weight"))
        fat_mass = to_float(cell(row, "fat_mass"))
        fat_pct = to_float(cell(row, "fat_pct"))

        if weight is None or fat_mass is None:
            print(f"  Skipping {date_str}: missing weight or fat_mass")
            continue

        era = classify_era(row)

        # Lean mass: use InBody value if available, else derive from weight - fat_mass
        lean_mass = to_float(cell(row, "lean_mass"))
        if lean_mass is None:
            lean_mass = round(weight - fat_mass, 1)

        rec = {
            "date": date_str,
            "weight_lbs": weight,
            "fat_mass_lbs": fat_mass,
            "fat_pct": fat_pct,
            "lean_mass_lbs": lean_mass,
            "smm_lbs": to_float(cell(row, "smm")) or "",
            "bmi": to_float(cell(row, "bmi")) or "",
            "era": era,
            # Extended
            "intracellular_lbs": to_float(cell(row, "intracellular")) or "",
            "extracellular_lbs": to_float(cell(row, "extracellular")) or "",
            "water_lbs": to_float(cell(row, "water")) or "",
            "visceral_fat_level": to_float(cell(row, "visceral")) or "",
            "ecw_tbw": to_float(cell(row, "ecw_tbw")) or "",
            "arm_muscle_r_lbs": to_float(cell(row, "arm_muscle_r")) or "",
            "arm_muscle_l_lbs": to_float(cell(row, "arm_muscle_l")) or "",
            "trunk_muscle_lbs": to_float(cell(row, "trunk_muscle")) or "",
            "leg_muscle_r_lbs": to_float(cell(row, "leg_muscle_r")) or "",
            "leg_muscle_l_lbs": to_float(cell(row, "leg_muscle_l")) or "",
            "arm_fat_r_lbs": to_float(cell(row, "arm_fat_r")) or "",
            "arm_fat_l_lbs": to_float(cell(row, "arm_fat_l")) or "",
            "trunk_fat_lbs": to_float(cell(row, "trunk_fat")) or "",
            "leg_fat_r_lbs": to_float(cell(row, "leg_fat_r")) or "",
            "leg_fat_l_lbs": to_float(cell(row, "leg_fat_l")) or "",
        }
        records.append(rec)

    # Sort by date
    records.sort(key=lambda r: r["date"])

    # Count eras
    eras = {}
    for r in records:
        eras[r["era"]] = eras.get(r["era"], 0) + 1

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ALL_FIELDS)
        w.writeheader()
        w.writerows(records)

    print(f"Wrote {len(records)} measurements to {OUT.name}")
    print(f"  Date range: {records[0]['date']} to {records[-1]['date']}")
    for era, count in sorted(eras.items()):
        print(f"  {era}: {count}")


if __name__ == "__main__":
    extract()
