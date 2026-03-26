"""Extract body composition measurements into a normalized CSV.

Sources:
- InBody exports: every ``InBody-*.csv`` in this directory
- Static BOD POD rows: ``Body composition.xlsx`` (BOD POD rows only)

The extractor prefers InBody rows when dates overlap, and for multiple InBody
measurements on the same date it keeps the richest row, then the latest
timestamp. Output remains ``composition/composition.csv`` so downstream
analysis continues to read the same file.
"""

import csv
from datetime import datetime
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parent
XLSX = ROOT / "Body composition.xlsx"
OUT = ROOT / "composition.csv"
INBODY_GLOB = "InBody-*.csv"
LB_PER_LITER_WATER = 2.2046226218

# Column indices in the XLSX (1-based, from openpyxl)
# Row 1 = category headers, Row 2 = column names, Row 3+ = data
XLSX_COL = {
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

CORE_FIELDS = [
    "date",
    "measured_at",
    "source",
    "device",
    "era",
    "weight_lbs",
    "fat_mass_lbs",
    "fat_pct",
    "lean_mass_lbs",
    "dry_lean_mass_lbs",
    "soft_lean_mass_lbs",
    "smm_lbs",
    "bmi",
    "score",
    "bmr_kj",
]

WATER_FIELDS = [
    "intracellular_lbs",
    "extracellular_lbs",
    "water_lbs",
    "intracellular_water_l",
    "extracellular_water_l",
    "total_body_water_l",
    "visceral_fat_level",
    "visceral_fat_area_cm2",
    "ecw_tbw",
    "ecw_ratio",
]

SEGMENTAL_FIELDS = [
    "arm_muscle_r_lbs",
    "arm_muscle_l_lbs",
    "trunk_muscle_lbs",
    "leg_muscle_r_lbs",
    "leg_muscle_l_lbs",
    "arm_fat_r_lbs",
    "arm_fat_l_lbs",
    "trunk_fat_lbs",
    "leg_fat_r_lbs",
    "leg_fat_l_lbs",
    "arm_ecw_ratio_r",
    "arm_ecw_ratio_l",
    "trunk_ecw_ratio",
    "leg_ecw_ratio_r",
    "leg_ecw_ratio_l",
]

ADDITIONAL_FIELDS = [
    "waist_hip_ratio",
    "waist_circumference_in",
    "upper_lower",
    "upper_lbs",
    "lower_lbs",
    "leg_muscle_level",
    "leg_lean_mass_lbs",
    "protein_lbs",
    "mineral_lbs",
    "bone_mineral_content_lbs",
    "body_cell_mass_lbs",
    "smi_kg_m2",
    "phase_angle_deg",
]

ALL_FIELDS = CORE_FIELDS + WATER_FIELDS + SEGMENTAL_FIELDS + ADDITIONAL_FIELDS


def xlsx_cell(row, col_name):
    return row[XLSX_COL[col_name] - 1]


def to_float(value):
    if value in (None, "", "-"):
        return None
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


def fmt(value, digits=None):
    if value is None:
        return ""
    if digits is not None:
        value = round(value, digits)
    return value


def liters_to_lbs(value):
    liters = to_float(value)
    if liters is None:
        return None
    return liters * LB_PER_LITER_WATER


def parse_inbody_datetime(raw_value):
    text = (raw_value or "").strip()
    return datetime.strptime(text, "%Y%m%d%H%M%S")


def base_record():
    return {field: "" for field in ALL_FIELDS}


def classify_inbody_era(record):
    has_full_segmental = bool(record["arm_fat_r_lbs"] or record["trunk_fat_lbs"])
    has_water_compartments = bool(record["intracellular_water_l"] or record["extracellular_water_l"])
    has_total_water_only = bool(record["total_body_water_l"])

    if has_full_segmental or has_water_compartments:
        return "inbody_full"
    if has_total_water_only:
        return "inbody_partial"
    return "inbody_summary"


def record_richness(record):
    ignored = {"date", "measured_at", "source", "device", "era"}
    return sum(1 for key, value in record.items() if key not in ignored and value != "")


def better_record(candidate, incumbent):
    cand_key = (
        record_richness(candidate),
        candidate["measured_at"],
        candidate["source"],
    )
    inc_key = (
        record_richness(incumbent),
        incumbent["measured_at"],
        incumbent["source"],
    )
    return cand_key > inc_key


def load_bodpod_rows():
    wb = openpyxl.load_workbook(XLSX, data_only=True)
    ws = wb.active

    records = []
    for row in ws.iter_rows(min_row=3, values_only=True):
        date_val = row[XLSX_COL["date"] - 1]
        if date_val is None:
            continue

        weight = to_float(xlsx_cell(row, "weight"))
        fat_mass = to_float(xlsx_cell(row, "fat_mass"))
        if weight is None or fat_mass is None:
            continue

        # Keep only the static BOD POD rows from the workbook.
        if to_float(xlsx_cell(row, "smm")) is not None or to_float(xlsx_cell(row, "water")) is not None:
            continue

        measured_at = (
            date_val.replace(hour=0, minute=0, second=0, microsecond=0)
            if hasattr(date_val, "replace")
            else datetime.strptime(str(date_val)[:10], "%Y-%m-%d")
        )
        rec = base_record()
        rec.update(
            {
                "date": measured_at.date().isoformat(),
                "measured_at": measured_at.isoformat(),
                "source": XLSX.name,
                "device": "BOD POD",
                "era": "bodpod",
                "weight_lbs": fmt(weight, 1),
                "fat_mass_lbs": fmt(fat_mass, 1),
                "fat_pct": fmt(to_float(xlsx_cell(row, "fat_pct")), 1),
                "lean_mass_lbs": fmt(weight - fat_mass, 1),
            }
        )
        records.append(rec)

    return records


def load_inbody_exports():
    records_by_date = {}
    paths = sorted(ROOT.glob(INBODY_GLOB))
    if not paths:
        raise FileNotFoundError(f"No files matched {INBODY_GLOB!r} in {ROOT}")

    for path in paths:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                measured_at = parse_inbody_datetime(row["date"])
                rec = base_record()

                rec.update(
                    {
                        "date": measured_at.date().isoformat(),
                        "measured_at": measured_at.isoformat(),
                        "source": path.name,
                        "device": (row.get("Measurement device.") or "").strip(),
                        "weight_lbs": fmt(to_float(row.get("Weight(lb)")), 1),
                        "fat_mass_lbs": fmt(to_float(row.get("Body Fat Mass(lb)")), 1),
                        "fat_pct": fmt(to_float(row.get("Percent Body Fat(%)")), 1),
                        "dry_lean_mass_lbs": fmt(to_float(row.get("Soft Lean Mass(lb)")), 1),
                        "soft_lean_mass_lbs": fmt(to_float(row.get("Soft Lean Mass(lb)")), 1),
                        "smm_lbs": fmt(to_float(row.get("Skeletal Muscle Mass(lb)")), 1),
                        "bmi": fmt(to_float(row.get("BMI(kg/m²)")), 1),
                        "score": fmt(to_float(row.get("InBody Score")), 1),
                        "bmr_kj": fmt(to_float(row.get("Basal Metabolic Rate(kJ)")), 0),
                        "intracellular_lbs": fmt(liters_to_lbs(row.get("Intracellular Water(L)")), 1),
                        "extracellular_lbs": fmt(liters_to_lbs(row.get("Extracellular Water(L)")), 1),
                        "water_lbs": fmt(liters_to_lbs(row.get("Total Body Water(L)")), 1),
                        "intracellular_water_l": fmt(to_float(row.get("Intracellular Water(L)")), 1),
                        "extracellular_water_l": fmt(to_float(row.get("Extracellular Water(L)")), 1),
                        "total_body_water_l": fmt(to_float(row.get("Total Body Water(L)")), 1),
                        "visceral_fat_level": fmt(to_float(row.get("Visceral Fat Level(Level)")), 1),
                        "visceral_fat_area_cm2": fmt(to_float(row.get("Visceral Fat Area(cm²)")), 1),
                        "ecw_tbw": fmt(to_float(row.get("ECW Ratio")), 3),
                        "ecw_ratio": fmt(to_float(row.get("ECW Ratio")), 3),
                        "arm_muscle_r_lbs": fmt(to_float(row.get("Right Arm Lean Mass(lb)")), 2),
                        "arm_muscle_l_lbs": fmt(to_float(row.get("Left Arm Lean Mass(lb)")), 2),
                        "trunk_muscle_lbs": fmt(to_float(row.get("Trunk Lean Mass(lb)")), 1),
                        "leg_muscle_r_lbs": fmt(to_float(row.get("Right Leg Lean Mass(lb)")), 2),
                        "leg_muscle_l_lbs": fmt(to_float(row.get("Left leg Lean Mass(lb)")), 2),
                        "arm_fat_r_lbs": fmt(to_float(row.get("Right Arm Fat Mass(lb)")), 1),
                        "arm_fat_l_lbs": fmt(to_float(row.get("Left Arm Fat Mass(lb)")), 1),
                        "trunk_fat_lbs": fmt(to_float(row.get("Trunk Fat Mass(lb)")), 1),
                        "leg_fat_r_lbs": fmt(to_float(row.get("Right Leg Fat Mass(lb)")), 1),
                        "leg_fat_l_lbs": fmt(to_float(row.get("Left Leg Fat Mass(lb)")), 1),
                        "arm_ecw_ratio_r": fmt(to_float(row.get("Right Arm ECW Ratio")), 3),
                        "arm_ecw_ratio_l": fmt(to_float(row.get("Left Arm ECW Ratio")), 3),
                        "trunk_ecw_ratio": fmt(to_float(row.get("Trunk ECW Ratio")), 3),
                        "leg_ecw_ratio_r": fmt(to_float(row.get("Right Leg ECW Ratio")), 3),
                        "leg_ecw_ratio_l": fmt(to_float(row.get("Left Leg ECW Ratio")), 3),
                        "waist_hip_ratio": fmt(to_float(row.get("Waist Hip Ratio")), 2),
                        "waist_circumference_in": fmt(to_float(row.get("Waist Circumference(inch)")), 1),
                        "upper_lower": fmt(to_float(row.get("Upper-Lower")), 2),
                        "upper_lbs": fmt(to_float(row.get("Upper")), 1),
                        "lower_lbs": fmt(to_float(row.get("Lower")), 1),
                        "leg_muscle_level": fmt(to_float(row.get("Leg Muscle Level(Level)")), 1),
                        "leg_lean_mass_lbs": fmt(to_float(row.get("Leg Lean Mass(lb)")), 1),
                        "protein_lbs": fmt(to_float(row.get("Protein(lb)")), 1),
                        "mineral_lbs": fmt(to_float(row.get("Mineral(lb)")), 1),
                        "bone_mineral_content_lbs": fmt(to_float(row.get("Bone Mineral Content(lb)")), 1),
                        "body_cell_mass_lbs": fmt(to_float(row.get("Body Cell Mass(lb)")), 1),
                        "smi_kg_m2": fmt(to_float(row.get("SMI(kg/m²)")), 1),
                        "phase_angle_deg": fmt(to_float(row.get("Whole Body Phase Angle(°)")), 1),
                    }
                )

                weight = to_float(row.get("Weight(lb)"))
                fat_mass = to_float(row.get("Body Fat Mass(lb)"))
                if weight is not None and fat_mass is not None:
                    rec["lean_mass_lbs"] = fmt(weight - fat_mass, 1)

                rec["era"] = classify_inbody_era(rec)

                incumbent = records_by_date.get(rec["date"])
                if incumbent is None or better_record(rec, incumbent):
                    records_by_date[rec["date"]] = rec

    return list(records_by_date.values())


def extract():
    inbody_records = load_inbody_exports()
    bodpod_records = load_bodpod_rows()

    records_by_date = {record["date"]: record for record in inbody_records}
    for record in bodpod_records:
        records_by_date.setdefault(record["date"], record)

    records = sorted(records_by_date.values(), key=lambda record: record["measured_at"])

    eras = {}
    for record in records:
        eras[record["era"]] = eras.get(record["era"], 0) + 1

    with OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ALL_FIELDS)
        writer.writeheader()
        writer.writerows(records)

    print(f"Read {len(inbody_records)} unique InBody dates from {INBODY_GLOB}")
    print(f"Read {len(bodpod_records)} BOD POD dates from {XLSX.name}")
    print(f"Wrote {len(records)} measurements to {OUT.name}")
    print(f"  Date range: {records[0]['date']} to {records[-1]['date']}")
    for era, count in sorted(eras.items()):
        print(f"  {era}: {count}")


if __name__ == "__main__":
    extract()
