#!/usr/bin/env python3
"""Extract canonical RMR measurements for patient 1."""

from __future__ import annotations

import csv
from pathlib import Path

import openpyxl

from fitmate_dump import TIPO_RMR, open_table, parse_rmr_blob, compute_rmr_stats


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "rmr_fitmate_user1.csv"
OUT = ROOT / "rmr.csv"
XLSX = ROOT / "RMR.xlsx"
PATIENT_ID = 1

FIELDS = [
    "date",
    "rmr_kcal",
    "device",
    "fasted",
    "vo2_mL_min",
    "ve_L_min",
    "rf_br_min",
    "feo2_pct",
    "cv_ve_pct",
    "cv_vo2_pct",
    "n_breaths",
    "duration_min",
]
def load_patient_name() -> str:
    for rec in open_table(str(ROOT), "ANAGRAFE.DBF"):
        row = dict(rec)
        if row["PROGRESS"] == PATIENT_ID:
            return f"{row['A_LASTNAME']}, {row['A_FRSTNAME']}"
    raise SystemExit("Patient 1 not found in ANAGRAFE.DBF")


def extract_patient_tests() -> list[dict[str, object]]:
    records = []
    for rec in open_table(str(ROOT), "TEST.DBF"):
        row = dict(rec)
        if row["T_ANAG"] != PATIENT_ID or row["T_TIPO"] != TIPO_RMR:
            continue
        parsed = parse_rmr_blob(row.get("T_TEST"))
        if parsed["rmr_kcal"] == 0:
            continue
        stats = compute_rmr_stats(parsed)
        rec = {
            "date": row["T_DATE"].isoformat(),
            "rmr_kcal": int(parsed["rmr_kcal"]),
            "device": "cosmed_fitmate",
            "fasted": True,
        }
        for k in ("vo2_mL_min", "ve_L_min", "rf_br_min", "feo2_pct",
                   "cv_ve_pct", "cv_vo2_pct",
                   "n_breaths", "duration_min"):
            rec[k] = stats.get(k, "")
        records.append(rec)
    records.sort(key=lambda row: row["date"])
    return records


def load_xlsx_rows() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    wb = openpyxl.load_workbook(XLSX, data_only=True)
    ws = wb.active
    lab_rows = []
    other_rows = []
    for date_val, rmr_val, device_val, *_ in ws.iter_rows(min_row=2, values_only=True):
        if date_val is None or rmr_val is None:
            continue
        rmr_kcal = int(float(rmr_val))
        if rmr_kcal == 0:
            continue
        device_str = "" if device_val is None else str(device_val).strip().lower()
        row = {
            "date": date_val.strftime("%Y-%m-%d"),
            "rmr_kcal": rmr_kcal,
            "device": "cosmed_lab" if device_str == "cosmed" else "cosmed_fitmate",
            "fasted": True,
        }
        if device_str == "cosmed":
            lab_rows.append(row)
        else:
            other_rows.append(row)
    lab_rows.sort(key=lambda row: row["date"])
    other_rows.sort(key=lambda row: row["date"])
    return lab_rows, other_rows


def compare_to_xlsx(dbf_rows: list[dict[str, object]], xlsx_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    dbf_map = {row["date"]: row["rmr_kcal"] for row in dbf_rows}
    xlsx_map = {row["date"]: row["rmr_kcal"] for row in xlsx_rows}

    all_dates = sorted(set(dbf_map) | set(xlsx_map))
    compare_rows = []
    for date in all_dates:
        dbf_val = dbf_map.get(date)
        xlsx_val = xlsx_map.get(date)
        if dbf_val is not None and xlsx_val is not None:
            if dbf_val == xlsx_val:
                status = "match"
            else:
                status = "value_mismatch"
        elif dbf_val is not None:
            status = "dbf_only"
        else:
            status = "xlsx_only"
        compare_rows.append(
            {
                "date": date,
                "dbf_rmr_kcal": "" if dbf_val is None else dbf_val,
                "xlsx_rmr_kcal": "" if xlsx_val is None else xlsx_val,
                "status": status,
            }
        )
    return compare_rows


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    patient_name = load_patient_name()
    dbf_rows = extract_patient_tests()
    lab_rows, xlsx_fitmate_rows = load_xlsx_rows()
    compare_rows = compare_to_xlsx(dbf_rows, xlsx_fitmate_rows)

    dbf_dates = {row["date"] for row in dbf_rows}
    manual_fitmate_rows = [row for row in xlsx_fitmate_rows if row["date"] not in dbf_dates]
    combined_rows = sorted(lab_rows + dbf_rows + manual_fitmate_rows, key=lambda row: row["date"])

    write_csv(OUT, combined_rows, FIELDS)

    statuses = {}
    for row in compare_rows:
        statuses[row["status"]] = statuses.get(row["status"], 0) + 1

    xlsx_dates = {row["date"] for row in xlsx_fitmate_rows}
    is_superset = xlsx_dates.issubset(dbf_dates) and statuses.get("value_mismatch", 0) == 0

    print(f"Patient: {patient_name} [1]")
    print(f"DBF rows: {len(dbf_rows)}")
    print(f"XLSX fitmate rows: {len(xlsx_fitmate_rows)}")
    print(f"Lab rows from XLSX: {len(lab_rows)}")
    print(f"Manual XLSX-only fitmate rows kept: {len(manual_fitmate_rows)}")
    print(f"Date range: {combined_rows[0]['date']} to {combined_rows[-1]['date']}")
    print(f"Superset of XLSX fitmate rows: {'yes' if is_superset else 'no'}")
    for status in ["match", "dbf_only", "xlsx_only", "value_mismatch"]:
        if statuses.get(status):
            print(f"  {status}: {statuses[status]}")
    print(f"Wrote: {OUT}")


if __name__ == "__main__":
    main()
