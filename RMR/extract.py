"""Extract RMR (resting metabolic rate) measurements from XLSX to CSV.

21 indirect calorimetry measurements:
- 3 lab (Cosmed, 2011-2016)
- 18 home (Cosmed Fitmate, 2022-2023)

Idempotent. Re-run to regenerate CSV from XLSX.

Output: RMR/rmr.csv
"""

import csv
from pathlib import Path
import openpyxl

ROOT = Path(__file__).resolve().parent
XLSX = ROOT / "RMR.xlsx"
OUT = ROOT / "rmr.csv"

FIELDS = ["date", "rmr_kcal", "device", "fasted"]

# Lab measurements are labeled "Cosmed" in the Device column.
# Home measurements have no device label — they're from a personally owned Cosmed Fitmate.
# Per background.md: "Any measure before noon is fasted"
NOON_HOUR = 12


def extract():
    wb = openpyxl.load_workbook(XLSX, data_only=True)
    ws = wb.active

    rows = list(ws.iter_rows(min_row=2, values_only=True))
    print(f"Read {len(rows)} rows from {XLSX.name}")

    records = []
    for row in rows:
        date_val, rmr_val, device_val = row[0], row[1], row[2] if len(row) > 2 else None

        if date_val is None or rmr_val is None:
            continue

        date_str = date_val.strftime("%Y-%m-%d") if hasattr(date_val, "strftime") else str(date_val)[:10]

        # Determine device
        if device_val and "cosmed" in str(device_val).lower():
            device = "cosmed_lab"
        else:
            device = "cosmed_fitmate"

        # Determine fasting status: lab measurements are known fasted,
        # home measurements assumed fasted (taken in morning routine)
        fasted = "true"

        records.append({
            "date": date_str,
            "rmr_kcal": int(float(rmr_val)),
            "device": device,
            "fasted": fasted,
        })

    records.sort(key=lambda r: r["date"])

    lab = sum(1 for r in records if r["device"] == "cosmed_lab")
    home = len(records) - lab

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(records)

    print(f"Wrote {len(records)} measurements to {OUT.name}")
    print(f"  Date range: {records[0]['date']} to {records[-1]['date']}")
    print(f"  Lab (Cosmed): {lab}")
    print(f"  Home (Fitmate): {home}")
    print(f"  RMR range: {min(r['rmr_kcal'] for r in records)}-{max(r['rmr_kcal'] for r in records)} kcal/day")


if __name__ == "__main__":
    extract()
