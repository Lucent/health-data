#!/usr/bin/env python3
"""Extract daily weight from manual-entry.xlsx + history.csv (WeightGurus).

manual-entry.xlsx: pre-WeightGurus historical data (1999–2018).
history.csv:       WeightGurus app export (2019+). Drop a fresh export any time.

Takes the first reading per day (most standardized: fasted, post-sleep).

Output: weight.csv — date, weight_lbs, time (if available)
"""

import csv
import openpyxl
from collections import defaultdict
from datetime import datetime
from pathlib import Path

WEIGHT_DIR = Path(__file__).parent


def main():
    by_day = defaultdict(list)

    # 1. Manual-entry xlsx (historical)
    xlsx_path = WEIGHT_DIR / "manual-entry.xlsx"
    if xlsx_path.exists():
        wb = openpyxl.load_workbook(xlsx_path, data_only=True)
        ws = wb["weight"]
        for row in ws.iter_rows(min_row=2, values_only=True):
            dt, weight = row[0], row[1]
            if dt is None or weight is None:
                continue
            date_str = dt.strftime("%Y-%m-%d")
            has_time = dt.hour != 0 or dt.minute != 0
            time_str = dt.strftime("%H:%M") if has_time else ""
            by_day[date_str].append((dt, weight, time_str))

    # 2. WeightGurus history.csv
    hist_path = WEIGHT_DIR / "history.csv"
    if hist_path.exists():
        with open(hist_path) as f:
            for row in csv.DictReader(f):
                dt = datetime.strptime(row["Date/Time"], "%b %d %Y %I:%M:%S %p")
                weight = float(row["Weight (lb)"])
                date_str = dt.strftime("%Y-%m-%d")
                time_str = dt.strftime("%H:%M")
                by_day[date_str].append((dt, weight, time_str))

    # Take first reading per day (earliest timestamp)
    results = []
    for date_str in sorted(by_day):
        entries = sorted(by_day[date_str], key=lambda x: x[0])
        dt, weight, time_str = entries[0]
        results.append({
            "date": date_str,
            "weight_lbs": weight,
            "time": time_str,
        })

    out_path = WEIGHT_DIR / "weight.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["date", "weight_lbs", "time"])
        writer.writeheader()
        writer.writerows(results)

    total_readings = sum(len(v) for v in by_day.values())
    print(f"Entries: {total_readings} readings -> {len(results)} days")
    print(f"Range: {results[0]['date']} to {results[-1]['date']}")
    print(f"Written: {out_path}")


if __name__ == "__main__":
    main()
