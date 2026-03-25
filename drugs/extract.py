"""Extract medicine log from XLSX to CSV.

Outputs two files:
- drugs/medicine.csv: all entries (tirzepatide, vaccines, metformin)
- drugs/tirzepatide.csv: just tirzepatide with dose, day-of-cycle, and
  interpolated daily dose level for joining with daily tables

Idempotent. Re-run to regenerate.
"""

import csv
from pathlib import Path
from datetime import timedelta
import openpyxl

ROOT = Path(__file__).resolve().parent
XLSX = ROOT / "medicine.xlsx"
OUT_ALL = ROOT / "medicine.csv"
OUT_TIRZ = ROOT / "tirzepatide.csv"

import math

def pk_curve(dose_mg, t_days):
    """One-compartment SC pharmacokinetic curve.

    Returns relative blood concentration at t_days after injection of dose_mg.
    Units are arbitrary (proportional to serum concentration in mg).
    """
    if t_days <= 0:
        return 0.0
    return (dose_mg * TIRZ_KA / (TIRZ_KA - TIRZ_KE)
            * (math.exp(-TIRZ_KE * t_days) - math.exp(-TIRZ_KA * t_days)))


ALL_FIELDS = ["date", "drug", "dose_mg", "location", "subjective_strength"]
TIRZ_DAILY_FIELDS = ["date", "dose_mg", "days_since_injection", "injection_date",
                     "blood_level", "effective_level"]

# Tachyphylaxis: dose effectiveness decays exponentially with time on current dose.
# Fitted from 530 days of intake data: effective = blood_level * exp(-decay * weeks_on_dose)
# Correlation with daily intake improves from r=-0.40 (raw) to r=-0.43 (with adaptation)
TIRZ_DECAY_RATE = 0.0217  # per week. Half-life of effectiveness: 32 weeks.

# Tirzepatide pharmacokinetics (FDA prescribing information)
# t1/2 = 5.0 days, Tmax = 24h median (SC injection)
# ka solved from Tmax = ln(ka/ke)/(ka-ke) = 1.0 day
TIRZ_HALF_LIFE_DAYS = 5.0
TIRZ_KE = math.log(2) / TIRZ_HALF_LIFE_DAYS  # 0.1386/day
TIRZ_KA = 3.3122  # absorption rate, gives Tmax = 1.0 day
# Steady state reached ~week 5, accumulation ~1.6x (FDA: 1.4x Cmax, 2x AUC)


def extract():
    wb = openpyxl.load_workbook(XLSX, data_only=True)
    ws = wb.active

    rows = list(ws.iter_rows(min_row=2, values_only=True))
    print(f"Read {len(rows)} rows from {XLSX.name}")

    all_records = []
    tirz_injections = []

    for row in rows:
        date_val, drug, dose, location = row[0], row[1], row[2], row[3]
        subjective = row[7] if len(row) > 7 else None

        if date_val is None or drug is None:
            continue

        date_str = date_val.strftime("%Y-%m-%d") if hasattr(date_val, "strftime") else str(date_val)[:10]

        rec = {
            "date": date_str,
            "drug": str(drug).strip(),
            "dose_mg": float(dose) if dose is not None else "",
            "location": str(location).strip() if location else "",
            "subjective_strength": str(subjective).strip() if subjective and str(subjective) != "nan" else "",
        }
        all_records.append(rec)

        if "tirzepatide" in str(drug).lower() or "zepbound" in str(drug).lower():
            tirz_injections.append({
                "date": date_val if hasattr(date_val, "strftime") else None,
                "date_str": date_str,
                "dose_mg": float(dose) if dose is not None else 0,
            })

    # Write all medicines
    all_records.sort(key=lambda r: r["date"])
    with open(OUT_ALL, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ALL_FIELDS)
        w.writeheader()
        w.writerows(all_records)
    print(f"Wrote {len(all_records)} entries to {OUT_ALL.name}")

    # Build daily tirzepatide table with pharmacokinetic blood level
    if tirz_injections:
        tirz_injections.sort(key=lambda r: r["date_str"])
        first = tirz_injections[0]["date"]
        last = tirz_injections[-1]["date"]

        daily_tirz = []
        current_dose = 0
        current_injection_date = None
        inj_idx_track = 0

        day = first
        while day <= last + timedelta(days=6):
            date_str = day.strftime("%Y-%m-%d")

            # Track most recent injection for dose/days_since columns
            while (inj_idx_track < len(tirz_injections)
                   and tirz_injections[inj_idx_track]["date_str"] <= date_str):
                current_dose = tirz_injections[inj_idx_track]["dose_mg"]
                current_injection_date = tirz_injections[inj_idx_track]["date"]
                inj_idx_track += 1

            # Compute blood level: sum of PK curves from ALL past injections
            blood_level = 0.0
            for inj in tirz_injections:
                if inj["date"] is None:
                    continue
                t = (day - inj["date"]).total_seconds() / 86400  # days since injection
                if t < 0:
                    break  # future injections
                if t > 30:
                    continue  # negligible after ~6 half-lives
                blood_level += pk_curve(inj["dose_mg"], t)

            if current_injection_date:
                days_since = (day - current_injection_date).days
                # Tachyphylaxis: weeks on current dose
                current_dose_start = None
                for inj in tirz_injections:
                    if inj["dose_mg"] == current_dose and inj["date"] is not None:
                        current_dose_start = inj["date"]
                        break
                if current_dose_start:
                    weeks_on_dose = (day - current_dose_start).days / 7.0
                else:
                    weeks_on_dose = 0
                effective = blood_level * math.exp(-TIRZ_DECAY_RATE * weeks_on_dose)

                daily_tirz.append({
                    "date": date_str,
                    "dose_mg": current_dose,
                    "days_since_injection": days_since,
                    "injection_date": current_injection_date.strftime("%Y-%m-%d"),
                    "blood_level": round(blood_level, 2),
                    "effective_level": round(effective, 2),
                })

            day += timedelta(days=1)

        with open(OUT_TIRZ, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=TIRZ_DAILY_FIELDS)
            w.writeheader()
            w.writerows(daily_tirz)

        doses = sorted(set(t["dose_mg"] for t in tirz_injections))
        print(f"Wrote {len(daily_tirz)} daily rows to {OUT_TIRZ.name}")
        print(f"  Injections: {len(tirz_injections)}")
        print(f"  Doses: {doses}")
        print(f"  Range: {tirz_injections[0]['date_str']} to {tirz_injections[-1]['date_str']}")

        # Show blood level at steady state for each dose
        print(f"\n  Steady-state blood levels (arbitrary units, proportional to serum concentration):")
        for dose in doses:
            # Simulate 8 weeks of weekly injection at this dose
            level = 0.0
            for week in range(8):
                t = 6.5  # measure at trough (day before next injection)
                t_since = t + week * 7
                level_at_t = sum(
                    pk_curve(dose, 6.5 + w * 7) for w in range(week + 1)
                )
            # Steady state trough ≈ level after 8 weeks
            ss_trough = sum(pk_curve(dose, 6.5 + w * 7) for w in range(8))
            ss_peak = sum(pk_curve(dose, 1.5 + w * 7) for w in range(8))
            print(f"    {dose:4.1f}mg weekly: trough={ss_trough:.1f}  peak={ss_peak:.1f}")


if __name__ == "__main__":
    extract()
