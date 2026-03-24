#!/usr/bin/env python3
"""Sanity check macro data at both daily and individual food item level.

Atwater factors: fat=9, carbs=4, protein=4.
Expected calories = fat*9 + carbs*4 + protein*4

Reports:
  1. Unique food items where macros don't add up, sorted by frequency
  2. Daily totals where reported vs computed calories diverge
"""

import csv
from collections import Counter
from pathlib import Path

INTAKE_DIR = Path(__file__).parent


def check_item(cal_s, fat_s, carbs_s, prot_s):
    """Return (cal, expected, diff, pct) or None if not checkable."""
    try:
        cal = int(cal_s) if cal_s else None
        fat = int(fat_s) if fat_s else 0
        carbs = int(carbs_s) if carbs_s else 0
        prot = int(prot_s) if prot_s else 0
    except ValueError:
        return None

    expected = fat * 9 + carbs * 4 + prot * 4
    if expected == 0 or cal is None:
        return None

    diff = cal - expected
    pct = (diff / cal * 100) if cal else 0
    return cal, expected, diff, pct


def main():
    # --- Individual food items ---
    foods_path = INTAKE_DIR / "intake_foods.csv"
    with open(foods_path) as f:
        foods = list(csv.DictReader(f))

    # Collect all outlier entries grouped by food name
    outlier_counts = Counter()     # food name -> count
    outlier_examples = {}          # food name -> (cal, expected, diff, pct) from first seen

    for r in foods:
        if r["meal"] == "TOTAL":
            continue
        result = check_item(r["calories"], r["fat_g"], r["carbs_g"], r["protein_g"])
        if result is None:
            continue
        cal, expected, diff, pct = result
        if abs(pct) > 50 and abs(diff) > 50:
            name = r["food"]
            outlier_counts[name] += 1
            if name not in outlier_examples:
                outlier_examples[name] = (cal, expected, diff, pct)

    total_entries = sum(outlier_counts.values())

    print(f"=== FOOD ITEMS WITH MACRO/CALORIE MISMATCH ===")
    print(f"(>50% and >50 cal off from fat*9 + carbs*4 + protein*4)")
    print(f"Total entries: {total_entries} of {len(foods)} ({total_entries/len(foods)*100:.1f}%)")
    print(f"Unique foods: {len(outlier_counts)}")
    print()
    print(f"{'count':>5s}  {'cal':>5s} {'expct':>5s} {'diff':>6s} {'pct':>6s}  food")
    for name, count in outlier_counts.most_common():
        cal, expected, diff, pct = outlier_examples[name]
        print(f"{count:>5d}  {cal:>5d} {expected:>5d} {diff:>+6d} {pct:>+5.0f}%  {name}")
    print()

    # --- Daily totals ---
    daily_path = INTAKE_DIR / "intake_daily.csv"
    with open(daily_path) as f:
        daily = list(csv.DictReader(f))

    daily_outliers = []
    total_off = 0
    total_days = 0

    for r in daily:
        result = check_item(r["calories"], r["fat_g"], r["carbs_g"], r["protein_g"])
        if result is None:
            continue
        cal, expected, diff, pct = result
        total_off += diff
        total_days += 1
        if abs(pct) > 20 and abs(diff) > 100:
            daily_outliers.append((r["date"], cal, expected, diff, pct))

    daily_outliers.sort(key=lambda x: -abs(x[4]))
    avg_off = total_off / total_days if total_days else 0

    print(f"=== DAILY TOTALS ({total_days} days) ===")
    print(f"Average daily surplus (reported - computed): {avg_off:+.0f} cal")
    print(f"  (positive = reported > fat*9+carbs*4+protein*4)")
    print(f"  (gap is typically alcohol, sugar alcohols, or rounding)")
    print()
    print(f"Outliers (>20% and >100 cal off): {len(daily_outliers)}")
    if daily_outliers:
        print(f"{'date':>12s}  {'reported':>8s}  {'computed':>8s}  {'diff':>6s}  {'pct':>6s}")
        for date, cal, exp, diff, pct in daily_outliers[:30]:
            print(f"{date:>12s}  {cal:>8d}  {exp:>8d}  {diff:>+6d}  {pct:>+5.1f}%")
        if len(daily_outliers) > 30:
            print(f"  ... and {len(daily_outliers) - 30} more")


if __name__ == "__main__":
    main()
