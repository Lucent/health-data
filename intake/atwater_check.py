#!/usr/bin/env python3
"""Validate calories against macros using Atwater factors (fat=9, carbs=4, protein=4).

Reports:
  1. Per-item: unique food items where macros don't match calories, grouped by
     food name with frequency counts. Excludes alcohol and sugar-alcohol products
     (which have calories from sources Atwater doesn't cover).
  2. Per-day: daily totals where reported vs computed calories diverge.
"""

import csv
import re
from collections import Counter
from pathlib import Path

INTAKE_DIR = Path(__file__).parent

ALCOHOL = re.compile(
    r"(?i)beer\b|ale\b|wine\b|lager\b|ipa\b|stout\b|porter\b|cocktail|negroni|"
    r"old fashion|martini|margarita|bourbon|whiskey|whisky|vodka|rum\b|gin\b|"
    r"tequila|sour\b|mule\b|mimosa|champagne|prosecco|corona|michelob|guinness|"
    r"stella|newcastle|dogfish|sam adams|chimay|abita|miller.*lite|miller.*high|"
    r"modelo|cider|sauvignon|cabernet|pinot|merlot|chardonnay|sangria|shandy|"
    r"pilsner|kolsch|bier\b|bira\b|amstel|cerveza|spritz|mojito|brut\b|lager|"
    r"saison|kölsch|flight\b|harder|requiem|gaffel|funkwerks|redhook|surly|"
    r"kross|cigar cocktail|generi.*cabernet|red wine|white wine|pinot grigio|"
    r"la tulipe|charles shaw|oyster bay|boulevard|tulipe|yuengling|allagash|"
    r"oskar blues|death by coconut|fat tire|dos equis|leinenkugel|brad.*silver"
)

SUGAR_ALCOHOL = re.compile(
    r"(?i)chocorite|choczero|enlightened|quest|mission.*carb balance|better bagel|"
    r"lily|atkins|fiber one|smart swap|no cow|built bar|erythritol|sugar alcohol|"
    r"minus sugar|minus erythritol|net carb|keto|betterbar|betterbrand"
)


def safe_int(v):
    try:
        return int(v) if v else 0
    except ValueError:
        return 0


def atwater(fat, carbs, prot):
    return fat * 9 + carbs * 4 + prot * 4


def main():
    with open(INTAKE_DIR / "intake_foods.csv") as f:
        foods = list(csv.DictReader(f))
    with open(INTAKE_DIR / "intake_daily.csv") as f:
        daily = list(csv.DictReader(f))

    # --- Per-item check ---
    outlier_counts = Counter()
    outlier_examples = {}
    skipped_alcohol = 0
    total_items = 0

    for r in foods:
        if r["meal"] == "TOTAL":
            continue
        total_items += 1
        cal = safe_int(r["calories"])
        fat = safe_int(r["fat_g"])
        carbs = safe_int(r["carbs_g"])
        prot = safe_int(r["protein_g"])
        if cal == 0:
            continue
        exp = atwater(fat, carbs, prot)
        if exp == 0:
            continue

        diff = cal - exp
        if ALCOHOL.search(r["food"]) or SUGAR_ALCOHOL.search(r["food"]):
            skipped_alcohol += 1
            continue
        if abs(diff) > 100:
            name = r["food"]
            outlier_counts[name] += 1
            if name not in outlier_examples:
                pct = diff / cal * 100
                outlier_examples[name] = (cal, exp, diff, pct)

    total_outlier_entries = sum(outlier_counts.values())
    print(f"=== PER-ITEM ATWATER CHECK ===")
    print(f"Items checked: {total_items} ({skipped_alcohol} alcohol/sugar-alcohol excluded)")
    print(f"Outliers (>100 cal off): {total_outlier_entries} entries, {len(outlier_counts)} unique foods")
    print(f"  ({total_outlier_entries/total_items*100:.1f}% of items)")
    print()
    if outlier_counts:
        print(f"{'count':>5s}  {'cal':>5s} {'expct':>5s} {'diff':>6s} {'pct':>6s}  food")
        for name, count in outlier_counts.most_common():
            cal, exp, diff, pct = outlier_examples[name]
            print(f"{count:>5d}  {cal:>5d} {exp:>5d} {diff:>+6d} {pct:>+5.0f}%  {name}")
        print()

    # --- Per-day check ---
    daily_outliers = []
    total_off = 0
    total_days = 0

    for r in daily:
        cal = safe_int(r["calories"])
        fat = safe_int(r["fat_g"])
        carbs = safe_int(r["carbs_g"])
        prot = safe_int(r["protein_g"])
        exp = atwater(fat, carbs, prot)
        if exp == 0 or cal == 0:
            continue
        diff = cal - exp
        pct = diff / cal * 100
        total_off += diff
        total_days += 1
        if abs(pct) > 20 and abs(diff) > 100:
            daily_outliers.append((r["date"], cal, exp, diff, pct))

    daily_outliers.sort(key=lambda x: -abs(x[4]))
    avg_off = total_off / total_days if total_days else 0

    print(f"=== PER-DAY ATWATER CHECK ({total_days} days) ===")
    print(f"Average daily surplus (reported - computed): {avg_off:+.0f} cal")
    print(f"  (positive = reported > fat*9+carbs*4+protein*4)")
    print(f"  (gap is typically alcohol, sugar alcohols, or rounding)")
    print()
    print(f"Outliers (>20% and >100 cal off): {len(daily_outliers)}")
    if daily_outliers:
        print(f"{'date':>12s}  {'reported':>8s}  {'computed':>8s}  {'diff':>6s}  {'pct':>6s}")
        for dt, cal, exp, diff, pct in daily_outliers[:30]:
            print(f"{dt:>12s}  {cal:>8d}  {exp:>8d}  {diff:>+6d}  {pct:>+5.1f}%")
        if len(daily_outliers) > 30:
            print(f"  ... and {len(daily_outliers) - 30} more")


if __name__ == "__main__":
    main()
