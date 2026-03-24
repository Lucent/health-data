#!/usr/bin/env python3
"""Verify intake data integrity. Run after every merge.

Checks:
  1. Every date's item calorie sum matches its daily total (±1)
  2. No date gaps in the full range
  3. Monthly checksums match (regenerates checksums.csv)
  4. Food macro errors >100 cal (excluding alcohol and sugar alcohol products)
  5. No duplicate items (same date+meal+food+calories)

Exit code 0 = all clean, 1 = failures found.
"""

import csv
import re
import sys
from collections import Counter, defaultdict
from datetime import date, timedelta
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


def load():
    with open(INTAKE_DIR / "intake_foods.csv") as f:
        foods = list(csv.DictReader(f))
    with open(INTAKE_DIR / "intake_daily.csv") as f:
        daily = {r["date"]: r for r in csv.DictReader(f)}
    return foods, daily


def safe_int(v):
    try:
        return int(v)
    except (ValueError, TypeError):
        return 0


def check_item_vs_total(foods, daily):
    """Every date's item calorie sum must match daily total within ±1."""
    items_by_date = defaultdict(list)
    for r in foods:
        items_by_date[r["date"]].append(r)

    failures = []
    for d in sorted(items_by_date):
        item_cal = sum(safe_int(r["calories"]) for r in items_by_date[d])
        daily_cal = safe_int(daily[d]["calories"]) if d in daily else None
        if daily_cal is None:
            failures.append((d, item_cal, "MISSING", "no daily total"))
        elif abs(item_cal - daily_cal) > 1:
            failures.append((d, item_cal, daily_cal, item_cal - daily_cal))
    return failures


def check_continuity(daily):
    """No gaps in the date range."""
    all_dates = sorted(daily.keys())
    start = date.fromisoformat(all_dates[0])
    end = date.fromisoformat(all_dates[-1])
    present = set(all_dates)
    missing = []
    d = start
    while d <= end:
        if d.isoformat() not in present:
            missing.append(d.isoformat())
        d += timedelta(days=1)
    return missing, all_dates[0], all_dates[-1], len(all_dates)


def check_checksums(foods, daily):
    """Regenerate monthly checksums and check all match."""
    food_by_month = defaultdict(lambda: {"items": 0, "cal": 0, "days": set()})
    for r in foods:
        m = r["date"][:7]
        food_by_month[m]["items"] += 1
        food_by_month[m]["days"].add(r["date"])
        food_by_month[m]["cal"] += safe_int(r["calories"])

    daily_by_month = defaultdict(lambda: {"days": 0, "cal": 0})
    for r in daily.values():
        m = r["date"][:7]
        daily_by_month[m]["days"] += 1
        daily_by_month[m]["cal"] += safe_int(r["calories"])

    months = sorted(set(food_by_month) | set(daily_by_month))
    failures = []
    rows = []
    for m in months:
        fi = food_by_month[m]["items"]
        fc = food_by_month[m]["cal"]
        dd = daily_by_month[m]["days"]
        dc = daily_by_month[m]["cal"]
        match = abs(fc - dc) <= dd
        rows.append((m, dd, fi, fc, dc, match))
        if not match:
            failures.append((m, fc, dc, fc - dc))

    # Write checksums
    with open(INTAKE_DIR / "checksums.csv", "w", newline="") as f:
        f.write("month,days,items,item_cal_sum,daily_cal_sum,checksums_match\n")
        for m, dd, fi, fc, dc, match in rows:
            f.write(f"{m},{dd},{fi},{fc},{dc},{'YES' if match else 'NO'}\n")

    return failures


def check_macros(foods):
    """Find food items where macros don't add up to calories (>100 cal off)."""
    outliers = []
    for r in foods:
        if r["meal"] == "TOTAL":
            continue
        if ALCOHOL.search(r["food"]) or SUGAR_ALCOHOL.search(r["food"]):
            continue
        cal = safe_int(r["calories"])
        fat = safe_int(r["fat_g"])
        carbs = safe_int(r["carbs_g"])
        prot = safe_int(r["protein_g"])
        if cal == 0:
            continue
        exp = fat * 9 + carbs * 4 + prot * 4
        if exp == 0:
            continue
        diff = cal - exp
        if abs(diff) > 100:
            outliers.append((r["date"], r["meal"], r["food"], cal, exp, diff))
    return outliers


def check_duplicates(foods):
    """Find exact duplicate items (same date+meal+food+calories)."""
    counts = Counter(
        (r["date"], r["meal"], r["food"], r["calories"])
        for r in foods if r["meal"] != "TOTAL"
    )
    return [(k, v) for k, v in counts.items() if v > 1]


def main():
    foods, daily = load()
    ok = True

    # 1. Item sum vs daily total
    failures = check_item_vs_total(foods, daily)
    if failures:
        ok = False
        print(f"FAIL: {len(failures)} dates with item sum != daily total")
        for d, isum, dtot, diff in failures[:10]:
            print(f"  {d}: items={isum} daily={dtot} diff={diff}")
    else:
        print(f"OK: All {len(daily)} dates: item sums match daily totals")

    # 2. Continuity
    missing, first, last, count = check_continuity(daily)
    if missing:
        ok = False
        print(f"FAIL: {len(missing)} missing dates in {first} to {last}")
        for d in missing[:10]:
            print(f"  {d}")
    else:
        print(f"OK: {count} consecutive days, {first} to {last}, zero gaps")

    # 3. Monthly checksums
    failures = check_checksums(foods, daily)
    if failures:
        ok = False
        print(f"FAIL: {len(failures)} months with checksum mismatch")
        for m, fc, dc, diff in failures:
            print(f"  {m}: item_sum={fc} daily_sum={dc} diff={diff:+d}")
    else:
        months = len(set(r["date"][:7] for r in daily.values()))
        print(f"OK: All {months} months checksum (written to checksums.csv)")

    # 4. Macro errors
    outliers = check_macros(foods)
    unique = len(set(o[2] for o in outliers))
    print(f"INFO: {len(outliers)} macro errors >100 cal ({unique} unique foods)")
    # Show post-2018 errors >200 cal (fixable)
    fixable = [o for o in outliers if o[0] >= "2018-01" and abs(o[5]) > 200]
    if fixable:
        print(f"  Fixable (post-2018, >200 cal):")
        for d, meal, food, cal, exp, diff in sorted(fixable, key=lambda x: -abs(x[5])):
            print(f"    {d} {meal}: cal={cal} exp={exp} diff={diff:+d}  {food[:50]}")

    # 5. Duplicates
    dupes = check_duplicates(foods)
    if dupes:
        print(f"WARN: {len(dupes)} duplicate items")
        for (d, meal, food, cal), count in dupes[:5]:
            print(f"  {d} {meal}: {food[:40]} (cal={cal}) x{count}")
    else:
        print(f"OK: No duplicate items")

    # Summary
    print()
    print(f"{'='*50}")
    print(f"Items: {len(foods)}  Days: {len(daily)}  {'ALL CLEAN' if ok else 'FAILURES FOUND'}")
    print(f"{'='*50}")

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
