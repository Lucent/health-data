#!/usr/bin/env python3
"""Fetch user's MFP food library (/food/mine) and write 1NF CSVs.

The /food/mine HTML page is a Next.js SPA whose `__NEXT_DATA__` script tag
embeds a `myFoods` query result. That payload contains every food in the
user's personal library — created or edited from public foods — with the
full nutritional panel (richer than /v2/foods/{id}: includes added_sugars,
sugar_alcohols, vitamin_d) and all serving sizes.

Outputs:
  intake/mfp_my_foods.csv          — PK mfp_id, one row per food + macros
  intake/mfp_my_foods_servings.csv — composite PK (mfp_id, serving_index)

These two tables LEFT JOIN to fully reconstruct any food. They are 1NF
(atomic columns, no nesting) so they can be hand-edited and re-used by
downstream processing.

Idempotent: re-running replaces the CSVs entirely from a fresh fetch.

Usage:
  python intake/mfp_mine.py            # fetch + write CSVs
  python intake/mfp_mine.py --offline  # parse existing /tmp dump (no request)
"""
import argparse
import csv
import datetime as dt
import json
import re
import sys
from pathlib import Path

import browser_cookie3
import requests

UA = "Mozilla/5.0 (X11; Linux x86_64; rv:149.0) Gecko/20100101 Firefox/149.0"
URL = "https://www.myfitnesspal.com/food/mine"

REPO = Path(__file__).parent
FOODS_PATH = REPO / "mfp_my_foods.csv"
SERVINGS_PATH = REPO / "mfp_my_foods_servings.csv"
CACHE_PATH = Path("/tmp/mfp_mine.html")


def paths_for_label(label):
    """If label given, save to mfp_my_foods_<label>.csv etc."""
    if not label:
        return FOODS_PATH, SERVINGS_PATH
    return (REPO / f"mfp_my_foods_{label}.csv",
            REPO / f"mfp_my_foods_{label}_servings.csv")

FOOD_COLS = [
    "mfp_id", "brand", "description",
    "public", "verified", "deleted",
    "country_code", "version", "user_id",
    "calories",
    "carbs_g", "net_carbs_g",
    "fat_g", "saturated_fat_g", "mono_fat_g", "poly_fat_g", "trans_fat_g",
    "protein_g", "cholest_mg", "sodium_mg", "potassium_mg",
    "sugar_g", "sugar_alcohols_g", "added_sugars_g", "fiber_g",
    "calcium_mg", "iron_mg", "vitamin_a", "vitamin_c", "vitamin_d",
    "grams_per_default_serving", "fdc_link", "fetched_at",
]
SERVING_COLS = [
    "mfp_id", "serving_index", "value", "unit", "gram_weight",
    "nutrition_multiplier",
]


def fetch_html():
    cj = browser_cookie3.firefox(domain_name="myfitnesspal.com")
    s = requests.Session()
    s.cookies = cj
    s.headers.update({"User-Agent": UA, "Accept": "text/html,*/*",
                      "Accept-Language": "en-US,en;q=0.5"})
    r = s.get(URL, timeout=20)
    r.raise_for_status()
    CACHE_PATH.write_text(r.text)
    print(f"  fetched {len(r.text)} bytes; cached {CACHE_PATH}", file=sys.stderr)
    return r.text


def parse_html(html):
    """Extract the `myFoods` query data from __NEXT_DATA__."""
    m = re.search(r'<script[^>]*id="__NEXT_DATA__"[^>]*>(.*?)</script>',
                  html, re.DOTALL)
    if not m:
        raise SystemExit("ERROR: __NEXT_DATA__ script tag not found in HTML")
    data = json.loads(m.group(1))
    queries = data["props"]["pageProps"]["dehydratedState"]["queries"]
    for q in queries:
        if q.get("queryKey") == ["myFoods"]:
            return q["state"]["data"]
    raise SystemExit("ERROR: 'myFoods' query not in dehydratedState")


def food_to_row(f, fetched_at):
    nut = f.get("nutritional_contents") or {}
    energy = nut.get("energy") or {}
    return {
        "mfp_id": f.get("id", ""),
        "brand": f.get("brand_name") or "",
        "description": f.get("description", ""),
        "public": f.get("public", ""),
        "verified": f.get("verified", ""),
        "deleted": f.get("deleted", ""),
        "country_code": f.get("country_code", ""),
        "version": f.get("version", ""),
        "user_id": f.get("user_id", ""),
        "calories": energy.get("value", ""),
        "carbs_g": nut.get("carbohydrates", ""),
        "net_carbs_g": nut.get("net_carbs", ""),
        "fat_g": nut.get("fat", ""),
        "saturated_fat_g": nut.get("saturated_fat", ""),
        "mono_fat_g": nut.get("monounsaturated_fat", ""),
        "poly_fat_g": nut.get("polyunsaturated_fat", ""),
        "trans_fat_g": nut.get("trans_fat", ""),
        "protein_g": nut.get("protein", ""),
        "cholest_mg": nut.get("cholesterol", ""),
        "sodium_mg": nut.get("sodium", ""),
        "potassium_mg": nut.get("potassium", ""),
        "sugar_g": nut.get("sugar", ""),
        "sugar_alcohols_g": nut.get("sugar_alcohols", ""),
        "added_sugars_g": nut.get("added_sugars", ""),
        "fiber_g": nut.get("fiber", ""),
        "calcium_mg": nut.get("calcium", ""),
        "iron_mg": nut.get("iron", ""),
        "vitamin_a": nut.get("vitamin_a", ""),
        "vitamin_c": nut.get("vitamin_c", ""),
        "vitamin_d": nut.get("vitamin_d", ""),
        "grams_per_default_serving": nut.get("grams", ""),
        "fdc_link": nut.get("nutrition_source_other", "") or "",
        "fetched_at": fetched_at,
    }


def serving_rows(f):
    out = []
    for sv in (f.get("serving_sizes") or []):
        out.append({
            "mfp_id": f.get("id", ""),
            "serving_index": sv.get("index", ""),
            "value": sv.get("value", ""),
            "unit": sv.get("unit", ""),
            "gram_weight": sv.get("gram_weight", ""),
            "nutrition_multiplier": sv.get("nutrition_multiplier", ""),
        })
    return out


def write_csvs(foods, foods_path=None, servings_path=None):
    foods_path = foods_path or FOODS_PATH
    servings_path = servings_path or SERVINGS_PATH
    fetched_at = dt.datetime.utcnow().isoformat(timespec="seconds") + "Z"
    food_rows = sorted((food_to_row(f, fetched_at) for f in foods),
                       key=lambda r: r["mfp_id"])
    s_rows = []
    for f in foods:
        s_rows.extend(serving_rows(f))
    s_rows.sort(key=lambda r: (r["mfp_id"], r["serving_index"]))

    with open(foods_path, "w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=FOOD_COLS)
        w.writeheader()
        w.writerows(food_rows)
    with open(servings_path, "w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=SERVING_COLS)
        w.writeheader()
        w.writerows(s_rows)
    print(f"  wrote {len(food_rows)} foods → {foods_path}", file=sys.stderr)
    print(f"  wrote {len(s_rows)} servings → {servings_path}", file=sys.stderr)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--offline", action="store_true",
                   help=f"parse {CACHE_PATH} instead of refetching")
    p.add_argument("--label", type=str, default=None,
                   help="save to mfp_my_foods_<label>.csv (default: no suffix). "
                        "Use when fetching /food/mine from a different logged-in "
                        "account — keeps libraries separate.")
    args = p.parse_args()

    if args.offline:
        if not CACHE_PATH.exists():
            print(f"ERROR: cache {CACHE_PATH} missing", file=sys.stderr)
            sys.exit(2)
        html = CACHE_PATH.read_text()
        print(f"  using cached {CACHE_PATH} ({len(html)} bytes)", file=sys.stderr)
    else:
        html = fetch_html()

    foods = parse_html(html)
    print(f"  parsed {len(foods)} foods from myFoods query", file=sys.stderr)

    # Report user_ids in the dataset so we know whose library this is
    user_ids = sorted({f.get("user_id") for f in foods if f.get("user_id")})
    print(f"  owning user_ids: {user_ids}", file=sys.stderr)

    foods_path, servings_path = paths_for_label(args.label)
    write_csvs(foods, foods_path, servings_path)


if __name__ == "__main__":
    main()
