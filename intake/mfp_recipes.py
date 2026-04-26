#!/usr/bin/env python3
"""Fetch user's MFP recipe box via /recipe_parser (HTML listing) + /v2/recipes
(JSON detail API). Writes 3 related 1NF CSV tables.

Two-phase pipeline, idempotent:
  Phase 1 (HTML pagination): enumerate every recipe_id for the logged-in user
  Phase 2 (v2 JSON API): for each recipe, GET /v2/recipes/<id> and extract
    - recipe-level macros (total batch — per-serving = total / servings)
    - ingredient list (each ingredient has food_id + full food object)

Outputs (per --label):
  intake/mfp_my_recipes_<label>.csv              — PK recipe_id
  intake/mfp_my_recipes_<label>_ingredients.csv  — composite PK (recipe_id, ingredient_index)
  intake/mfp_foods_from_recipes_<label>.csv      — deduped ingredient foods
                                                   (same schema as mfp_my_foods.csv)

The ingredient-food dictionary is a free bonus: one /v2/recipes call yields
the full food panel for every ingredient. These get merged into the
mine cache automatically (mfp_enrich.py globs all mfp_my_foods*.csv).

Usage:
  python intake/mfp_recipes.py --label caitlins              # phase 1 + phase 2
  python intake/mfp_recipes.py --label caitlins --listing-only
  python intake/mfp_recipes.py --label caitlins --details-only --limit 5
  python intake/mfp_recipes.py --label caitlins --force-refetch  # redo phase 2
"""
import argparse
import csv
import datetime as dt
import json
import random
import re
import sys
import time
from pathlib import Path

import browser_cookie3
import requests

UA = "Mozilla/5.0 (X11; Linux x86_64; rv:149.0) Gecko/20100101 Firefox/149.0"
LIST_URL = "https://www.myfitnesspal.com/recipe_parser"
RECIPE_API = "https://api.myfitnesspal.com/v2/recipes/{recipe_id}"
AUTH_URL = "https://www.myfitnesspal.com/user/auth_token?refresh=true"

RECIPE_FIELDS = [
    "nutritional_contents", "ingredients",
    "country_code", "user_id", "version",
]

REPO = Path(__file__).parent

RECIPE_COLS = [
    "recipe_id", "version", "name", "user_id", "public", "client_id",
    "servings", "prep_time_minutes", "country_code",
    "cooking_methods", "courses", "cuisines", "categories", "dietary_needs",
    "calories_total", "grams_total",
    "carbs_g", "net_carbs_g",
    "fat_g", "saturated_fat_g", "mono_fat_g", "poly_fat_g", "trans_fat_g",
    "protein_g", "cholest_mg", "sodium_mg", "potassium_mg",
    "sugar_g", "sugar_alcohols_g", "added_sugars_g", "fiber_g",
    "calcium_mg", "iron_mg", "vitamin_a", "vitamin_c", "vitamin_d",
    "fdc_link", "source_url", "fetched_at",
]

INGREDIENT_COLS = [
    "recipe_id", "ingredient_index",
    "raw_text",
    "ingredient_food_id", "ingredient_food_version",
    "ingredient_food_brand", "ingredient_food_description",
    "ingredient_food_verified",
    "serving_value", "serving_unit", "nutrition_multiplier",
    "ingredient_servings",
]

# Mirror the mfp_my_foods.csv schema for ingredient-food dump
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


def polite_sleep():
    time.sleep(random.uniform(1.0, 4.0))


def make_session():
    cj = browser_cookie3.firefox(domain_name="myfitnesspal.com")
    s = requests.Session()
    s.cookies = cj
    s.headers.update({"User-Agent": UA,
                      "Accept": "application/json, text/html, */*",
                      "Accept-Language": "en-US,en;q=0.5"})
    return s


def _user_id_from_token(access_token):
    """Access token base64-decodes to 'a:mfp-main-js:<user_id>::mfp-js:<ts>...'."""
    import base64
    try:
        decoded = base64.b64decode(access_token + "===").decode("utf-8", errors="replace")
        parts = decoded.split(":")
        # Layout: ['a', 'mfp-main-js', '<user_id>', '', 'mfp-js', '<ts>', ...]
        if len(parts) >= 3 and parts[2].isdigit():
            return parts[2]
    except Exception:
        pass
    return None


def refresh_auth(s, user_id=""):
    r = s.get(AUTH_URL, timeout=20)
    r.raise_for_status()
    ct = r.headers.get("Content-Type", "")
    if not ct.startswith("application/json"):
        raise RuntimeError("auth_token non-JSON response; session likely expired")
    auth = r.json()
    token = auth["access_token"]
    s.headers["Authorization"] = f"Bearer {token}"
    s.headers["mfp-client-id"] = "mfp-main-js"
    uid = user_id or _user_id_from_token(token)
    if not uid:
        raise RuntimeError("couldn't extract user_id from access_token")
    s.headers["mfp-user-id"] = str(uid)
    return uid


def fetch_listing_page(s, page):
    params = {} if page == 1 else {"page": page, "sort_order": "recent"}
    r = s.get(LIST_URL, params=params, timeout=20)
    r.raise_for_status()
    return r.text


def parse_listing(html):
    """Return list of {recipe_id, version_id, name, calories_per_serving, servings}."""
    out = []
    for m in re.finditer(r'<li class="recipe-info".*?</li>', html, re.DOTALL):
        block = m.group(0)
        rid_m = re.search(r'/recipe/view/(\d+)', block)
        version_m = re.search(r'/recipe/copy/\d+\?version=(\d+)', block)
        name_m = re.search(r'<a\s+title="([^"]+)"\s+href="/recipe/view/\d+"', block)
        cal_m = re.search(r'<i class="icon-fire"></i>\s*(\d+)', block)
        serv_m = re.search(r'icon-user"></i>\s*(\d+)\s*servings?', block)
        if not rid_m:
            continue
        out.append({
            "recipe_id": rid_m.group(1),
            "version": version_m.group(1) if version_m else "",
            "name": name_m.group(1).strip() if name_m else "",
            "calories_per_serving_listing": cal_m.group(1) if cal_m else "",
            "servings_listing": serv_m.group(1) if serv_m else "",
        })
    return out


def fetch_recipe_v2(s, recipe_id):
    qs = "&".join(f"fields%5B%5D={f}" for f in RECIPE_FIELDS)
    r = s.get(RECIPE_API.format(recipe_id=recipe_id) + "?" + qs, timeout=25)
    if not r.ok:
        return None, f"HTTP {r.status_code}: {r.text[:200]}"
    j = r.json()
    return j.get("item"), None


# --- row builders ---

def _safe(v):
    return "" if v is None else v


def nutrition_row(nut):
    nut = nut or {}
    energy = nut.get("energy") or {}
    return {
        "calories": _safe(energy.get("value")),
        "carbs_g": _safe(nut.get("carbohydrates")),
        "net_carbs_g": _safe(nut.get("net_carbs")),
        "fat_g": _safe(nut.get("fat")),
        "saturated_fat_g": _safe(nut.get("saturated_fat")),
        "mono_fat_g": _safe(nut.get("monounsaturated_fat")),
        "poly_fat_g": _safe(nut.get("polyunsaturated_fat")),
        "trans_fat_g": _safe(nut.get("trans_fat")),
        "protein_g": _safe(nut.get("protein")),
        "cholest_mg": _safe(nut.get("cholesterol")),
        "sodium_mg": _safe(nut.get("sodium")),
        "potassium_mg": _safe(nut.get("potassium")),
        "sugar_g": _safe(nut.get("sugar")),
        "sugar_alcohols_g": _safe(nut.get("sugar_alcohols")),
        "added_sugars_g": _safe(nut.get("added_sugars")),
        "fiber_g": _safe(nut.get("fiber")),
        "calcium_mg": _safe(nut.get("calcium")),
        "iron_mg": _safe(nut.get("iron")),
        "vitamin_a": _safe(nut.get("vitamin_a")),
        "vitamin_c": _safe(nut.get("vitamin_c")),
        "vitamin_d": _safe(nut.get("vitamin_d")),
    }


def recipe_to_row(item, listing_info, fetched_at):
    """Build row for mfp_my_recipes_<label>.csv."""
    nut = item.get("nutritional_contents") or {}
    row = {
        "recipe_id": item.get("id", ""),
        "version": item.get("version", ""),
        "name": item.get("name", listing_info.get("name", "")),
        "user_id": item.get("user_id", ""),
        "public": item.get("public", ""),
        "client_id": item.get("client_id", ""),
        "servings": item.get("servings", ""),
        "prep_time_minutes": item.get("prep_time_minutes", ""),
        "country_code": item.get("country_code", ""),
        "cooking_methods": ";".join(item.get("cooking_methods") or []),
        "courses": ";".join(item.get("courses") or []),
        "cuisines": ";".join(item.get("cuisines") or []),
        "categories": ";".join(item.get("categories") or []),
        "dietary_needs": ";".join(item.get("dietary_needs") or []),
        "grams_total": _safe(nut.get("grams")),
        "fdc_link": _safe(nut.get("nutrition_source_other")),
        "source_url": f"https://www.myfitnesspal.com/recipe/view/{item.get('id','')}",
        "fetched_at": fetched_at,
    }
    nr = nutrition_row(nut)
    # recipe nutrition is total-for-recipe; rename calories → calories_total
    nr["calories_total"] = nr.pop("calories")
    row.update(nr)
    return row


def ingredients_rows(item):
    """Build rows for mfp_my_recipes_<label>_ingredients.csv."""
    out = []
    for idx, ing in enumerate(item.get("ingredients") or []):
        food = ing.get("food") or {}
        sv = ing.get("serving_size") or {}
        out.append({
            "recipe_id": item.get("id", ""),
            "ingredient_index": idx,
            "raw_text": ing.get("raw_text", ""),
            "ingredient_food_id": food.get("id", ""),
            "ingredient_food_version": food.get("version", ""),
            "ingredient_food_brand": food.get("brand_name") or "",
            "ingredient_food_description": food.get("description", ""),
            "ingredient_food_verified": food.get("verified", ""),
            "serving_value": sv.get("value", ""),
            "serving_unit": sv.get("unit", ""),
            "nutrition_multiplier": sv.get("nutrition_multiplier", ""),
            "ingredient_servings": ing.get("servings", ""),
        })
    return out


def ingredient_food_to_row(food, fetched_at):
    """Convert ingredient's food to a mfp_my_foods-shaped row."""
    nut = food.get("nutritional_contents") or {}
    row = {
        "mfp_id": food.get("id", ""),
        "brand": food.get("brand_name") or "",
        "description": food.get("description", ""),
        "public": food.get("public", ""),
        "verified": food.get("verified", ""),
        "deleted": food.get("deleted", ""),
        "country_code": food.get("country_code", ""),
        "version": food.get("version", ""),
        "user_id": food.get("user_id", ""),
        "grams_per_default_serving": _safe(nut.get("grams")),
        "fdc_link": _safe(nut.get("nutrition_source_other")),
        "fetched_at": fetched_at,
    }
    row.update(nutrition_row(nut))  # adds calories, carbs_g, etc.
    return row


# --- persistence ---

def load_csv(path, key):
    if not path.exists():
        return {}
    with open(path) as f:
        return {r[key]: r for r in csv.DictReader(f)}


def save_csv(path, columns, rows_by_key):
    tmp = path.with_suffix(".csv.tmp")
    with open(tmp, "w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=columns)
        w.writeheader()
        for k in sorted(rows_by_key):
            w.writerow({c: rows_by_key[k].get(c, "") for c in columns})
    tmp.replace(path)


def save_ingredients_csv(path, rows):
    """Ingredients have composite key (recipe_id, ingredient_index)."""
    tmp = path.with_suffix(".csv.tmp")
    with open(tmp, "w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=INGREDIENT_COLS)
        w.writeheader()
        rows.sort(key=lambda r: (r["recipe_id"], int(r.get("ingredient_index") or 0)))
        for r in rows:
            w.writerow({c: r.get(c, "") for c in INGREDIENT_COLS})
    tmp.replace(path)


# --- main ---

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--label", type=str, required=True,
                   help="write to mfp_my_recipes_<label>.csv etc.")
    p.add_argument("--max-pages", type=int, default=None)
    p.add_argument("--listing-only", action="store_true")
    p.add_argument("--details-only", action="store_true")
    p.add_argument("--limit", type=int, default=None,
                   help="phase 2: cap number of recipe detail fetches (testing)")
    p.add_argument("--force-refetch", action="store_true",
                   help="redo phase 2 for rows that already have macros")
    args = p.parse_args()

    recipes_path = REPO / f"mfp_my_recipes_{args.label}.csv"
    ingredients_path = REPO / f"mfp_my_recipes_{args.label}_ingredients.csv"
    foods_path = REPO / f"mfp_foods_from_recipes_{args.label}.csv"

    s = make_session()
    recipes = load_csv(recipes_path, "recipe_id")
    ingredients_existing = {}
    if ingredients_path.exists():
        with open(ingredients_path) as f:
            for r in csv.DictReader(f):
                ingredients_existing[(r["recipe_id"], r["ingredient_index"])] = r
    foods = load_csv(foods_path, "mfp_id")

    print(f"  existing: {len(recipes)} recipes, "
          f"{len(ingredients_existing)} ingredient rows, "
          f"{len(foods)} ingredient foods", file=sys.stderr)

    fetched_at = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")

    # Phase 1: HTML pagination
    if not args.details_only:
        page = 1
        new_count = 0
        while True:
            if args.max_pages and page > args.max_pages:
                break
            print(f"  listing page {page}...", file=sys.stderr)
            try:
                html = fetch_listing_page(s, page)
            except Exception as e:
                print(f"    ERROR: {e}", file=sys.stderr)
                break
            batch = parse_listing(html)
            if not batch:
                print(f"    empty — stopping", file=sys.stderr)
                break
            new_here = 0
            for b in batch:
                rid = b["recipe_id"]
                if rid in recipes:
                    # update listing-derived fields only
                    recipes[rid]["name"] = b.get("name", recipes[rid].get("name", ""))
                    # preserve details if already fetched
                else:
                    recipes[rid] = {
                        "recipe_id": rid,
                        "version": b["version"],
                        "name": b["name"],
                        "fetched_at": fetched_at,
                    }
                    new_here += 1
            new_count += new_here
            print(f"    {len(batch)} shown ({new_here} new)", file=sys.stderr)
            save_csv(recipes_path, RECIPE_COLS, recipes)
            if f"page={page + 1}" not in html:
                print(f"    no next-page link — stopping", file=sys.stderr)
                break
            page += 1
            polite_sleep()
        print(f"\n  phase 1 done: {len(recipes)} total, {new_count} new", file=sys.stderr)

    if args.listing_only:
        return

    # Phase 2: v2 API details
    uid = refresh_auth(s)
    print(f"  auth refreshed (user_id={uid})", file=sys.stderr)

    needs_details = []
    for rid, r in recipes.items():
        if args.force_refetch or not r.get("calories_total"):
            needs_details.append(rid)
    if args.limit:
        needs_details = needs_details[:args.limit]
    print(f"\n  phase 2: {len(needs_details)} recipes need detail fetch", file=sys.stderr)

    for i, rid in enumerate(needs_details):
        polite_sleep()
        listing = recipes.get(rid, {})
        try:
            item, err = fetch_recipe_v2(s, rid)
        except Exception as e:
            print(f"  [{i+1}] {rid}: exception {type(e).__name__}: {e}", file=sys.stderr)
            continue
        if err:
            print(f"  [{i+1}] {rid}: {err}", file=sys.stderr)
            continue

        now_ts = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
        row = recipe_to_row(item, listing, now_ts)
        recipes[rid] = row

        # Ingredients rows
        ing_rows = ingredients_rows(item)
        # Replace existing ingredient rows for this recipe
        keys_to_drop = [k for k in ingredients_existing if k[0] == rid]
        for k in keys_to_drop:
            del ingredients_existing[k]
        for ir in ing_rows:
            ingredients_existing[(rid, str(ir["ingredient_index"]))] = ir

        # Ingredient foods dictionary (dedup by mfp_id; keep most recent fetched_at)
        for ing in (item.get("ingredients") or []):
            food = ing.get("food")
            if not food or not food.get("id"):
                continue
            frow = ingredient_food_to_row(food, now_ts)
            fid = frow["mfp_id"]
            prev = foods.get(fid)
            if prev is None or (prev.get("fetched_at", "") < frow["fetched_at"]):
                foods[fid] = frow

        # Persist after every recipe
        save_csv(recipes_path, RECIPE_COLS, recipes)
        save_ingredients_csv(ingredients_path, list(ingredients_existing.values()))
        save_csv(foods_path, FOOD_COLS, foods)

        name = item.get("name", "?")
        cal = (row.get("calories_total") or 0) and (float(row["calories_total"]) /
              (float(row.get("servings") or 1) or 1))
        ing_ct = len(ing_rows)
        print(f"  [{i+1}/{len(needs_details)}] {rid} {name!r}  "
              f"{cal:.0f}cal/serv  {ing_ct} ingredients", file=sys.stderr)


if __name__ == "__main__":
    main()
