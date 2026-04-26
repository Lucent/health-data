#!/usr/bin/env python3
"""MFP food-dictionary enrichment — read-only, idempotent, 1NF.

Builds intake/mfp_food_dict.csv: one row per distinct `food` string in
intake_foods.csv with MFP-derived per-entry macros including saturated /
mono / poly / trans fat, potassium, calcium, iron, vitamins A/C/D, and
USDA FDC link when available.

Analysis pipeline LEFT JOINs on:
    intake/intake_foods.csv.food  ←→  intake/mfp_food_dict.csv.diary_food_name
Macros are stored per diary entry (already scaled).

Lookup order (per food):
  1. Mine cache (intake/mfp_my_foods.csv): exact (brand, description) match
     after annotation-stripping. Zero MFP requests. Calorie-ratio verified.
  2. MFP public search: build candidate list, score by name similarity.
     Try exact-match candidates first (tier 1), then approximate candidates
     ordered by (verified-first, score-desc) (tier 2 = USDA-verified +
     scalar-verified; tier 3 = scalar-verified only). Up to
     MAX_CANDIDATES_TO_VERIFY detail fetches per food.

Verification: scale = diary_calories / mfp_base_calories, then every other
diary macro must equal mfp_base_macro × scale within tolerance:
    |diff| <= max(1.0, 0.06 * expected)

Idempotent. Re-runs skip status='resolved*'. Re-process via:
  --force NAME      (one diary food)
  --status STATUS   (whole bucket)
  --mfp-id ID       (manual override: skip lookup entirely)

Auth: Firefox cookies via browser_cookie3, UA spoofed to Firefox 149.0.
NEVER mutates the diary. Read-only calls only.
"""
import argparse
import csv
import datetime as dt
import random
import re
import sys
import time
from collections import Counter
from pathlib import Path

import browser_cookie3
import cloudscraper
import myfitnesspal
import requests

UA = "Mozilla/5.0 (X11; Linux x86_64; rv:149.0) Gecko/20100101 Firefox/149.0"
USERNAME = "Lucent"

REPO = Path(__file__).parent
DICT_PATH = REPO / "mfp_food_dict.csv"
INTAKE_PATH = REPO / "intake_foods.csv"
MINE_FOODS_PATH = REPO / "mfp_my_foods.csv"
MINE_SERVINGS_PATH = REPO / "mfp_my_foods_servings.csv"

RATE_MIN, RATE_MAX = 1.0, 4.0
MAX_CANDIDATES_TO_VERIFY = 10

# Score thresholds for accepting an approximate match
APPROX_SCORE_THRESHOLD = 0.55  # min similarity to even consider as candidate

STATUS_RESOLVED = "resolved"                       # tier 1: exact match
STATUS_RESOLVED_APPROX_VERIFIED = "resolved_approx_verified"  # tier 2
STATUS_RESOLVED_APPROX = "resolved_approx"         # tier 3
STATUS_RESOLVED_FROM_MINE = "resolved_from_mine"   # mine-cache hit
STATUS_RESOLVED_OVERRIDE = "resolved_override"     # manual --mfp-id, trusted
STATUS_NOT_FOUND = "match_not_found"
STATUS_PARSE = "parse_fail"
STATUS_NO_CALS = "diary_zero_calories"

RESOLVED_STATUSES = {
    STATUS_RESOLVED, STATUS_RESOLVED_APPROX_VERIFIED,
    STATUS_RESOLVED_APPROX, STATUS_RESOLVED_FROM_MINE,
    STATUS_RESOLVED_OVERRIDE,
}

COLUMNS = [
    "diary_food_name", "status", "match_tier", "match_score",
    "mfp_id", "mfp_brand", "mfp_name", "mfp_verified", "source",
    "diary_qty", "diary_unit_string", "scale_factor",
    "atwater_drift",  # MFP-side Atwater consistency (downstream quality filter)
    "calories", "carbs_g", "fat_g",
    "saturated_fat_g", "mono_fat_g", "poly_fat_g", "trans_fat_g",
    "protein_g", "cholest_mg", "sodium_mg", "potassium_mg",
    "sugar_g", "sugar_alcohols_g", "added_sugars_g", "fiber_g",
    "calcium_mg", "iron_mg", "vitamin_a", "vitamin_c", "vitamin_d",
    "fdc_link", "notes", "fetched_at",
]


def polite_sleep():
    time.sleep(random.uniform(RATE_MIN, RATE_MAX))


# --- Diary food-name parser ---

FRAC_RE = re.compile(r"^(\d+)/(\d+)$")


def parse_qty(s):
    s = s.strip()
    m = FRAC_RE.match(s)
    if m:
        return float(m.group(1)) / float(m.group(2))
    try:
        return float(s)
    except ValueError:
        return None


def parse_diary_name(s):
    """Parse 'Brand - Name, qty unit_string' → dict.

    Brand and name keep user annotations like '(Corrected Per Bag)'.
    unit_string keeps everything after the qty.
    """
    out = {"brand": None, "name": None, "qty": None, "unit_string": None, "err": None}
    m = None
    for match in re.finditer(r",\s*([\d./]+)\s+(.+)$", s):
        m = match
    if not m:
        out["err"] = "no ', qty unit' suffix"
        return out
    qty = parse_qty(m.group(1))
    if qty is None:
        out["err"] = f"bad qty {m.group(1)!r}"
        return out
    food_part = s[:m.start()].strip()
    if " - " in food_part:
        brand, _, name = food_part.partition(" - ")
        brand, name = brand.strip(), name.strip()
    else:
        brand, name = None, food_part
    out.update(brand=brand, name=name, qty=qty,
               unit_string=m.group(2).strip())
    return out


# --- Annotation stripping + name similarity ---

ANNOTATION_RE = re.compile(
    r"\(\s*\*?\s*(?:Corrected|Minus\s+Erythritol|Minus\s+Sugar\s+Alcohol(?:s)?|"
    r"Without\s+Erythritol|New|Old|Updated|Fixed|Per\s+Bag|"
    r"Corrected\s+Per\s+Bag)\b[^)]*\*?\s*\)",
    re.IGNORECASE,
)
TRAILING_NOTE_RE = re.compile(r"\s+\d{1,2}/\d{1,2}/\d{2,4}\b")  # "4/25/13" date suffix


def strip_annotations(s):
    if not s:
        return ""
    s = ANNOTATION_RE.sub("", s)
    s = TRAILING_NOTE_RE.sub("", s)
    return re.sub(r"\s+", " ", s).strip()


WORD_RE = re.compile(r"[A-Za-z0-9%]+")


def tokens(s):
    s = (s or "").lower()
    # treat singular/plural endings as same
    raw = WORD_RE.findall(s)
    return {t[:-1] if t.endswith("s") and len(t) > 3 else t for t in raw}


def jaccard(a, b):
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def norm(s):
    return (s or "").strip().lower()


def score_match(diary_brand, diary_name, mfp_brand, mfp_name):
    """Score (0..1) how well an MFP entry matches the diary food.

    1.00 — exact (brand, name) string match
    0.85 — exact match after annotation-stripping
    0.70..0.85 — same brand, near-token-equal name
    0.40..0.70 — name token overlap, brand related
    0.00 — reject
    """
    db, dn = norm(diary_brand), norm(diary_name)
    mb, mn = norm(mfp_brand), norm(mfp_name)
    dn_strip = norm(strip_annotations(diary_name))
    mn_strip = norm(strip_annotations(mfp_name))

    # Brand alignment
    if db and mb:
        if db == mb:
            brand_aff = 1.0
        elif db in mb or mb in db:
            brand_aff = 0.85
        else:
            # different brands — but if names match, still possible
            brand_aff = 0.0
    elif not db and not mb:
        brand_aff = 1.0
    elif (db == dn) and not mb:
        # diary "Avocado - Avocado" with MFP empty brand
        brand_aff = 1.0
    elif not db and mb:
        # diary has no brand but MFP does; weak alignment
        brand_aff = 0.5
    else:
        brand_aff = 0.5

    # Name alignment
    if dn == mn:
        name_aff = 1.0
    elif dn_strip == mn_strip and dn_strip:
        name_aff = 0.85
    elif dn_strip == mn or dn == mn_strip:
        name_aff = 0.85
    else:
        td, tm = tokens(dn_strip), tokens(mn_strip)
        if not td or not tm:
            name_aff = 0.0
        elif tm.issubset(td) or td.issubset(tm):
            # Subset match — one is contained in the other. High-confidence
            # identity signal even when token counts differ (diary verbose,
            # MFP terse, or vice versa).
            smaller, larger = sorted((len(td), len(tm)))
            name_aff = 0.55 + 0.35 * (smaller / larger)  # 0.55..0.9
        else:
            sim = jaccard(td, tm)
            name_aff = sim if sim >= 0.5 else 0.0

    # Combined score; penalize zero brand_aff hard
    if brand_aff == 0.0:
        return 0.0
    if name_aff == 0.0:
        return 0.0
    if name_aff == 1.0 and brand_aff == 1.0:
        return 1.0
    return round(0.5 * brand_aff + 0.5 * name_aff, 3)


# --- MFP client ---

def make_client():
    cj = browser_cookie3.firefox(domain_name="myfitnesspal.com")

    def _session(sess=None, **kw):
        s = sess if sess is not None else requests.Session()
        s.headers.update({
            "User-Agent": UA,
            "Accept": "application/json, text/html, */*",
            "Accept-Language": "en-US,en;q=0.5",
        })
        return s
    cloudscraper.create_scraper = _session
    myfitnesspal.Client._get_user_metadata = lambda self: {"username": USERNAME}
    return myfitnesspal.Client(cookiejar=cj)


# --- Scalar (calorie-ratio) verification ---

# Structural macros: tightly coupled to calories, MFP rarely mis-stores these.
# If they don't scale-verify, it's not the same food. Required.
STRUCTURAL_MACROS = [
    ("carbs_g", "carbohydrates", 0.5),
    ("fat_g", "fat", 0.5),
    ("protein_g", "protein", 0.5),
]

# Incidental macros: MFP edits these over years (sodium rounded, cholesterol
# measured or not, fiber added later, sugar split from added_sugar later).
# Drift here is normal; tolerate up to MAX_INCIDENTAL_MISMATCH failures.
INCIDENTAL_MACROS = [
    ("cholest_mg", "cholesterol", 1.0),
    ("sodium_mg", "sodium", 5.0),
    ("sugars_g", "sugar", 1.0),
    ("fiber_g", "fiber", 0.5),
]
MAX_INCIDENTAL_MISMATCH = 2

# Legacy MACRO_MAP kept for back-compat
MACRO_MAP = STRUCTURAL_MACROS + INCIDENTAL_MACROS


def safe_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _check_macro(diary_row, raw_nut, scale, diary_field, mfp_key, tol_abs, tol_pct=0.06):
    """Returns (checked_bool, problem_tuple_or_None)."""
    dv = safe_float(diary_row.get(diary_field))
    if dv is None:
        return False, None
    mv = safe_float(raw_nut.get(mfp_key))
    if mv is None:
        if dv != 0:
            return True, (diary_field, dv, "mfp=None")
        return True, None
    expected = mv * scale
    allow = max(tol_abs, abs(expected) * tol_pct)
    # MFP rounds trace values to 0 — tolerate small diary nonzero against MFP zero
    if mv == 0 and abs(dv) <= 5:
        allow = max(allow, abs(dv))
    if abs(dv - expected) > allow:
        return True, (diary_field, dv, round(expected, 3))
    return True, None


def verify_by_calorie_ratio(diary_row, raw_nut):
    """Returns (scale, structural_checked, structural_problems, incidental_problems)."""
    base_cal = safe_float((raw_nut.get("energy") or {}).get("value"))
    diary_cal = safe_float(diary_row.get("calories"))
    if base_cal is None or base_cal == 0:
        return None, 0, [("calories", diary_cal, "mfp_base=0")], []
    if diary_cal is None:
        return None, 0, [("calories", "missing", base_cal)], []
    scale = diary_cal / base_cal

    s_problems = []
    s_checked = 1  # calories implicit
    for df, mk, tol in STRUCTURAL_MACROS:
        ck, prob = _check_macro(diary_row, raw_nut, scale, df, mk, tol)
        if ck:
            s_checked += 1
        if prob:
            s_problems.append(prob)

    i_problems = []
    for df, mk, tol in INCIDENTAL_MACROS:
        _, prob = _check_macro(diary_row, raw_nut, scale, df, mk, tol)
        if prob:
            i_problems.append(prob)

    return scale, s_checked, s_problems, i_problems


def verify_ok(s_checked, s_problems, i_problems,
              max_structural=0, max_incidental=MAX_INCIDENTAL_MISMATCH):
    """The food is a scale-verified identity if:
      - At least 2 structural fields checked
      - All structural macros match within tolerance (max_structural=0 default)
      - At most max_incidental incidental-macro mismatches (MFP DB drift OK)
    """
    return (s_checked >= 2
            and len(s_problems) <= max_structural
            and len(i_problems) <= max_incidental)


def scaled(v, scale):
    if v is None or scale is None:
        return ""
    return round(float(v) * scale, 3)


def build_row(food_name, parsed, source, status, tier, score,
              mfp_id, brand, name, verified, scale, raw, notes,
              diary_row=None):
    row = {c: "" for c in COLUMNS}
    row["diary_food_name"] = food_name
    row["fetched_at"] = dt.datetime.utcnow().isoformat(timespec="seconds") + "Z"
    row["diary_qty"] = parsed["qty"]
    row["diary_unit_string"] = parsed["unit_string"]
    row["status"] = status
    row["match_tier"] = tier
    row["match_score"] = score
    row["source"] = source
    row["notes"] = notes
    if mfp_id is not None:
        row["mfp_id"] = mfp_id
        row["mfp_brand"] = brand or ""
        row["mfp_name"] = name or ""
        row["mfp_verified"] = "" if verified is None else str(verified)

    # Per-macro drift protection: if a scaled MFP macro differs from the
    # diary's logged value by more than 60% relative, leave that column
    # blank. Downstream LEFT JOIN + COALESCE will then use diary's
    # authoritative value for that macro while still benefiting from the
    # other enrichment columns (sat/mono/poly/potassium/vitamins).
    def _blanked_or_scaled(scaled_val, diary_field):
        if scaled_val is None or scaled_val == "" or diary_row is None:
            return scaled_val
        dv = safe_float(diary_row.get(diary_field))
        if dv is None:
            return scaled_val
        try:
            sv = float(scaled_val)
        except (TypeError, ValueError):
            return scaled_val
        # Relative drift vs diary (diary is authoritative; 1.0 floor to avoid
        # div-by-zero on sparse 0-valued macros). 25% tolerance: beyond that,
        # the scaled enrichment value differs too much from what the user
        # actually logged; store blank so downstream COALESCE uses diary.
        denom = max(abs(dv), 1.0)
        if abs(sv - dv) / denom > 0.25:
            return ""  # fall back to diary
        return scaled_val

    if scale is not None and raw is not None:
        row["scale_factor"] = round(scale, 6)
        drift = atwater_drift(raw)
        row["atwater_drift"] = "" if drift is None else round(drift, 4)
        # Structural + incidental macros that overlap with diary get drift
        # protection. Micronutrients (no diary counterpart) scale unconditionally.
        row.update(
            calories=_blanked_or_scaled(scaled((raw.get("energy") or {}).get("value"), scale), "calories"),
            carbs_g=_blanked_or_scaled(scaled(raw.get("carbohydrates"), scale), "carbs_g"),
            fat_g=_blanked_or_scaled(scaled(raw.get("fat"), scale), "fat_g"),
            saturated_fat_g=scaled(raw.get("saturated_fat"), scale),
            mono_fat_g=scaled(raw.get("monounsaturated_fat"), scale),
            poly_fat_g=scaled(raw.get("polyunsaturated_fat"), scale),
            trans_fat_g=scaled(raw.get("trans_fat"), scale),
            protein_g=_blanked_or_scaled(scaled(raw.get("protein"), scale), "protein_g"),
            cholest_mg=_blanked_or_scaled(scaled(raw.get("cholesterol"), scale), "cholest_mg"),
            sodium_mg=_blanked_or_scaled(scaled(raw.get("sodium"), scale), "sodium_mg"),
            potassium_mg=scaled(raw.get("potassium"), scale),
            sugar_g=_blanked_or_scaled(scaled(raw.get("sugar"), scale), "sugars_g"),
            sugar_alcohols_g=scaled(raw.get("sugar_alcohols"), scale),
            added_sugars_g=scaled(raw.get("added_sugars"), scale),
            fiber_g=_blanked_or_scaled(scaled(raw.get("fiber"), scale), "fiber_g"),
            calcium_mg=scaled(raw.get("calcium"), scale),
            iron_mg=scaled(raw.get("iron"), scale),
            vitamin_a=scaled(raw.get("vitamin_a"), scale),
            vitamin_c=scaled(raw.get("vitamin_c"), scale),
            vitamin_d=scaled(raw.get("vitamin_d"), scale),
            fdc_link=raw.get("nutrition_source_other", "") or "",
        )
    return row


# --- Mine cache ---

MINE_FIELD_MAP = {
    # mine CSV column → raw_nut key (matches /v2/foods schema)
    "calories": ("energy", "value"),
    "carbs_g": "carbohydrates",
    "net_carbs_g": "net_carbs",
    "fat_g": "fat",
    "saturated_fat_g": "saturated_fat",
    "mono_fat_g": "monounsaturated_fat",
    "poly_fat_g": "polyunsaturated_fat",
    "trans_fat_g": "trans_fat",
    "protein_g": "protein",
    "cholest_mg": "cholesterol",
    "sodium_mg": "sodium",
    "potassium_mg": "potassium",
    "sugar_g": "sugar",
    "sugar_alcohols_g": "sugar_alcohols",
    "added_sugars_g": "added_sugars",
    "fiber_g": "fiber",
    "calcium_mg": "calcium",
    "iron_mg": "iron",
    "vitamin_a": "vitamin_a",
    "vitamin_c": "vitamin_c",
    "vitamin_d": "vitamin_d",
    "grams_per_default_serving": "grams",
}


def _recipes_as_foods(path):
    """Yield recipe rows converted into food-shaped dicts for mine cache.

    Recipes store total-batch macros; we divide by `servings` to get
    per-serving macros, matching the per-default-serving convention of
    mfp_my_foods. Uses recipe_id as mfp_id.
    """
    if not path.exists():
        return
    with open(path) as f:
        for r in csv.DictReader(f):
            rid = r.get("recipe_id")
            servings = safe_float(r.get("servings"))
            if not rid or not servings or servings <= 0:
                continue
            raw = {}
            # Recipe columns are totals for the batch; divide by servings.
            def _div(val):
                v = safe_float(val)
                return None if v is None else v / servings
            mapping = [
                ("calories_total", ("energy", "value")),
                ("carbs_g", "carbohydrates"),
                ("net_carbs_g", "net_carbs"),
                ("fat_g", "fat"),
                ("saturated_fat_g", "saturated_fat"),
                ("mono_fat_g", "monounsaturated_fat"),
                ("poly_fat_g", "polyunsaturated_fat"),
                ("trans_fat_g", "trans_fat"),
                ("protein_g", "protein"),
                ("cholest_mg", "cholesterol"),
                ("sodium_mg", "sodium"),
                ("potassium_mg", "potassium"),
                ("sugar_g", "sugar"),
                ("sugar_alcohols_g", "sugar_alcohols"),
                ("added_sugars_g", "added_sugars"),
                ("fiber_g", "fiber"),
                ("calcium_mg", "calcium"),
                ("iron_mg", "iron"),
                ("vitamin_a", "vitamin_a"),
                ("vitamin_c", "vitamin_c"),
                ("vitamin_d", "vitamin_d"),
            ]
            for csv_col, key in mapping:
                v = _div(r.get(csv_col, ""))
                if v is None:
                    continue
                if isinstance(key, tuple):
                    raw.setdefault(key[0], {})[key[1]] = v
                    if key[0] == "energy":
                        raw["energy"]["unit"] = "calories"
                else:
                    raw[key] = v
            # Recipes are per-serving (grams=1 equivalent)
            raw["grams"] = 1.0
            raw["nutrition_source_other"] = r.get("fdc_link", "")
            yield {
                "mfp_id": rid,
                "brand": "",  # recipes have no brand
                "description": r.get("name", ""),
                "verified": False,
                "source_file": path.name,
                "fetched_at": r.get("fetched_at", ""),
                "raw": raw,
                "is_recipe": True,
            }


def load_mine():
    """Load every foods table + recipe table into one mine cache, deduped.

    Sources unioned:
      - mfp_my_foods*.csv       — user's own food library (/food/mine dumps)
      - mfp_foods_from_recipes_*.csv — ingredient foods discovered via /v2/recipes
      - mfp_my_recipes_*.csv    — recipes themselves, treated as foods keyed
                                  by recipe_id with per-serving macros
    When the same id appears in multiple files, the most recent fetched_at
    wins. Recipe-ids and mfp-ids live in the same integer namespace; no collision.
    """
    paths = sorted(
        list(REPO.glob("mfp_my_foods*.csv"))
        + list(REPO.glob("mfp_foods_from_recipes_*.csv"))
    )
    # exclude servings/ingredients CSVs
    paths = [p for p in paths
             if "servings" not in p.name and "ingredients" not in p.name]
    if not paths and not list(REPO.glob("mfp_my_recipes_*.csv")):
        return []
    by_id = {}
    for p in paths:
        with open(p) as f:
            for r in csv.DictReader(f):
                raw = {}
                for csv_col, key in MINE_FIELD_MAP.items():
                    v = safe_float(r.get(csv_col, ""))
                    if v is None:
                        continue
                    if isinstance(key, tuple):
                        raw.setdefault(key[0], {})[key[1]] = v
                        if key[0] == "energy":
                            raw["energy"]["unit"] = "calories"
                    else:
                        raw[key] = v
                raw["nutrition_source_other"] = r.get("fdc_link", "")
                food = {
                    "mfp_id": r["mfp_id"],
                    "brand": r["brand"],
                    "description": r["description"],
                    "verified": r.get("verified") == "True",
                    "source_file": p.name,
                    "fetched_at": r.get("fetched_at", ""),
                    "raw": raw,
                }
                existing = by_id.get(r["mfp_id"])
                if existing is None or food["fetched_at"] > existing["fetched_at"]:
                    by_id[r["mfp_id"]] = food

    # Also load recipes (name-keyed for matching; recipe_id used as mfp_id)
    for p in sorted(REPO.glob("mfp_my_recipes_*.csv")):
        if "ingredients" in p.name:
            continue
        for food in _recipes_as_foods(p):
            existing = by_id.get(food["mfp_id"])
            if existing is None or food["fetched_at"] > existing["fetched_at"]:
                by_id[food["mfp_id"]] = food
    return list(by_id.values())


def find_in_mine(mine, parsed, diary_row):
    """Return (score, tier, food, scale, raw, notes) or None.

    For recipes (is_recipe=True) with exact name match, allow larger
    macro drift — recipe owners edit ingredient quantities over years
    while the recipe identity remains stable. Catches the "slightly
    modified recipe" case without letting truly-different recipes sneak
    through (structural drift remains bounded).
    """
    best = None
    for f in mine:
        sc = score_match(parsed["brand"], parsed["name"], f["brand"], f["description"])
        if sc < APPROX_SCORE_THRESHOLD:
            continue
        scale, s_checked, s_probs, i_probs = verify_by_calorie_ratio(diary_row, f["raw"])
        # Recipes at exact name match: allow larger macro drift to catch
        # owner edits over time. BUT require the scale to roughly equal
        # diary_qty (recipes are per-serving). A reformulation that doubles
        # the recipe's calories would push scale far off diary_qty — reject.
        # For recipes with exact name match, accept reformulation drift.
        # Recipe owners edit ingredient quantities often; the diary-log is a
        # snapshot of that vintage. Structural drift (≤1 mismatch allowed)
        # catches genuinely-different recipes that happen to share a name —
        # e.g. Sag Paneer that differs on 3+ macros will still reject.
        # scale and diary_qty are stored separately so downstream can filter
        # reformulation-drift matches if it wants only-current-version macros.
        is_edited_recipe = f.get("is_recipe") and sc >= 1.0
        max_struct = 1 if is_edited_recipe else 0
        max_inc = 4 if is_edited_recipe else MAX_INCIDENTAL_MISMATCH
        if not verify_ok(s_checked, s_probs, i_probs,
                         max_structural=max_struct, max_incidental=max_inc):
            continue
        tier = 1 if sc >= 1.0 else 2
        notes = (f"mine cache{'  recipe' if f.get('is_recipe') else ''}; "
                 f"{s_checked} structural ({len(s_probs)} off), "
                 f"{len(i_probs)} incidental off (score={sc}, scale={scale:.4f})")
        if s_probs or i_probs:
            notes += "; " + "; ".join(
                f"{p[0]}:{p[1]}!={p[2]}" for p in (s_probs + i_probs))
        cand = (sc, tier, f, scale, f["raw"], notes)
        if best is None or sc > best[0]:
            best = cand
            if sc >= 1.0:
                break
    return best


# --- Per-food processing ---

def fetch_and_verify(client, hit, diary_row, retry_on_error=True):
    """Fetch details with one retry on transient API failure."""
    polite_sleep()
    last_exc = None
    for attempt in range(2 if retry_on_error else 1):
        try:
            det = client.get_food_item_details(hit.mfp_id)
            raw = det.details or {}
            scale, s_checked, s_probs, i_probs = verify_by_calorie_ratio(diary_row, raw)
            return det, raw, scale, s_checked, s_probs, i_probs
        except Exception as e:
            last_exc = e
            time.sleep(random.uniform(4.0, 8.0))
    raise last_exc


def batch_fetch_candidates(client, mfp_ids, chunk_size=30):
    """Fetch full /v2/foods details for many mfp_ids in chunks of chunk_size.

    Returns {str(mfp_id): item_dict}. Uses the library's authed session.
    A single API call returns up to chunk_size items — far cheaper than
    iterating get_food_item_details once per candidate.
    """
    from urllib.parse import urljoin
    out = {}
    fields = ("&fields%5B%5D=nutritional_contents"
              "&fields%5B%5D=serving_sizes"
              "&fields%5B%5D=brand_name"
              "&fields%5B%5D=description"
              "&fields%5B%5D=verified")
    ids = [str(i) for i in mfp_ids]
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i + chunk_size]
        query = "&".join(f"ids%5B%5D={mid}" for mid in chunk)
        url = urljoin(client.BASE_API_URL, "/v2/foods") + "?" + query + fields
        polite_sleep()
        try:
            r = client._get_request_for_url(url, send_token=True)
            if r.ok:
                for item in r.json().get("items", []):
                    out[str(item["id"])] = item
        except Exception:
            continue
    return out


# --- Scale-cleanness: reward scales that match simple rational fractions ---
# For "1 bagel" diary entry, the true scale is a small-denominator rational
# (1/1 for an exact match, 1/2 if MFP is per-2-bagel, 2/1 if MFP is per-half,
# etc.). A scale like 0.963 is a red flag: it signals we have the wrong SKU.

def scale_cleanness(scale, max_denom=6):
    """Returns 1.0 for perfectly clean scale (N/D, small N/D), 0.0 for non-rational.

    Linearly decays across 5% relative difference from nearest simple rational.
    """
    if scale is None or scale <= 0:
        return 0.0
    from fractions import Fraction
    try:
        frac = Fraction(scale).limit_denominator(max_denom)
        closest = float(frac)
    except (ValueError, ZeroDivisionError):
        return 0.0
    if closest == 0:
        return 0.0
    rel_diff = abs(scale - closest) / scale
    return max(0.0, 1.0 - rel_diff * 20.0)  # 1.0 at 0% diff, 0.0 at 5%+


# Unit synonyms for serving-match detection. Intentionally small — not a
# general normalization layer, just grouping the common collisions MFP uses.
UNIT_SYNONYMS = {
    "cup": ["cup", "cups", "cup(s)"],
    "tbsp": ["tbsp", "tablespoon", "tablespoons", "tbs", "tbsp(s)"],
    "tsp": ["tsp", "teaspoon", "teaspoons", "tsp(s)"],
    "oz": ["oz", "ounce", "ounces", "oz(s)"],
    "fl oz": ["fl oz", "fluid ounce", "fluid ounces"],
    "g": ["g", "gram", "grams"],
    "ml": ["ml", "milliliter", "milliliters"],
    "slice": ["slice", "slices", "slice(s)"],
    "piece": ["piece", "pieces"],
    "bagel": ["bagel", "bagels"],
    "bar": ["bar", "bars"],
    "shell": ["shell", "shells"],
    "serving": ["serving", "servings"],
    "leaf": ["leaf", "leaves"],
}


def _canon_unit(u):
    if not u:
        return ""
    u = u.strip().lower()
    for canon, alts in UNIT_SYNONYMS.items():
        if u in alts:
            return canon
    return u


def atwater_drift(raw_nut):
    """Compute |stated_cal - Atwater(macros)| / stated_cal.

    Atwater: 4*(carbs-fiber) + 9*fat + 4*protein + 2*fiber + 2.4*sugar_alcohols.
    Returns None if stated_cal <= 0. Low drift = clean entry; high drift
    (>20%) suggests user-error or unrecorded alcohol/sugar-alcohols.
    """
    cal = safe_float((raw_nut.get("energy") or {}).get("value"))
    if cal is None or cal <= 0:
        return None
    c = safe_float(raw_nut.get("carbohydrates")) or 0.0
    f = safe_float(raw_nut.get("fat")) or 0.0
    p = safe_float(raw_nut.get("protein")) or 0.0
    fib = safe_float(raw_nut.get("fiber")) or 0.0
    sa = safe_float(raw_nut.get("sugar_alcohols")) or 0.0
    # Fiber isn't fully absorbed — give it 2 cal/g (standard US labeling)
    predicted = 4.0 * max(0.0, c - fib) + 9.0 * f + 4.0 * p + 2.0 * fib + 2.4 * sa
    return abs(cal - predicted) / cal


def atwater_cleanness(raw_nut):
    """Score [0, 1]. 1.0 when stated cals match macros via Atwater; decays
    linearly to 0.0 at 20% drift. For tie-breaking, not hard-rejection.
    """
    drift = atwater_drift(raw_nut)
    if drift is None:
        return 0.5  # neutral — can't evaluate
    return max(0.0, 1.0 - drift * 5.0)


def unit_serving_match(raw_nut, serving_sizes, diary_qty, diary_unit_string, cal_scale):
    """Check if any MFP serving option matches the diary's unit.

    For a match, the serving's implied calories-per-diary-portion must agree
    with the calorie-ratio scale within 3%. Returns a cleanness in [0, 1]
    where 1.0 is a perfect unit-and-cal match. Returns None if no candidate
    serving had a unit match (not applicable).

    Example: diary "1 cup Crispix" at 110 cal, MFP entry with a "cup"
    serving (value=1, multiplier=1.0) and default_cal=110 → unit_scale=1.0,
    matches cal_scale=1.0, cleanness=1.0. Whereas an MFP entry with "3/4 cup"
    base (value=0.75, mult=1.0) at default 150 cal → its "cup" serving (if
    present) has value=1.0, mult=1.333 → unit_scale=1.333, cal_scale=0.733.
    Ratio 1.82 off → cleanness=0.
    """
    if not diary_unit_string or not serving_sizes:
        return None
    # Extract the first "unit" word from the diary unit string (strip trailing parens)
    import re
    diary_unit_raw = re.split(r"[\s(]", diary_unit_string.strip(), 1)[0]
    target = _canon_unit(diary_unit_raw)
    if not target:
        return None

    best_cleanness = None
    for sv in serving_sizes:
        sv_canon = _canon_unit(sv.get("unit", ""))
        if sv_canon != target:
            continue
        sv_value = safe_float(sv.get("value")) or 1.0
        sv_mult = safe_float(sv.get("nutrition_multiplier")) or 1.0
        if sv_value == 0:
            continue
        # Unit-derived scale: to scale MFP default macros to diary portion
        unit_scale = diary_qty * sv_mult / sv_value
        if unit_scale <= 0:
            continue
        if cal_scale is None or cal_scale <= 0:
            continue
        ratio_diff = abs(unit_scale - cal_scale) / max(unit_scale, cal_scale)
        # Decay across 10% diff (MFP "1 cup" multipliers are user-entered
        # approximations like 0.75188 where the real value is 0.75).
        cleanness = max(0.0, 1.0 - ratio_diff * 10)
        if best_cleanness is None or cleanness > best_cleanness:
            best_cleanness = cleanness
    return best_cleanness


def process_food(client, food_name, diary_row, mine, override_id=None):
    parsed = parse_diary_name(food_name)
    if parsed["err"]:
        return build_row(food_name, parsed, "", STATUS_PARSE, "", "",
                         None, None, None, None, None, None, parsed["err"])
    if safe_float(diary_row.get("calories")) in (None, 0.0):
        return build_row(food_name, parsed, "", STATUS_NO_CALS, "", "",
                         None, None, None, None, None, None,
                         f"diary calories={diary_row.get('calories')!r}")

    # Manual override path. Human has identified the food; we trust it.
    # We still record any macro drift in notes for transparency, but accept
    # the match regardless (scale gets used verbatim, so minor reformulation
    # differences just mean the enrichment carries small structural drift).
    if override_id:
        try:
            polite_sleep()
            det = client.get_food_item_details(int(override_id))
            raw = det.details or {}
            scale, s_checked, s_probs, i_probs = verify_by_calorie_ratio(diary_row, raw)
            cleanness = scale_cleanness(scale)
            notes = (f"manual override (trusted); scale={scale:.4f} "
                     f"cleanness={cleanness:.2f}; {s_checked} structural "
                     f"({len(s_probs)} off), {len(i_probs)} incidental off")
            if s_probs or i_probs:
                notes += "; " + "; ".join(
                    f"{p[0]}:{p[1]}!={p[2]}" for p in (s_probs + i_probs))
            return build_row(food_name, parsed, "override",
                             STATUS_RESOLVED_OVERRIDE, 1, 1.0,
                             int(override_id), det.brand, det.name, det.verified,
                             scale, raw, notes, diary_row=diary_row)
        except Exception as e:
            return build_row(food_name, parsed, "override", STATUS_NOT_FOUND,
                             "", "", int(override_id), None, None, None, None,
                             None, f"override fetch error: {type(e).__name__}: {e}")

    # 1) Mine cache check first — zero MFP requests
    if mine:
        m = find_in_mine(mine, parsed, diary_row)
        if m:
            sc, tier, f, scale, raw, notes = m
            status = (STATUS_RESOLVED if tier == 1 else
                      STATUS_RESOLVED_FROM_MINE)
            return build_row(food_name, parsed, "mine", status, tier, sc,
                             f["mfp_id"], f["brand"], f["description"],
                             f["verified"], scale, raw, notes, diary_row=diary_row)

    # 2) MFP search — with retry using a compacted query if first misses
    def _build_queries():
        """Ordered list of query strings to try. Stops at first verifying match."""
        qs = []
        name = parsed["name"]
        brand = parsed["brand"]
        if brand and brand.lower() != name.lower():
            qs.append(f"{brand} {name}")
        else:
            qs.append(name)
        # One general fallback: brand + last 3 name tokens. This drops
        # category prefixes without special-casing any particular punctuation.
        # "Salad Dressing: Catalina, Light" tokens → last 3 = "Dressing Catalina Light".
        # For truly problematic foods, write an explicit override (see
        # mfp_food_overrides.csv) rather than extending this heuristic.
        name_tokens = re.findall(r"[A-Za-z0-9%]+", name)
        if len(name_tokens) > 3:
            last3 = " ".join(name_tokens[-3:])
            q = f"{brand} {last3}" if brand else last3
            if q.lower() not in {x.lower() for x in qs}:
                qs.append(q)
        return qs

    queries = _build_queries()
    all_results = []
    seen = set()
    scored = []
    last_error = None
    for qi, query in enumerate(queries):
        # Retry on ConnectionError up to 3x with exponential backoff
        results = None
        for attempt in range(3):
            try:
                results = client.get_food_search_results(query)
                break
            except Exception as e:
                last_error = f"search error on query {qi} ({query!r}): {type(e).__name__}: {e}"
                err_type = type(e).__name__
                if err_type in ("ConnectionError", "Timeout", "ConnectTimeout",
                                "ReadTimeout", "MaxRetryError"):
                    if attempt < 2:
                        time.sleep(random.uniform(5.0, 15.0) * (attempt + 1))
                        continue
                break
        if results is None:
            continue
        for r in results:
            if r.mfp_id in seen:
                continue
            seen.add(r.mfp_id)
            sc = score_match(parsed["brand"], parsed["name"], r.brand, r.name)
            if sc >= APPROX_SCORE_THRESHOLD:
                scored.append((sc, r))
        all_results.extend(results)
        # If first query already found exact-name candidates, don't bother refetching
        if any(sc >= 1.0 for sc, _ in scored):
            break
        # After first query, if we have plenty of high-score candidates, stop
        if len(scored) >= MAX_CANDIDATES_TO_VERIFY and qi > 0:
            break
        if qi + 1 < len(queries):
            polite_sleep()

    results = all_results  # for diagnostics

    if not scored and last_error:
        return build_row(food_name, parsed, "search", STATUS_NOT_FOUND,
                         "", "", None, None, None, None, None, None,
                         last_error)

    if not scored:
        top = [f"'{r.brand}|{r.name}'" for r in results[:5]]
        return build_row(food_name, parsed, "search", STATUS_NOT_FOUND,
                         "", "", None, None, None, None, None, None,
                         f"no candidate ≥ {APPROX_SCORE_THRESHOLD}; top5={top}")

    # Two passes:
    #   Pass A: score == 1.0 (exact identity), verify each → return tier 1
    #   Pass B: remaining, ordered by (score desc), verify → tier 2 if
    #           candidate.verified else tier 3. Track best.
    exact = [(s, r) for s, r in scored if s >= 1.0]
    approx = sorted([(s, r) for s, r in scored if s < 1.0], key=lambda x: -x[0])

    fetched = 0
    tried = []
    fetch_errors = []

    # Calorie-aware candidate hunting:
    # 1. Batch-fetch full details for ALL scored candidates (up to 100 —
    #    MFP's search ceiling). One /v2/foods call per 30-id chunk.
    # 2. KEEP deleted entries — old SKUs from 2011-era are often flagged
    #    deleted by MFP and those are often the best match for old logs.
    # 3. For each, compute scale, unit-match cleanness, Atwater cleanness.
    # 4. Rank by (structural_verified, cleanness, atwater, USDA, name_score).
    # 5. Accept the top structural-verified candidate.
    BATCH_CAP = 100
    by_name = (exact + approx)[:BATCH_CAP]
    candidate_ids = [c.mfp_id for _, c in by_name]

    batch = batch_fetch_candidates(client, candidate_ids)
    tried = list(batch.keys())

    ranked = []  # tuples for sorting
    rows_by_id = {}
    for sc, hit in by_name:
        item = batch.get(str(hit.mfp_id))
        if not item:
            continue
        raw = item.get("nutritional_contents") or {}
        try:
            scale, s_checked, s_probs, i_probs = verify_by_calorie_ratio(
                diary_row, raw)
        except Exception:
            continue
        rational_cleanness = scale_cleanness(scale)
        # Unit-based cleanness — if MFP entry exposes a serving with the
        # diary's unit AND its implied scale agrees with the calorie scale,
        # that's a strong identity signal. Use it when available.
        serving_sizes = item.get("serving_sizes") or []
        unit_cleanness = unit_serving_match(
            raw, serving_sizes, parsed["qty"], parsed["unit_string"], scale)
        # Combined cleanness: use unit signal when present, else rational.
        if unit_cleanness is not None:
            cleanness = max(unit_cleanness, rational_cleanness)
        else:
            cleanness = rational_cleanness
        # Atwater consistency: does the entry's stated cals match macros?
        # Low drift = high-quality source. Tie-breaker for ranking.
        atwater = atwater_cleanness(raw)
        verified = item.get("verified", False)

        # Is this a per-gram/per-100g entry? MFP's `grams` field in
        # nutritional_contents tells us: 1 = per-serving base, 100 = per-100g.
        # For per-100g/per-gram entries, the scale is just the gram weight
        # ratio and doesn't need to be a simple rational.
        grams_base = safe_float(raw.get("grams"))
        is_gram_based = grams_base is not None and grams_base >= 50

        # Structural tolerance: normally 0 mismatches. Relax to 1 mismatch
        # when exact name + any identity signal supports "same food, slight
        # drift". Identity signals (any one qualifies):
        #  - cleanness >= 0.9 (scale is a simple rational)
        #  - unit_cleanness >= 0.8 (MFP has a serving matching diary unit,
        #    and its implied scale agrees with calorie scale within 2%)
        #  - is_gram_based (per-100g USDA entry; scale is raw gram ratio)
        #  - USDA-verified AND cleanness >= 0.7 (wider margin for trusted data)
        # Drift-protection guard below prevents a wildly-off single macro
        # from corrupting downstream daily totals (stored as blank → diary wins).
        high_confidence = sc >= 1.0 and (
            cleanness >= 0.9
            or is_gram_based
            or (unit_cleanness is not None and unit_cleanness >= 0.8)
            or (verified and cleanness >= 0.7)
        )
        max_struct = 1 if high_confidence else 0
        structural_ok = len(s_probs) <= max_struct and s_checked >= 2
        incidental_ok = len(i_probs) <= MAX_INCIDENTAL_MISMATCH
        accept = structural_ok and incidental_ok

        # Rank key: prefer accept, then cleanness, then USDA, then Atwater,
        # then name_score. Atwater is a tiebreaker — lets us pick cleaner
        # entries when name/cleanness tie.
        ranked.append((
            -int(accept),             # True=-1, False=0 (sort asc, True first)
            -cleanness,               # higher = cleaner
            -int(verified),           # USDA first
            -atwater,                 # cleaner Atwater first
            -sc,                      # name score
            hit.mfp_id,
            sc, hit, item, scale, s_checked, s_probs, i_probs, raw,
        ))

    ranked.sort(key=lambda t: t[:6])

    # Return first accepted (accept=True sorts to the head via accept_neg=-1).
    for rec in ranked:
        accept_neg, cn_neg, ver_neg, atw_neg, sc_neg, mid, sc, hit, item, scale, s_checked, s_probs, i_probs, raw = rec
        if accept_neg == 0:
            break  # hit the first non-accepted — rest are also non-accepted
        verified = bool(-ver_neg)
        cleanness = -cn_neg

        prob_notes = ""
        if s_probs or i_probs:
            prob_notes = "; " + "; ".join(
                f"{p[0]}:{p[1]}!={p[2]}" for p in (s_probs + i_probs))

        # Determine tier
        if sc >= 1.0:
            status = STATUS_RESOLVED
            tier = 1
        elif verified:
            status = STATUS_RESOLVED_APPROX_VERIFIED
            tier = 2
            if len(i_probs) > 0:
                # For approximate names, demand tight incidental match too
                # (tier 2 budget is stricter)
                continue
        else:
            status = STATUS_RESOLVED_APPROX
            tier = 3
            if len(i_probs) > 0:
                continue

        notes = (f"tier{tier}, score={sc}, scale={scale:.4f} "
                 f"(cleanness={cleanness:.2f}), {s_checked} structural + "
                 f"{len(i_probs)} incidental off" + prob_notes)
        brand = item.get("brand_name") or ""
        name = item.get("description", hit.name)
        return build_row(food_name, parsed, "search",
                         status, tier, sc,
                         hit.mfp_id, brand, name, verified,
                         scale, raw, notes, diary_row=diary_row)

    # Nothing scalar-verified — return the closest candidate diagnostic
    diag = ", ".join(
        f"({sc}|{r.brand}|{r.name})" for sc, r in (exact + approx)[:5])
    return build_row(food_name, parsed, "search", STATUS_NOT_FOUND,
                     "", "", None, None, None, None, None, None,
                     f"tried {fetched} candidates [{','.join(tried)}]; "
                     f"none scalar-verified. Top scored: {diag}; errors={fetch_errors}")


# --- Persistence ---

def load_dict():
    if not DICT_PATH.exists():
        return {}
    with open(DICT_PATH) as f:
        return {r["diary_food_name"]: r for r in csv.DictReader(f)}


def save_dict(d):
    tmp = DICT_PATH.with_suffix(".csv.tmp")
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        for name in sorted(d):
            w.writerow({c: d[name].get(c, "") for c in COLUMNS})
    tmp.replace(DICT_PATH)


def load_intake():
    counts = Counter()
    sample = {}
    with open(INTAKE_PATH) as f:
        for r in csv.DictReader(f):
            counts[r["food"]] += 1
            sample.setdefault(r["food"], r)
    return counts, sample


# --- CLI ---

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--limit", type=int, default=10)
    p.add_argument("--force", type=str, default=None)
    p.add_argument("--status", type=str, default=None)
    p.add_argument("--mfp-id", type=str, default=None)
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    existing = load_dict()
    counts, samples = load_intake()
    mine = load_mine()
    print(f"Mine cache: {len(mine)} foods", file=sys.stderr)

    if args.force:
        if args.force not in counts:
            print(f"ERROR: {args.force!r} not in intake_foods.csv", file=sys.stderr)
            sys.exit(2)
        targets = [args.force]
    elif args.status:
        targets = sorted(
            (n for n, r in existing.items() if r.get("status") == args.status),
            key=lambda n: -counts.get(n, 0))[:args.limit]
    else:
        # Default: only process foods that have NEVER been attempted.
        # Previously-tried rows (match_not_found, parse_fail, diary_zero_calories)
        # are skipped to avoid wasting requests re-hunting the same persistent
        # misses. To retry those, use --status match_not_found explicitly.
        targets = []
        for name, _ in counts.most_common():
            if name not in existing:
                targets.append(name)
            if len(targets) >= args.limit:
                break

    print(f"Targets: {len(targets)}", file=sys.stderr)
    for t in targets:
        print(f"  [{counts[t]:5d}]  {t}", file=sys.stderr)

    if args.dry_run:
        for t in targets:
            print(f"\n--- {t!r}\n  parsed: {parse_diary_name(t)}")
        return

    client = make_client()
    print(f"[auth] user_id={client.user_id} username={USERNAME}", file=sys.stderr)

    # Session stats: running totals + rolling window
    total_rows = sum(counts.values())
    resolved_rows_before = sum(
        counts.get(n, 0) for n, r in existing.items()
        if r.get("status", "").startswith("resolved"))
    session_resolved = 0
    session_not_found = 0
    session_other = 0
    recent_window = []  # list of bool: True=resolved, False=not

    def _coverage_pct():
        total_resolved_now = sum(
            counts.get(n, 0) for n, r in existing.items()
            if r.get("status", "").startswith("resolved"))
        return total_resolved_now / total_rows * 100

    for idx, t in enumerate(targets, 1):
        print(f"\n=== [{idx}/{len(targets)}] {t!r} (count={counts[t]}) ===")
        try:
            row = process_food(client, t, samples[t], mine,
                               override_id=args.mfp_id)
        except Exception as e:
            row = {c: "" for c in COLUMNS}
            row["diary_food_name"] = t
            row["status"] = "error"
            row["notes"] = f"{type(e).__name__}: {e}"
            row["fetched_at"] = dt.datetime.utcnow().isoformat(timespec="seconds") + "Z"
            print(f"  UNEXPECTED {row['notes']}")

        existing[t] = row
        save_dict(existing)

        # Track per-food outcome
        is_resolved = row["status"].startswith("resolved")
        if is_resolved:
            session_resolved += 1
        elif row["status"] == "match_not_found":
            session_not_found += 1
        else:
            session_other += 1
        recent_window.append(is_resolved)
        if len(recent_window) > 10:
            recent_window.pop(0)

        # Track genuine API errors for session-expiry / network detection.
        # match_not_found without 'search error' in notes just means the food
        # isn't in MFP (legitimate miss). match_not_found WITH search error in
        # notes is a connection/auth failure that should halt the batch.
        notes_lower = row.get("notes", "").lower()
        is_api_error = (row["status"] == "error"
                        or "search error" in notes_lower
                        or "connectionerror" in notes_lower
                        or "max retries" in notes_lower)

        print(f"  status   = {row['status']}  tier={row.get('match_tier')}  "
              f"score={row.get('match_score')}  source={row.get('source')}")
        if row.get("mfp_id"):
            print(f"  mfp_id   = {row['mfp_id']}  ({row['mfp_brand']} | {row['mfp_name']})  "
                  f"verified={row.get('mfp_verified')}")
            print(f"  scale    = {row['scale_factor']}")
        if row.get("notes"):
            print(f"  notes    = {row['notes']}")

        # Live progress line
        n_done = idx
        hit_rate = session_resolved / n_done * 100
        recent_rate = (sum(recent_window) / len(recent_window) * 100) if recent_window else 0
        print(f"  PROGRESS {n_done}/{len(targets)}  "
              f"session: {session_resolved}R {session_not_found}NF {session_other}O  "
              f"hit: {hit_rate:.0f}%  last-{len(recent_window)}: {recent_rate:.0f}%  "
              f"coverage: {_coverage_pct():.2f}%",
              flush=True)

        # Auto-halt on consecutive API errors (genuine session expiry).
        # A run of match_not_found is legitimate — many diary foods are
        # custom / not in MFP's public DB.
        if is_api_error:
            recent_api_errors = getattr(main, "_api_err_count", 0) + 1
        else:
            recent_api_errors = 0
        main._api_err_count = recent_api_errors
        if recent_api_errors >= 3:
            print(f"\n!!! {recent_api_errors} consecutive API errors — session may "
                  "have expired. Halting. Log in to MFP and rerun.",
                  file=sys.stderr, flush=True)
            break

        if row.get("source") != "mine":
            polite_sleep()


if __name__ == "__main__":
    main()
