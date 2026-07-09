"""DA4. Binge food clusters with stronger normalization.

Adds to DA3:
  - Apostrophe stripping (universal): "boar's head" → "boars head"
  - Trader Joe / Trader Joes / Trader Joe's all map to "trader joes"
  - Boarshead / Boar's Head map to "boars head"
  - Ben & Jerry's variants → "ben jerrys"
  - Strip "Per Bag", "Per Container", "ultra-filtered", "ultra-pasteurized"
  - Strip "*correct*", asterisks
  - Generic-food aliases for items where the brand drops in/out:
      * Heavy cream (any brand) → "heavy cream"
      * Maple syrup (pure, any brand) → "maple syrup" — but EXCLUDE ChocZero sugar-free
      * Slate chocolate milk (all variants) → "slate chocolate milk"
  - Other normalizations stay the same as DA3
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
INTAKE_MULTIPLIER = 1.09


def food_key(food):
    """Reduce a food name to a normalized cluster key (more aggressive than DA3)."""
    s = str(food).lower().strip()

    # Strip portion: everything after the LAST comma
    if "," in s:
        s = s.rsplit(",", 1)[0].strip()

    # Remove apostrophes universally
    s = s.replace("'", "").replace("’", "")

    # Remove parentheticals
    s = re.sub(r"\([^)]*\)", "", s)

    # Remove date stamps
    s = re.sub(r"\b(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)[a-z]*\s*\d{2,4}",
                "", s)

    # Reformulation / version markers
    s = re.sub(
        r"\b(old|new|corrected|approx|limited\s+edition|updated|original|"
        r"per\s+bag|per\s+container|v\d+|correct|chickn|chikn|"
        r"chick\s*n|chik\s*n)\b",
        "", s)

    # Product-irrelevant modifiers
    s = re.sub(
        r"\b(ultra[-\s]?filtered|ultra[-\s]?pasteurized|"
        r"grass[-\s]?fed|cage[-\s]?free|extra\s+virgin|"
        r"unsweetened|low[-\s]?fat|reduced[-\s]?fat|fat[-\s]?free)\b",
        "", s)

    # Brand normalizations
    s = re.sub(r"morning\s*star", "morningstar", s)
    s = re.sub(r"trader\s+joe[s]?", "trader joes", s)
    s = re.sub(r"boar[s]?\s*head", "boars head", s)
    s = re.sub(r"ben\s*&?\s*jerry[s]?", "ben jerrys", s)
    s = re.sub(r"hannaford\s*-\s*", "", s)  # generic when brand=Hannaford

    # Strip special chars / leading dashes
    s = re.sub(r"[*~]", "", s)
    s = re.sub(r"^\s*-\s*", "", s)
    s = re.sub(r"\s*-\s*$", "", s)

    # Collapse whitespace
    s = re.sub(r"\s+", " ", s).strip()

    # ─── Specific aliases (after the above scrubbing) ───

    # Slate chocolate milk (4 named variants → 1)
    if re.search(r"\bslate\b", s) and "milk" in s and (
            "chocolate" in s or "dark" in s):
        return "slate chocolate milk"

    # Heavy cream (any brand, any case)
    if re.search(r"^(?:.*-\s*)?heavy\s+cream\s*$", s):
        return "heavy cream"

    # Pure maple syrup (any brand) but NOT sugar-free / ChocZero / Choczero
    if "maple syrup" in s and "choczero" not in s and "sugar-free" not in s and "sugar free" not in s and "honest" not in s:
        return "maple syrup"

    # Generic butter (unsalted variants) — keep kerrygold separate as a deliberate brand choice
    if re.search(r"^butter[\s,]*unsalted\s*$", s) or re.search(r"^unsalted\s+butter\s*$", s):
        return "butter unsalted"
    if re.search(r"^butter\s*$", s):
        return "butter"

    # Generic olive oil
    if re.search(r"^olive\s+oil\s*$", s):
        return "olive oil"

    # Avocado generic
    if re.search(r"^avocado(\s*-\s*avocado)?\s*$", s):
        return "avocado"

    return s


def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    intake["calories"] = intake["calories"] * INTAKE_MULTIPLIER
    kf = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    foods = pd.read_csv(ROOT / "intake" / "intake_foods.csv", parse_dates=["date"])

    df = intake[["date", "calories"]].merge(
        kf[["date", "fat_mass_lbs", "tdee"]], on="date")
    df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df = df.sort_values("date").reset_index(drop=True)
    df["surplus"] = df["calories"] - df["tdee"]
    df["binge"] = (df["surplus"] > 1000).astype(int)
    df["year"] = df["date"].dt.year

    pre = df[(df["effective_level"] == 0) & (df["year"] >= 2014)].copy()
    print(f"Days analyzed: {len(pre)} ({pre['date'].min().date()} → {pre['date'].max().date()})")
    print(f"Binges in window: {pre['binge'].sum()} ({pre['binge'].mean()*100:.1f}%)")

    foods = foods.merge(pre[["date", "binge", "year"]], on="date", how="inner")
    foods["fuzzy_key"] = foods["food"].apply(food_key)

    # Deduplicate within day
    foods_dedupe = foods.drop_duplicates(subset=["date", "fuzzy_key"])

    grp = foods_dedupe.groupby("fuzzy_key")
    agg = grp.agg(
        n=("binge", "size"),
        binge_count=("binge", "sum"),
        first_date=("date", "min"),
        last_date=("date", "max"),
        n_years=("year", "nunique"),
    ).reset_index()
    agg["binge_share"] = agg["binge_count"] / agg["n"]
    agg["span_days"] = (agg["last_date"] - agg["first_date"]).dt.days
    overall = pre["binge"].mean()
    agg["lift"] = agg["binge_share"] / overall

    long_term = agg[(agg["n"] >= 50) & (agg["n_years"] >= 3) & (agg["span_days"] >= 700)].copy()
    print(f"\nClusters passing filter (n≥50, ≥3 years, span≥700d): {len(long_term)}")

    name_map = (foods.groupby("fuzzy_key")["food"]
                 .agg(lambda s: sorted(set(s)))
                 .to_dict())
    n_orig_per_cluster = (foods.groupby("fuzzy_key")["food"].nunique().to_dict())

    print("\n" + "=" * 110)
    print("TOP TRIGGERS (lift > 1.5)")
    print("=" * 110)
    triggers = long_term[long_term["lift"] > 1.5].sort_values("lift", ascending=False)
    for _, r in triggers.iterrows():
        print(f"\n  KEY: '{r['fuzzy_key']}'")
        print(f"    n={int(r['n'])}, binge%={r['binge_share']*100:.1f}%, "
              f"lift={r['lift']:.2f}, years={int(r['n_years'])}, "
              f"span={r['span_days']/365:.1f}y "
              f"({r['first_date'].strftime('%Y-%m')} → {r['last_date'].strftime('%Y-%m')}), "
              f"{n_orig_per_cluster.get(r['fuzzy_key'], 1)} orig variants")

    print("\n" + "=" * 110)
    print("HELPFUL CLUSTERS (lift ≤ 0.40)")
    print("=" * 110)
    helpful = long_term[long_term["lift"] <= 0.40].sort_values("lift").head(25)
    for _, r in helpful.iterrows():
        print(f"\n  KEY: '{r['fuzzy_key']}'")
        print(f"    n={int(r['n'])}, binge%={r['binge_share']*100:.1f}%, "
              f"lift={r['lift']:.2f}, years={int(r['n_years'])}, "
              f"span={r['span_days']/365:.1f}y, "
              f"{n_orig_per_cluster.get(r['fuzzy_key'], 1)} orig variants")

    print("\n" + "=" * 110)
    print("CLUSTERS APPEARING IN ≥ 7 YEARS (most cross-era)")
    print("=" * 110)
    staples = long_term[long_term["n_years"] >= 7].sort_values("lift", ascending=False)
    print(f"\n  {'lift':>5} {'years':>6} {'n':>5} {'binge%':>7} {'orig variants':>13}  {'cluster':<60}")
    for _, r in staples.iterrows():
        nv = n_orig_per_cluster.get(r["fuzzy_key"], 1)
        print(f"  {r['lift']:>5.2f} {int(r['n_years']):>6} {int(r['n']):>5} "
              f"{r['binge_share']*100:>6.1f}% {nv:>13}  {r['fuzzy_key'][:60]:<60}")

    # Show specific consolidations that came from the new aliases
    print("\n" + "=" * 110)
    print("SPECIFIC ALIAS CONSOLIDATIONS (verify they merged correctly):")
    print("=" * 110)
    for key in ["slate chocolate milk", "heavy cream", "maple syrup",
                "boars head - chipotle gouda cheese", "boars head - chipotle gouda",
                "avocado"]:
        if key in name_map:
            print(f"\n  '{key}' merged {len(name_map[key])} variants:")
            for nm in name_map[key][:10]:
                print(f"    • {nm[:80]}")
            if len(name_map[key]) > 10:
                print(f"    …and {len(name_map[key]) - 10} more")
        else:
            close_keys = [k for k in name_map if key.split()[0] in k]
            if close_keys:
                print(f"\n  '{key}' not found exactly. Nearby keys:")
                for k in close_keys[:5]:
                    print(f"    • '{k}' (n_variants={len(name_map[k])})")

    out_path = ROOT / "analysis" / "DA4_binge_foods_v2.csv"
    long_term.sort_values("lift", ascending=False).to_csv(out_path, index=False)
    print(f"\nSaved {out_path}")


if __name__ == "__main__":
    main()
