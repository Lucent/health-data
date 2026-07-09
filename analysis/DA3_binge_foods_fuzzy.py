"""DA3. Binge food analysis with fuzzy grouping across reformulations and portions.

Same as DA2 but consolidates entries that differ only by:
  - Portion size (last comma-separated segment)
  - Parenthetical notes ("(corrected)", "(minus erythritol)", "(old)")
  - Date stamps ("Oct 2015")
  - "old" / "new" / "(updated)" markers
  - Asterisks, tildes
  - "morning star" vs "morningstar"

Each cluster's contents are printed so we can verify the grouping makes sense
before drawing conclusions.
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
INTAKE_MULTIPLIER = 1.09


def food_key(food):
    """Reduce a food name to a normalized cluster key."""
    s = str(food).lower().strip()
    # Strip portion: everything after the LAST comma
    if "," in s:
        s = s.rsplit(",", 1)[0].strip()
    # Remove parenthetical notes
    s = re.sub(r"\([^)]*\)", "", s)
    # Remove month-year date stamps
    s = re.sub(r"\b(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)[a-z]*\s*\d{2,4}",
                "", s)
    # Remove old/new/corrected/limited markers
    s = re.sub(r"\b(old|new|corrected|approx|limited edition|updated|original)\b",
                "", s)
    # Normalize "morning star" / "morningstar"
    s = re.sub(r"morning\s*star", "morningstar", s)
    # Strip special chars
    s = re.sub(r"[*~]", "", s)
    # Strip trailing dashes / commas / colons
    s = re.sub(r"[\-,:;]+$", "", s).strip()
    # Collapse whitespace
    s = re.sub(r"\s+", " ", s).strip()
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
    print(f"Days analyzed: {len(pre)} ({pre['date'].min().date()} to {pre['date'].max().date()})")
    print(f"Total binges in window: {pre['binge'].sum()} ({pre['binge'].mean()*100:.1f}%)")

    # Apply fuzzy key
    foods = foods.merge(pre[["date", "binge", "year"]], on="date", how="inner")
    foods["fuzzy_key"] = foods["food"].apply(food_key)

    # For each fuzzy_key, count distinct DAYS with that food (deduplicate within a day)
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

    # Same long-term filter
    long_term = agg[(agg["n"] >= 50) & (agg["n_years"] >= 3) & (agg["span_days"] >= 700)].copy()
    print(f"\nFuzzy clusters passing filter (n≥50, ≥3 years, span≥700d): {len(long_term)}")

    # Map keys back to constituent food names for transparency
    name_map = (foods.groupby("fuzzy_key")["food"]
                 .agg(lambda s: sorted(set(s)))
                 .to_dict())

    # ── Top triggers ──
    print("\n" + "=" * 110)
    print("TOP TRIGGER CLUSTERS (lift > 1.5, sorted by lift)")
    print("=" * 110)
    triggers = long_term[long_term["lift"] > 1.5].sort_values("lift", ascending=False)
    for _, r in triggers.iterrows():
        print(f"\n  KEY: '{r['fuzzy_key']}'")
        print(f"    n={int(r['n'])}, binge%={r['binge_share']*100:.1f}%, "
              f"lift={r['lift']:.2f}, years={int(r['n_years'])}, "
              f"span={r['span_days']/365:.1f}y "
              f"({r['first_date'].strftime('%Y-%m')} → {r['last_date'].strftime('%Y-%m')})")
        names = name_map.get(r["fuzzy_key"], [])
        print(f"    Variants ({len(names)}):")
        for nm in names[:8]:
            print(f"      • {nm[:80]}")
        if len(names) > 8:
            print(f"      …and {len(names) - 8} more")

    # ── Helpful clusters ──
    print("\n" + "=" * 110)
    print("HELPFUL CLUSTERS (lift ≤ 0.40, sorted by ascending lift)")
    print("=" * 110)
    helpful = long_term[long_term["lift"] <= 0.40].sort_values("lift").head(20)
    for _, r in helpful.iterrows():
        print(f"\n  KEY: '{r['fuzzy_key']}'")
        print(f"    n={int(r['n'])}, binge%={r['binge_share']*100:.1f}%, "
              f"lift={r['lift']:.2f}, years={int(r['n_years'])}, "
              f"span={r['span_days']/365:.1f}y "
              f"({r['first_date'].strftime('%Y-%m')} → {r['last_date'].strftime('%Y-%m')})")
        names = name_map.get(r["fuzzy_key"], [])
        print(f"    Variants ({len(names)}):")
        for nm in names[:8]:
            print(f"      • {nm[:80]}")
        if len(names) > 8:
            print(f"      …and {len(names) - 8} more")

    # ── Cross-era staples (most years) ──
    print("\n" + "=" * 110)
    print("CLUSTERS APPEARING IN ≥ 7 YEARS (most cross-era), sorted by lift descending")
    print("=" * 110)
    staples = long_term[long_term["n_years"] >= 7].sort_values("lift", ascending=False)
    print(f"\n  {'lift':>5} {'years':>6} {'n':>5} {'binge%':>7} {'cluster':<60}")
    for _, r in staples.iterrows():
        print(f"  {r['lift']:>5.2f} {int(r['n_years']):>6} {int(r['n']):>5} "
              f"{r['binge_share']*100:>6.1f}% {r['fuzzy_key'][:60]:<60}")

    out_path = ROOT / "analysis" / "DA3_binge_foods_fuzzy.csv"
    long_term.sort_values("lift", ascending=False).to_csv(out_path, index=False)
    print(f"\nSaved {out_path}")


if __name__ == "__main__":
    main()
