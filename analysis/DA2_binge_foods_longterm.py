"""DA2. Binge-associated foods, restricted to foods that span a long period.

Filters:
  - Exclude 2011-2013 (initial willpower-era rapid loss). Restrict to 2014+ pre-tirz.
  - Foods must appear in ≥ 3 distinct calendar years AND n ≥ 50 days.
  - Binge definition: intake > TDEE + 1000 cal (AG canonical).

Lift = food's binge_rate / overall_binge_rate. >1 = trigger; <1 = helpful (or just neutral).
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
INTAKE_MULTIPLIER = 1.09


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
    df["month_key"] = df["date"].dt.to_period("M").astype(str)

    # Exclude 2011-2013 willpower era + tirz era
    pre = df[(df["effective_level"] == 0) & (df["year"] >= 2014)].copy()
    print(f"Days analyzed: {len(pre)} ({pre['date'].min().date()} to {pre['date'].max().date()})")
    print(f"Total binges in this window: {pre['binge'].sum()} "
          f"({pre['binge'].mean()*100:.1f}% rate)")

    # Merge food log with daily binge flag (only days in our window)
    foods = foods.merge(pre[["date", "binge", "year", "month_key"]],
                         on="date", how="inner")

    # Per-food summary: n_days, distinct years, distinct months, binge_rate
    grp = foods.groupby("food")
    agg = grp.agg(
        n=("binge", "size"),
        binge_count=("binge", "sum"),
        first_date=("date", "min"),
        last_date=("date", "max"),
        n_years=("year", "nunique"),
        n_months=("month_key", "nunique"),
    ).reset_index()
    agg["binge_share"] = agg["binge_count"] / agg["n"]
    agg["span_days"] = (agg["last_date"] - agg["first_date"]).dt.days
    overall = pre["binge"].mean()
    agg["lift"] = agg["binge_share"] / overall

    # Filter: n ≥ 50 AND distinct years ≥ 3 AND span ≥ 700 days (≈ 2 years)
    long_term = agg[(agg["n"] >= 50) & (agg["n_years"] >= 3) & (agg["span_days"] >= 700)].copy()

    print(f"\nFoods passing filter (n≥50, ≥3 distinct years, span≥700d): {len(long_term)}")

    # Less strict: n ≥ 30, ≥ 2 years, span ≥ 400 days — to catch borderline triggers
    medium = agg[(agg["n"] >= 30) & (agg["n_years"] >= 2) & (agg["span_days"] >= 400)].copy()
    print(f"Foods passing relaxed filter (n≥30, ≥2 distinct years, span≥400d): {len(medium)}")
    print("\n" + "=" * 110)
    print("RELAXED-FILTER TRIGGERS (lift > 1.2, sorted by lift):")
    print("=" * 110)
    rel_trig = medium[medium["lift"] > 1.2].sort_values("lift", ascending=False)
    if len(rel_trig) == 0:
        print("  None.")
    else:
        print(f"  {'food':<60} {'n':>5} {'binge%':>7} {'lift':>5} {'years':>6} {'span_y':>7} "
              f"{'first':<11}{'last':<11}")
        for _, r in rel_trig.head(40).iterrows():
            food_short = r["food"][:58]
            print(f"  {food_short:<60} {int(r['n']):>5} {r['binge_share']*100:>6.1f}% "
                  f"{r['lift']:>5.2f} {int(r['n_years']):>6} {r['span_days']/365:>7.1f} "
                  f"{r['first_date'].strftime('%Y-%m'):<8}    {r['last_date'].strftime('%Y-%m'):<8}")

    # ── Top trigger foods ──
    print("\n" + "=" * 110)
    print("TOP TRIGGER FOODS (lift > 1.5, sorted by lift)")
    print("=" * 110)
    triggers = long_term[long_term["lift"] >= 1.5].sort_values("lift", ascending=False)
    print(f"  {'food':<60} {'n':>5} {'binge%':>7} {'lift':>5} {'years':>6} {'span_y':>7} "
          f"{'first':<11}{'last':<11}")
    for _, r in triggers.head(40).iterrows():
        food_short = r["food"][:58]
        print(f"  {food_short:<60} {int(r['n']):>5} {r['binge_share']*100:>6.1f}% "
              f"{r['lift']:>5.2f} {int(r['n_years']):>6} {r['span_days']/365:>7.1f} "
              f"{r['first_date'].strftime('%Y-%m'):<8}    {r['last_date'].strftime('%Y-%m'):<8}")

    # ── Helpful foods (low binge association) ──
    print("\n" + "=" * 110)
    print("HELPFUL FOODS (lift ≤ 0.50, sorted by ascending lift)")
    print("=" * 110)
    helpful = long_term[long_term["lift"] <= 0.50].sort_values("lift")
    print(f"  {'food':<60} {'n':>5} {'binge%':>7} {'lift':>5} {'years':>6} {'span_y':>7} "
          f"{'first':<11}{'last':<11}")
    for _, r in helpful.head(40).iterrows():
        food_short = r["food"][:58]
        print(f"  {food_short:<60} {int(r['n']):>5} {r['binge_share']*100:>6.1f}% "
              f"{r['lift']:>5.2f} {int(r['n_years']):>6} {r['span_days']/365:>7.1f} "
              f"{r['first_date'].strftime('%Y-%m'):<8}    {r['last_date'].strftime('%Y-%m'):<8}")

    # ── Span: foods with most years (cross-era confirmation) ──
    print("\n" + "=" * 110)
    print("FOODS WITH MOST YEARS (cross-era staples), sorted by lift")
    print("=" * 110)
    staples = long_term[long_term["n_years"] >= 6].sort_values("lift", ascending=False)
    print(f"  {'food':<60} {'n':>5} {'binge%':>7} {'lift':>5} {'years':>6}")
    for _, r in staples.iterrows():
        food_short = r["food"][:58]
        print(f"  {food_short:<60} {int(r['n']):>5} {r['binge_share']*100:>6.1f}% "
              f"{r['lift']:>5.2f} {int(r['n_years']):>6}")

    # Save artifact
    out_path = ROOT / "analysis" / "DA2_binge_foods.csv"
    long_term.sort_values("lift", ascending=False).to_csv(out_path, index=False)
    print(f"\nSaved {out_path}")


if __name__ == "__main__":
    main()
