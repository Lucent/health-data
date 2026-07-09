"""DA5. Foods that exert downward (or upward) pressure on the day's total intake.

For each fuzzy cluster (DA4 normalization):
  - mean intake on days consumed                      (raw)
  - mean intake on same-era days NOT consumed         (within ±30d of each food day)
  - within-era residual: mean (intake_day − 30d centered rolling mean)
       Negative = consuming this food correlates with below-era-average intake
       Positive = consuming this food correlates with above-era-average intake

The within-era residual is the cleanest "pressure" metric: it controls for slow
era drift in baseline intake, so a food that always appears in low-intake phases
doesn't get credit just for that.

Sort by ascending residual (most-downward-pressure first).
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
INTAKE_MULTIPLIER = 1.09


def food_key(food):
    """Same fuzzy normalizer as DA4."""
    s = str(food).lower().strip()
    if "," in s:
        s = s.rsplit(",", 1)[0].strip()
    s = s.replace("'", "").replace("’", "")
    s = re.sub(r"\([^)]*\)", "", s)
    s = re.sub(r"\b(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)[a-z]*\s*\d{2,4}", "", s)
    s = re.sub(r"\b(old|new|corrected|approx|limited\s+edition|updated|original|"
                r"per\s+bag|per\s+container|v\d+|correct|chickn|chikn|"
                r"chick\s*n|chik\s*n)\b", "", s)
    s = re.sub(r"\b(ultra[-\s]?filtered|ultra[-\s]?pasteurized|"
                r"grass[-\s]?fed|cage[-\s]?free|extra\s+virgin|"
                r"unsweetened|low[-\s]?fat|reduced[-\s]?fat|fat[-\s]?free)\b", "", s)
    s = re.sub(r"morning\s*star", "morningstar", s)
    s = re.sub(r"trader\s+joe[s]?", "trader joes", s)
    s = re.sub(r"boar[s]?\s*head", "boars head", s)
    s = re.sub(r"ben\s*&?\s*jerry[s]?", "ben jerrys", s)
    s = re.sub(r"hannaford\s*-\s*", "", s)
    s = re.sub(r"[*~]", "", s)
    s = re.sub(r"^\s*-\s*", "", s)
    s = re.sub(r"\s*-\s*$", "", s)
    s = re.sub(r"\s+", " ", s).strip()
    if re.search(r"\bslate\b", s) and "milk" in s and ("chocolate" in s or "dark" in s):
        return "slate chocolate milk"
    if re.search(r"^(?:.*-\s*)?heavy\s+cream\s*$", s):
        return "heavy cream"
    if "maple syrup" in s and "choczero" not in s and "sugar-free" not in s and "sugar free" not in s and "honest" not in s:
        return "maple syrup"
    if re.search(r"^butter[\s,]*unsalted\s*$", s) or re.search(r"^unsalted\s+butter\s*$", s):
        return "butter unsalted"
    if re.search(r"^butter\s*$", s):
        return "butter"
    if re.search(r"^olive\s+oil\s*$", s):
        return "olive oil"
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

    pre = df[(df["effective_level"] == 0) & (df["year"] >= 2014)].copy().reset_index(drop=True)
    print(f"Days: {len(pre)} ({pre['date'].min().date()} → {pre['date'].max().date()})")
    print(f"Mean intake: {pre['calories'].mean():.0f} cal/day, "
          f"median: {pre['calories'].median():.0f}")

    # Within-era baseline: 30-day centered rolling mean intake
    pre["intake_30dc"] = pre["calories"].rolling(30, min_periods=15, center=True).mean()
    pre["intake_resid_30d"] = pre["calories"] - pre["intake_30dc"]
    print(f"30-d centered rolling mean range: {pre['intake_30dc'].min():.0f} - "
          f"{pre['intake_30dc'].max():.0f}")

    # Food log — rename per-item calories to avoid conflict with daily total
    foods = foods.rename(columns={"calories": "item_calories"})
    pre_merge = pre[["date", "binge", "year", "calories", "intake_resid_30d"]].rename(
        columns={"calories": "day_calories"})
    foods = foods.merge(pre_merge, on="date", how="inner")
    foods["fuzzy_key"] = foods["food"].apply(food_key)
    foods_dedupe = foods.drop_duplicates(subset=["date", "fuzzy_key"])

    # Per-cluster aggregates
    grp = foods_dedupe.groupby("fuzzy_key")
    agg = grp.agg(
        n=("binge", "size"),
        binge_count=("binge", "sum"),
        first_date=("date", "min"),
        last_date=("date", "max"),
        n_years=("year", "nunique"),
        mean_intake=("day_calories", "mean"),
        mean_resid=("intake_resid_30d", "mean"),
    ).reset_index()
    agg["binge_share"] = agg["binge_count"] / agg["n"]
    agg["span_days"] = (agg["last_date"] - agg["first_date"]).dt.days
    overall_binge = pre["binge"].mean()
    agg["lift"] = agg["binge_share"] / overall_binge

    long_term = agg[(agg["n"] >= 50) & (agg["n_years"] >= 3) & (agg["span_days"] >= 700)].copy()
    print(f"\nClusters passing filter (n≥50, ≥3 years, span≥700d): {len(long_term)}")

    # ── Top DOWNWARD-PRESSURE foods (most negative residual) ──
    print("\n" + "=" * 120)
    print("TOP DOWNWARD-PRESSURE FOODS (most negative within-era residual)")
    print("Days with this food average X cal LOWER than the surrounding 30-day window")
    print("=" * 120)
    print(f"  {'cluster':<55} {'n':>5} {'mean_intake':>11} {'resid':>8} "
          f"{'binge%':>7} {'years':>6}")
    down = long_term.sort_values("mean_resid").head(30)
    for _, r in down.iterrows():
        print(f"  {r['fuzzy_key'][:55]:<55} {int(r['n']):>5} "
              f"{r['mean_intake']:>11.0f} {r['mean_resid']:>+8.0f} "
              f"{r['binge_share']*100:>6.1f}% {int(r['n_years']):>6}")

    # ── Top UPWARD-PRESSURE foods (most positive residual) ──
    print("\n" + "=" * 120)
    print("TOP UPWARD-PRESSURE FOODS (most positive within-era residual)")
    print("Days with this food average X cal HIGHER than the surrounding 30-day window")
    print("=" * 120)
    print(f"  {'cluster':<55} {'n':>5} {'mean_intake':>11} {'resid':>8} "
          f"{'binge%':>7} {'years':>6}")
    up = long_term.sort_values("mean_resid", ascending=False).head(30)
    for _, r in up.iterrows():
        print(f"  {r['fuzzy_key'][:55]:<55} {int(r['n']):>5} "
              f"{r['mean_intake']:>11.0f} {r['mean_resid']:>+8.0f} "
              f"{r['binge_share']*100:>6.1f}% {int(r['n_years']):>6}")

    # ── Compare residual vs binge-rate as ranking signal ──
    print("\n" + "=" * 120)
    print("FOODS WHERE INTAKE-RESIDUAL AND BINGE-RATE AGREE (both pull down OR both pull up)")
    print("=" * 120)
    long_term["binge_dir"] = np.sign(long_term["lift"] - 1.0)
    long_term["resid_dir"] = np.sign(long_term["mean_resid"])
    agree = long_term[long_term["binge_dir"] == long_term["resid_dir"]]
    print(f"  Both pull DOWN (resid<0 AND lift<1):")
    for _, r in agree[(agree["resid_dir"] < 0)].sort_values("mean_resid").head(20).iterrows():
        print(f"    {r['fuzzy_key'][:55]:<55} resid={r['mean_resid']:>+5.0f}  "
              f"lift={r['lift']:.2f}  binge%={r['binge_share']*100:.1f}%  n={int(r['n'])}")
    print(f"\n  Both pull UP (resid>0 AND lift>1):")
    for _, r in agree[(agree["resid_dir"] > 0)].sort_values("mean_resid", ascending=False).head(20).iterrows():
        print(f"    {r['fuzzy_key'][:55]:<55} resid={r['mean_resid']:>+5.0f}  "
              f"lift={r['lift']:.2f}  binge%={r['binge_share']*100:.1f}%  n={int(r['n'])}")

    # ── Foods where signals DISAGREE (interesting cases) ──
    print("\n" + "=" * 120)
    print("INTERESTING DISAGREEMENTS")
    print("(food doesn't trigger binges but raises intake, or doesn't help binges but lowers intake)")
    print("=" * 120)
    disagree = long_term[long_term["binge_dir"] != long_term["resid_dir"]]
    print(f"\n  Lowers intake but DOESN'T reduce binges (resid<0, lift>1):")
    for _, r in disagree[(disagree["resid_dir"] < 0) & (disagree["binge_dir"] > 0)].sort_values("mean_resid").iterrows():
        print(f"    {r['fuzzy_key'][:55]:<55} resid={r['mean_resid']:>+5.0f}  "
              f"lift={r['lift']:.2f}  binge%={r['binge_share']*100:.1f}%  n={int(r['n'])}")
    print(f"\n  Raises intake but DOESN'T cause binges (resid>0, lift<1):")
    for _, r in disagree[(disagree["resid_dir"] > 0) & (disagree["binge_dir"] < 0)].sort_values("mean_resid", ascending=False).head(15).iterrows():
        print(f"    {r['fuzzy_key'][:55]:<55} resid={r['mean_resid']:>+5.0f}  "
              f"lift={r['lift']:.2f}  binge%={r['binge_share']*100:.1f}%  n={int(r['n'])}")

    out_path = ROOT / "analysis" / "DA5_intake_pressure.csv"
    long_term.sort_values("mean_resid").to_csv(out_path, index=False)
    print(f"\nSaved {out_path}")


if __name__ == "__main__":
    main()
