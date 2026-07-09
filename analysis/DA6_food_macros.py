"""DA6. What's the actual difference between trigger and protective food clusters?

Same fuzzy clustering as DA4/DA5. For each cluster, compute mean macro composition
per food log entry (the serving size as logged), then test which macro features
distinguish triggers from protectors.

Macro features per cluster:
  - fat % of calories
  - carb % of calories
  - protein % of calories
  - sugar grams / 100 cal
  - fiber grams / 100 cal
  - sodium mg / 100 cal
  - sugar fraction of carbs
  - per-serving cal density

Filter relaxed to n≥30, ≥2 years, span≥400d to get a larger sample for
correlation analysis.
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
INTAKE_MULTIPLIER = 1.09


def food_key(food):
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
    df["intake_30dc"] = df["calories"].rolling(30, min_periods=15, center=True).mean()
    df["intake_resid_30d"] = df["calories"] - df["intake_30dc"]

    pre = df[(df["effective_level"] == 0) & (df["year"] >= 2014)].copy()
    foods_in = foods.copy()
    foods_in = foods_in.rename(columns={"calories": "item_calories"})
    foods_in = foods_in.merge(
        pre[["date", "binge", "year", "intake_resid_30d"]], on="date", how="inner")
    foods_in["fuzzy_key"] = foods_in["food"].apply(food_key)

    # Per-cluster macro means (averaged across each food log entry)
    macro_grp = foods_in.groupby("fuzzy_key")
    macros = macro_grp.agg(
        n_entries=("food", "size"),
        cal=("item_calories", "mean"),
        carbs=("carbs_g", "mean"),
        fat=("fat_g", "mean"),
        protein=("protein_g", "mean"),
        sugar=("sugars_g", "mean"),
        fiber=("fiber_g", "mean"),
        sodium=("sodium_mg", "mean"),
    ).reset_index()
    # Density features
    macros["fat_pct"] = macros["fat"] * 9 / macros["cal"]
    macros["carb_pct"] = macros["carbs"] * 4 / macros["cal"]
    macros["protein_pct"] = macros["protein"] * 4 / macros["cal"]
    macros["sugar_per_100cal"] = macros["sugar"] * 100 / macros["cal"]
    macros["fiber_per_100cal"] = macros["fiber"] * 100 / macros["cal"]
    macros["sodium_per_100cal"] = macros["sodium"] * 100 / macros["cal"]
    macros["sugar_frac_of_carb"] = macros["sugar"] / macros["carbs"].replace(0, np.nan)

    # Per-cluster daily aggregation (deduped by date)
    foods_dedupe = foods_in.drop_duplicates(subset=["date", "fuzzy_key"])
    daily_grp = foods_dedupe.groupby("fuzzy_key")
    daily = daily_grp.agg(
        n=("binge", "size"),
        binge_count=("binge", "sum"),
        first_date=("date", "min"),
        last_date=("date", "max"),
        n_years=("year", "nunique"),
        mean_resid=("intake_resid_30d", "mean"),
    ).reset_index()
    daily["binge_share"] = daily["binge_count"] / daily["n"]
    daily["span_days"] = (daily["last_date"] - daily["first_date"]).dt.days
    overall_binge = pre["binge"].mean()
    daily["lift"] = daily["binge_share"] / overall_binge

    # Merge
    full = daily.merge(macros, on="fuzzy_key", how="left")

    # Filter: relaxed to get more clusters
    long_term = full[(full["n"] >= 30) & (full["n_years"] >= 2) & (full["span_days"] >= 400)].copy()
    print(f"Clusters (n≥30, ≥2 years, span≥400d): {len(long_term)}")

    # ── Trigger vs protective summary by quadrant ──
    long_term["category"] = "neutral"
    long_term.loc[(long_term["lift"] > 1.3) & (long_term["mean_resid"] > 50),
                    "category"] = "trigger"
    long_term.loc[(long_term["lift"] < 0.7) & (long_term["mean_resid"] < -30),
                    "category"] = "protective"
    long_term.loc[(long_term["lift"] > 1.3) & (long_term["mean_resid"] <= 50),
                    "category"] = "binge_only"
    long_term.loc[(long_term["lift"] < 0.7) & (long_term["mean_resid"] >= -30),
                    "category"] = "calm_high_intake"

    print("\nCategory counts:")
    print(long_term["category"].value_counts())

    # ── Mean macros per category ──
    print("\n" + "=" * 110)
    print("MEAN MACROS BY CATEGORY")
    print("=" * 110)
    cols_show = ["cal", "fat_pct", "carb_pct", "protein_pct",
                  "sugar_per_100cal", "fiber_per_100cal", "sodium_per_100cal"]
    by_cat = long_term.groupby("category")[cols_show].mean()
    by_cat["n_clusters"] = long_term.groupby("category").size()
    print(by_cat.round(3))

    # ── Correlations: macros vs lift, macros vs residual ──
    print("\n" + "=" * 110)
    print("CORRELATIONS BETWEEN MACRO FEATURES AND BINGE LIFT / INTAKE RESIDUAL")
    print("=" * 110)
    for c in ["cal", "fat_pct", "carb_pct", "protein_pct",
               "sugar_per_100cal", "fiber_per_100cal", "sodium_per_100cal",
               "sugar_frac_of_carb"]:
        v = long_term.dropna(subset=[c])
        r_lift = float(np.corrcoef(v[c], v["lift"])[0, 1])
        r_resid = float(np.corrcoef(v[c], v["mean_resid"])[0, 1])
        print(f"  {c:<25}  r(lift) = {r_lift:+.3f}   r(resid) = {r_resid:+.3f}   n={len(v)}")

    # ── Show triggers and protective with macros ──
    print("\n" + "=" * 130)
    print("TRIGGER CLUSTERS (lift > 1.5) — what are their macro profiles?")
    print("=" * 130)
    print(f"  {'cluster':<48} {'n':>4} {'cal':>5} {'fat%':>5} {'carb%':>6} "
          f"{'prot%':>6} {'sug/100':>7} {'fib/100':>7} {'Na/100':>7} {'lift':>5}")
    triggers = long_term[long_term["lift"] > 1.5].sort_values("lift", ascending=False)
    for _, r in triggers.iterrows():
        print(f"  {r['fuzzy_key'][:48]:<48} {int(r['n']):>4} {r['cal']:>5.0f} "
              f"{r['fat_pct']*100:>5.0f} {r['carb_pct']*100:>6.0f} "
              f"{r['protein_pct']*100:>6.0f} {r['sugar_per_100cal']:>7.2f} "
              f"{r['fiber_per_100cal']:>7.2f} {r['sodium_per_100cal']:>7.0f} "
              f"{r['lift']:>5.2f}")

    print("\n" + "=" * 130)
    print("PROTECTIVE CLUSTERS (lift < 0.5)")
    print("=" * 130)
    print(f"  {'cluster':<48} {'n':>4} {'cal':>5} {'fat%':>5} {'carb%':>6} "
          f"{'prot%':>6} {'sug/100':>7} {'fib/100':>7} {'Na/100':>7} {'lift':>5}")
    prot = long_term[long_term["lift"] < 0.5].sort_values("lift")
    for _, r in prot.head(40).iterrows():
        print(f"  {r['fuzzy_key'][:48]:<48} {int(r['n']):>4} {r['cal']:>5.0f} "
              f"{r['fat_pct']*100:>5.0f} {r['carb_pct']*100:>6.0f} "
              f"{r['protein_pct']*100:>6.0f} {r['sugar_per_100cal']:>7.2f} "
              f"{r['fiber_per_100cal']:>7.2f} {r['sodium_per_100cal']:>7.0f} "
              f"{r['lift']:>5.2f}")

    # ── Plot ──
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    cat_colors = {"trigger": "tab:red", "protective": "tab:green",
                   "binge_only": "tab:orange", "calm_high_intake": "tab:blue",
                   "neutral": "lightgray"}

    ax = axes[0, 0]
    for cat, c in cat_colors.items():
        sub = long_term[long_term["category"] == cat]
        ax.scatter(sub["lift"], sub["mean_resid"], c=c, label=cat, alpha=0.7,
                    s=40, edgecolor="black", linewidth=0.5)
    ax.axvline(1.0, color="black", lw=0.5)
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlabel("binge lift (binge_rate / 0.077)")
    ax.set_ylabel("intake residual (cal/d, vs 30d mean)")
    ax.set_title("Lift vs intake residual — color by category")
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_xscale("log")

    ax = axes[0, 1]
    sc = ax.scatter(long_term["fat_pct"] * 100, long_term["lift"],
                     c=long_term["sugar_per_100cal"], cmap="viridis",
                     s=40, alpha=0.7, edgecolor="black", linewidth=0.5)
    plt.colorbar(sc, ax=ax, label="sugar g/100 cal")
    ax.set_xlabel("fat % of calories")
    ax.set_ylabel("binge lift")
    ax.set_title("Lift vs fat fraction (color = sugar density)")
    ax.set_yscale("log")
    ax.grid(alpha=0.3)

    ax = axes[1, 0]
    sc = ax.scatter(long_term["protein_pct"] * 100, long_term["lift"],
                     c=long_term["fiber_per_100cal"], cmap="viridis",
                     s=40, alpha=0.7, edgecolor="black", linewidth=0.5)
    plt.colorbar(sc, ax=ax, label="fiber g/100 cal")
    ax.set_xlabel("protein % of calories")
    ax.set_ylabel("binge lift")
    ax.set_title("Lift vs protein fraction (color = fiber density)")
    ax.set_yscale("log")
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    sc = ax.scatter(long_term["sodium_per_100cal"], long_term["mean_resid"],
                     c=long_term["fat_pct"] * 100, cmap="viridis",
                     s=40, alpha=0.7, edgecolor="black", linewidth=0.5)
    plt.colorbar(sc, ax=ax, label="fat %")
    ax.set_xlabel("sodium mg / 100 cal")
    ax.set_ylabel("intake residual")
    ax.set_title("Sodium vs intake residual (color = fat %)")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    plot_path = ROOT / "analysis" / "DA6_food_macros.png"
    fig.savefig(plot_path, dpi=140)
    print(f"\nSaved {plot_path}")

    long_term.to_csv(ROOT / "analysis" / "DA6_food_macros.csv", index=False)


if __name__ == "__main__":
    main()
