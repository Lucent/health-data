"""DA. Descriptive: what does binge frequency / excess intake correlate with?

Step back from set-point modelling. Just look at binges and excess as
descriptive objects across the 15 years.

Tests:
  1. Definition sweep: how does binge rate change with the threshold?
  2. Day-of-week pattern
  3. Month / season pattern
  4. Year-over-year and diet-epoch
  5. Inter-arrival distribution: are binges Poisson or clustered?
  6. AR structure: does binge yesterday predict binge today?
  7. Excess intake (continuous) vs steps, sleep, temperature, day-of-week,
     after-restriction (with phase fixed effects to remove era-level confound)
  8. What foods appear disproportionately on binge days
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
INTAKE_MULTIPLIER = 1.09


def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    intake["calories"] = intake["calories"] * INTAKE_MULTIPLIER
    kf = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    epochs = pd.read_csv(ROOT / "intake" / "diet_epochs.csv", parse_dates=["start", "end"])
    steps_df = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])
    sleep_df = pd.read_csv(ROOT / "steps-sleep" / "sleep.csv", parse_dates=["date"])

    df = intake[["date", "calories"]].merge(
        kf[["date", "fat_mass_lbs", "tdee"]], on="date")
    df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df = df.merge(steps_df[["date", "steps"]], on="date", how="left")
    df = df.merge(sleep_df[["date", "sleep_hours"]], on="date", how="left")
    df = df.sort_values("date").reset_index(drop=True)
    df["surplus"] = df["calories"] - df["tdee"]
    df["weekday"] = df["date"].dt.dayofweek  # Monday=0
    df["month"] = df["date"].dt.month
    df["year"] = df["date"].dt.year

    pre = df[df["effective_level"] == 0].copy().reset_index(drop=True)
    print(f"Pre-tirz days: {len(pre)}")

    # ── 1. Definition sweep ──
    print("\n" + "=" * 80)
    print("1. BINGE DEFINITION SWEEP")
    print("=" * 80)
    print(f"  {'threshold':>12} {'n':>5} {'rate':>6}")
    for t in [500, 750, 1000, 1500, 2000, 2500, 3000]:
        n = int((pre["surplus"] > t).sum())
        print(f"  >TDEE+{t:>5} cal  {n:>5} {n/len(pre)*100:>5.1f}%")
    pre["binge"] = (pre["surplus"] > 1000).astype(int)
    print(f"\nUsing > TDEE+1000 cal (AG convention): {pre['binge'].sum()} binges, "
          f"{pre['binge'].mean()*100:.1f}% rate")

    # ── 2. Day-of-week ──
    print("\n" + "=" * 80)
    print("2. DAY-OF-WEEK PATTERN")
    print("=" * 80)
    print(f"  {'day':<10} {'n':>4} {'binge%':>7} {'mean surp':>10} {'median surp':>11}")
    days = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]
    for i, d in enumerate(days):
        g = pre[pre["weekday"] == i]
        print(f"  {d:<10} {len(g):>4} {g['binge'].mean()*100:>6.1f}% "
              f"{g['surplus'].mean():>+10.0f} {g['surplus'].median():>+11.0f}")

    # ── 3. Month / season ──
    print("\n" + "=" * 80)
    print("3. MONTH / SEASON")
    print("=" * 80)
    months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    print(f"  {'month':<6} {'n':>5} {'binge%':>7} {'mean surp':>10}")
    for i in range(1, 13):
        g = pre[pre["month"] == i]
        print(f"  {months[i-1]:<6} {len(g):>5} {g['binge'].mean()*100:>6.1f}% "
              f"{g['surplus'].mean():>+10.0f}")

    # ── 4. Year-over-year + epochs ──
    print("\n" + "=" * 80)
    print("4. YEAR & EPOCH")
    print("=" * 80)
    print(f"  {'year':<6} {'n':>5} {'binge%':>7} {'mean surp':>10} {'mean intake':>12}")
    for y in sorted(pre["year"].unique()):
        g = pre[pre["year"] == y]
        if len(g) < 50:
            continue
        print(f"  {y:<6} {len(g):>5} {g['binge'].mean()*100:>6.1f}% "
              f"{g['surplus'].mean():>+10.0f} {g['calories'].mean():>12.0f}")

    print("\nBy diet epoch (top 15 by binge rate):")
    epoch_results = []
    for _, e in epochs.iterrows():
        m = (pre["date"] >= e["start"]) & (pre["date"] <= e["end"])
        g = pre[m]
        if len(g) < 14:
            continue
        epoch_results.append({
            "label": e["label"],
            "start": e["start"],
            "n_days": len(g),
            "binge_rate": g["binge"].mean(),
            "mean_intake": g["calories"].mean(),
            "mean_surp": g["surplus"].mean(),
        })
    er = pd.DataFrame(epoch_results).sort_values("binge_rate", ascending=False)
    print(f"  {'epoch':<28} {'start':<12} {'n':>5} {'binge%':>7} "
          f"{'intake':>7} {'surp':>7}")
    for _, r in er.head(15).iterrows():
        print(f"  {r['label']:<28} {r['start'].strftime('%Y-%m-%d'):<12} "
              f"{int(r['n_days']):>5} {r['binge_rate']*100:>6.1f}% "
              f"{r['mean_intake']:>7.0f} {r['mean_surp']:>+7.0f}")

    # ── 5. Inter-arrival ──
    print("\n" + "=" * 80)
    print("5. INTER-ARRIVAL TIMES BETWEEN BINGES")
    print("=" * 80)
    binge_dates = pre[pre["binge"] == 1]["date"].values
    if len(binge_dates) > 1:
        gaps = np.diff(binge_dates).astype("timedelta64[D]").astype(int)
        print(f"  n binges: {len(binge_dates)}")
        print(f"  Mean gap: {gaps.mean():.1f} days")
        print(f"  Median gap: {np.median(gaps):.0f} days")
        print(f"  IQR: {np.percentile(gaps, 25):.0f} - {np.percentile(gaps, 75):.0f}")
        print(f"  % of gaps == 1 (consecutive day binges): {np.mean(gaps == 1)*100:.0f}%")
        print(f"  % of gaps ≤ 3 days: {np.mean(gaps <= 3)*100:.0f}%")
        print(f"  % of gaps ≥ 30 days: {np.mean(gaps >= 30)*100:.0f}%")
        # Compare to Poisson expectation
        lam = 1 / gaps.mean()  # rate per day
        # Under Poisson, gap distribution is exponential with mean 1/λ
        p_gap_1 = 1 - np.exp(-lam)
        print(f"\n  Under Poisson with same mean: P(gap=1) = {p_gap_1*100:.1f}%")
        print(f"  Observed P(gap=1): {np.mean(gaps == 1)*100:.1f}% "
              f"({'CLUSTERED' if np.mean(gaps == 1) > p_gap_1 + 0.05 else 'POISSON-LIKE'})")

    # ── 6. AR structure ──
    print("\n" + "=" * 80)
    print("6. BINGE AUTOCORRELATION")
    print("=" * 80)
    for lag in [1, 2, 3, 7, 14, 30]:
        # Probability of binge on day t given binge on day t-lag
        binge_lag = pre["binge"].shift(lag)
        v = ~binge_lag.isna()
        p_yes = pre.loc[v & (binge_lag == 1), "binge"].mean()
        p_no = pre.loc[v & (binge_lag == 0), "binge"].mean()
        print(f"  lag={lag:>2}d:  P(binge | binge {lag}d ago) = {p_yes*100:.1f}%, "
              f"P(binge | no binge {lag}d ago) = {p_no*100:.1f}%, "
              f"ratio = {p_yes/p_no if p_no > 0 else float('nan'):.2f}")

    # ── 7. Daily covariate correlations (with phase FE) ──
    print("\n" + "=" * 80)
    print("7. DAILY COVARIATES vs BINGE PROBABILITY (with epoch fixed effects)")
    print("=" * 80)
    pre["epoch_id"] = -1
    for ei, (_, e) in enumerate(epochs.iterrows()):
        m = (pre["date"] >= e["start"]) & (pre["date"] <= e["end"])
        pre.loc[m, "epoch_id"] = ei

    # Simple within-epoch correlations
    valid = pre.dropna(subset=["binge", "steps", "sleep_hours"])
    print(f"\n  Marginal correlations (no controls):")
    for col in ["steps", "sleep_hours", "weekday"]:
        if col in valid.columns and valid[col].notna().any():
            v = valid.dropna(subset=[col])
            r = float(np.corrcoef(v[col], v["binge"])[0, 1])
            print(f"    binge vs {col}: r = {r:+.4f}  (n={len(v)})")

    # By weekday with confidence
    print(f"\n  Mean ± SE binge rate by weekday:")
    for i, d in enumerate(days):
        g = valid[valid["weekday"] == i]
        if len(g) > 100:
            p = g["binge"].mean()
            se = np.sqrt(p * (1 - p) / len(g))
            print(f"    {d}: {p*100:.2f}% ± {se*100:.2f}%  (n={len(g)})")

    # Lagged restriction → binge
    print(f"\n  After-restriction effect:")
    pre["surplus_7d"] = pre["surplus"].rolling(7, min_periods=4).mean().shift(1)
    pre["surplus_14d"] = pre["surplus"].rolling(14, min_periods=7).mean().shift(1)
    pre["surplus_30d"] = pre["surplus"].rolling(30, min_periods=15).mean().shift(1)
    for col in ["surplus_7d", "surplus_14d", "surplus_30d"]:
        v = pre.dropna(subset=[col, "binge"])
        # Bin into quintiles
        v["bin"] = pd.qcut(v[col], 5, labels=False, duplicates="drop")
        print(f"\n    {col} quintiles → binge rate:")
        for b in sorted(v["bin"].dropna().unique()):
            g = v[v["bin"] == b]
            print(f"      Q{int(b)+1} (mean {col}={g[col].mean():>+6.0f}): "
                  f"binge rate = {g['binge'].mean()*100:.2f}% (n={len(g)})")

    # ── 8. Food correlates ──
    print("\n" + "=" * 80)
    print("8. WHAT FOODS APPEAR ON BINGE DAYS DISPROPORTIONATELY?")
    print("=" * 80)
    foods = pd.read_csv(ROOT / "intake" / "intake_foods.csv", parse_dates=["date"])
    foods = foods.merge(pre[["date", "binge"]], on="date", how="inner")
    food_binge = (foods.groupby("food")["binge"]
                   .agg(["mean", "count"])
                   .rename(columns={"mean": "binge_share", "count": "n"})
                   .query("n >= 30"))
    overall = pre["binge"].mean()
    food_binge["lift"] = food_binge["binge_share"] / overall
    print(f"\nOverall binge rate: {overall*100:.1f}%")
    print(f"\nTop 25 foods most associated with binge days (lift > 1.5, n ≥ 30):")
    top = food_binge.sort_values("lift", ascending=False).head(25)
    for f, r in top.iterrows():
        print(f"  {f[:50]:<52}  binge_share={r['binge_share']*100:>5.1f}%  "
              f"lift={r['lift']:>4.2f}  n={int(r['n'])}")

    print(f"\nTop 25 foods LEAST associated with binge days (lift < 0.5):")
    bot = food_binge[food_binge["lift"] < 0.5].sort_values("lift").head(25)
    for f, r in bot.iterrows():
        print(f"  {f[:50]:<52}  binge_share={r['binge_share']*100:>5.1f}%  "
              f"lift={r['lift']:>4.2f}  n={int(r['n'])}")

    # ── Plot ──
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    ax = axes[0, 0]
    by_year = pre.groupby("year")["binge"].mean() * 100
    ax.bar(by_year.index, by_year.values)
    ax.set_xlabel("Year")
    ax.set_ylabel("Binge rate (%)")
    ax.set_title("Binge rate by year (pre-tirz)")
    ax.grid(alpha=0.3)

    ax = axes[0, 1]
    by_dow = pre.groupby("weekday")["binge"].mean() * 100
    ax.bar(range(7), by_dow.values)
    ax.set_xticks(range(7))
    ax.set_xticklabels(days)
    ax.set_ylabel("Binge rate (%)")
    ax.set_title("Binge rate by day of week")
    ax.grid(alpha=0.3)

    ax = axes[1, 0]
    by_month = pre.groupby("month")["binge"].mean() * 100
    ax.bar(range(1, 13), by_month.values)
    ax.set_xticks(range(1, 13))
    ax.set_xticklabels(months)
    ax.set_ylabel("Binge rate (%)")
    ax.set_title("Binge rate by month")
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    if len(binge_dates) > 1:
        ax.hist(np.clip(gaps, 0, 50), bins=50, alpha=0.7)
        ax.set_xlabel("days between consecutive binges")
        ax.set_ylabel("count")
        ax.set_title(f"Inter-arrival distribution (clipped at 50d)\nmedian = {np.median(gaps):.0f}d, mean = {gaps.mean():.1f}d")
        ax.grid(alpha=0.3)

    fig.tight_layout()
    plot_path = ROOT / "analysis" / "DA_binge_correlates.png"
    fig.savefig(plot_path, dpi=140)
    print(f"\nSaved {plot_path}")


if __name__ == "__main__":
    main()
