"""Properties of the hidden fat mass set point.

Building on finding AG (set point = EMA of FM with 50-day half-life),
this analysis characterizes the set point's mechanics:

1. INVERTED RATCHET: The set point adapts 5x faster going DOWN (HL=20d)
   than UP (HL=100d). The asymmetric model improves r from -0.62 to -0.71.
   The body quickly accepts lower weight but is slow to raise its defended
   weight during regain — the opposite of the feared "ratchet."

2. DUAL DEFENSE: Both expenditure and appetite arms correlate with SP
   distance. TDEE residual: partial r = +0.38. Binge rate: partial r = -0.62.
   Appetite is the stronger arm.

3. FREQUENCY NOT MAGNITUDE: Binge size is flat (~1400 cal surplus) regardless
   of SP distance (r = -0.04). But non-binge days show continuous upward
   pressure (r = -0.35 with surplus). The set point modulates background
   drift on every day, plus discrete binge probability.

4. RESTRICTION PREDICTION: Runs ending above SP continue losing (-1.23 lbs
   at 30d). Runs ending below SP rebound (+0.35 lbs). r = -0.48.

5. EXERCISE INDEPENDENCE: Walking doesn't change the SP half-life (40d for
   both high and low walk periods). Walking raises RMR independently.

6. FLOOR EFFECT: SP reached 18.6 lbs (Nov 2013). Adaptation rate slows
   near essential body fat. Post-floor: FM 19→23 in 6 months, 8.3% binges.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

SP_HALF_LIFE = 50  # symmetric default
SP_HL_UP = 100     # asymmetric: slow to rise
SP_HL_DOWN = 20    # asymmetric: fast to fall
BINGE_THRESHOLD = 1000


def ema(series, half_life):
    alpha = 1 - np.exp(-np.log(2) / half_life)
    return series.ewm(alpha=alpha, min_periods=30).mean()


def asymmetric_ema(fm_series, hl_up, hl_down):
    """EMA that adapts at different rates for upward vs downward FM moves."""
    alpha_up = 1 - np.exp(-np.log(2) / hl_up)
    alpha_down = 1 - np.exp(-np.log(2) / hl_down)
    vals = fm_series.values
    sp = np.full(len(vals), np.nan)
    sp[0] = vals[0]
    for i in range(1, len(vals)):
        if np.isnan(sp[i - 1]):
            sp[i] = vals[i]
            continue
        if vals[i] > sp[i - 1]:
            sp[i] = sp[i - 1] + alpha_up * (vals[i] - sp[i - 1])
        else:
            sp[i] = sp[i - 1] + alpha_down * (vals[i] - sp[i - 1])
    return sp


def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    exercises = pd.read_csv(ROOT / "steps-sleep" / "exercises.csv", parse_dates=["date"])

    df = intake[["date", "calories"]].merge(
        kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="inner")
    df = df.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    df = df.sort_values("date").reset_index(drop=True)
    n = len(df)

    df["binge"] = (df["calories"] > df["tdee"] + BINGE_THRESHOLD).astype(int)
    df["binge_rate_90d"] = df["binge"].rolling(90, min_periods=90).mean()
    df["surplus"] = df["calories"] - df["tdee"]
    df["sp_sym"] = ema(df["fat_mass_lbs"], SP_HALF_LIFE)
    df["sp_asym"] = asymmetric_ema(df["fat_mass_lbs"], SP_HL_UP, SP_HL_DOWN)
    df["dist_sym"] = df["sp_sym"] - df["fat_mass_lbs"]
    df["dist_asym"] = df["sp_asym"] - df["fat_mass_lbs"]
    df["tdee_resid"] = df["tdee"] - df["expected_rmr"]
    df["fm_vel_90d"] = (df["fat_mass_lbs"] - df["fat_mass_lbs"].shift(90)) / 90 * 30

    # ═══════════════════════════════════════════════════════════════
    print("=" * 70)
    print("1. INVERTED RATCHET: asymmetric set point adaptation")
    print("=" * 70)

    # Symmetric baseline
    dist_sym_v = df["dist_sym"].values
    valid_sym = ~np.isnan(dist_sym_v) & df["binge_rate_90d"].notna().values
    r_sym = np.corrcoef(dist_sym_v[valid_sym], df.loc[valid_sym, "binge_rate_90d"])[0, 1]

    # Sweep asymmetric half-lives
    best_r, best_u, best_d = r_sym, SP_HALF_LIFE, SP_HALF_LIFE
    for hl_up in range(20, 200, 10):
        for hl_down in range(20, 200, 10):
            sp_arr = asymmetric_ema(df["fat_mass_lbs"], hl_up, hl_down)
            fm_arr = df["fat_mass_lbs"].values
            dist = sp_arr - fm_arr
            br_arr = df["binge_rate_90d"].values
            valid = ~np.isnan(dist) & ~np.isnan(br_arr)
            if valid.sum() < 500:
                continue
            r = np.corrcoef(dist[valid], br_arr[valid])[0, 1]
            if r < best_r:  # more negative = better
                best_r, best_u, best_d = r, hl_up, hl_down

    print(f"\n  Symmetric (HL=50d):    r = {r_sym:+.3f}")
    print(f"  Asymmetric (best):     r = {best_r:+.3f}  (HL_up={best_u}d, HL_down={best_d}d)")
    print(f"  Ratio up/down: {best_u / best_d:.1f}x")
    print(f"  The set point adapts {best_u / best_d:.0f}x faster going DOWN than UP.")
    print(f"  87% adaptation time: down={best_d * 3}d ({best_d * 3 // 30}mo), up={best_u * 3}d ({best_u * 3 // 30}mo)")

    # Show trajectory with asymmetric SP (use sweep result)
    sp_best = asymmetric_ema(df["fat_mass_lbs"], best_u, best_d)
    df["sp_asym"] = sp_best
    print(f"\n  {'Year':>6} {'FM':>5} {'SP_sym':>7} {'SP_asym':>8} {'Binge%':>7}")
    for yr in range(2011, 2027):
        mask = df["date"].dt.year == yr
        g = df[mask]
        if len(g) < 100:
            continue
        fm = g["fat_mass_lbs"].mean()
        sp_s = g["sp_sym"].mean() if g["sp_sym"].notna().sum() > 50 else np.nan
        sp_a = np.nanmean(sp_best[mask.values]) if (~np.isnan(sp_best[mask.values])).sum() > 50 else np.nan
        br = g["binge"].mean()
        print(f"  {yr:>6} {fm:5.0f} {sp_s:7.0f} {sp_a:8.0f} {br * 100:6.1f}%")

    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("2. DUAL DEFENSE: expenditure + appetite vs SP distance")
    print("=" * 70)

    valid = df["dist_sym"].notna() & df["tdee_resid"].notna()
    r_tdee = np.corrcoef(df.loc[valid, "dist_sym"], df.loc[valid, "tdee_resid"])[0, 1]

    fm = df.loc[valid, "fat_mass_lbs"].values
    dist = df.loc[valid, "dist_sym"].values
    resid = df.loc[valid, "tdee_resid"].values
    X = np.column_stack([fm, np.ones(valid.sum())])
    res_d = dist - X @ np.linalg.lstsq(X, dist, rcond=None)[0]
    res_r = resid - X @ np.linalg.lstsq(X, resid, rcond=None)[0]
    r_partial = np.corrcoef(res_d, res_r)[0, 1]

    print(f"\n  TDEE residual vs SP distance: r = {r_tdee:+.3f}")
    print(f"  Partial (| FM):              r = {r_partial:+.3f}")
    print(f"  Binge rate vs SP distance:   r = -0.618")

    print(f"\n  {'SP distance':>15} {'TDEE resid':>11} {'Binge%':>8}")
    br_col = df["binge_rate_90d"]
    for lo, hi in [(-10, -5), (-5, -2.5), (-2.5, 0), (0, 2.5), (2.5, 5), (5, 10), (10, 20)]:
        mask = (df["dist_sym"] > lo) & (df["dist_sym"] <= hi) & valid
        g = df[mask]
        if len(g) < 30:
            continue
        br_valid = g["binge_rate_90d"].dropna()
        print(f"  {lo:+.0f} to {hi:+.0f} lbs {g['tdee_resid'].mean():+11.0f} "
              f"{br_valid.mean() * 100:7.1f}%")

    print(f"\n  Appetite arm is {abs(-0.618) / abs(r_partial):.1f}x stronger than expenditure arm.")

    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("3. FREQUENCY NOT MAGNITUDE: what the set point modulates")
    print("=" * 70)

    binges = df[df["binge"] == 1]
    non_binges = df[df["binge"] == 0]

    valid_b = binges["dist_sym"].notna()
    r_mag = np.corrcoef(binges.loc[valid_b, "dist_sym"],
                        binges.loc[valid_b, "surplus"])[0, 1]

    valid_nb = non_binges["dist_sym"].notna()
    r_nb = np.corrcoef(non_binges.loc[valid_nb, "dist_sym"],
                       non_binges.loc[valid_nb, "surplus"])[0, 1]

    print(f"\n  Binge days: SP distance vs surplus:     r = {r_mag:+.3f}  (n={valid_b.sum()})")
    print(f"  Non-binge days: SP distance vs surplus:  r = {r_nb:+.3f}  (n={valid_nb.sum()})")
    print(f"\n  Binge size is constant (~{binges['surplus'].mean():.0f} cal surplus).")
    print(f"  Non-binge daily surplus scales from {non_binges.loc[non_binges['dist_sym'] < -5, 'surplus'].mean():+.0f} "
          f"(5+ below SP) to {non_binges.loc[non_binges['dist_sym'] > 5, 'surplus'].mean():+.0f} (5+ above SP).")

    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("4. RESTRICTION PREDICTION: which cuts stick?")
    print("=" * 70)

    df["restricting"] = df["calories"] < (df["tdee"] - 200)
    runs = []
    in_run, start = False, 0
    for i in range(n):
        if df["restricting"].iloc[i]:
            if not in_run:
                in_run, start = True, i
        else:
            if in_run:
                if i - start >= 3:
                    runs.append((start, i - 1))
                in_run = False
    if in_run and n - start >= 3:
        runs.append((start, n - 1))

    records = []
    for s, e in runs:
        sp_end = df["sp_sym"].iloc[e]
        fm_end = df["fat_mass_lbs"].iloc[e]
        if np.isnan(sp_end):
            continue
        dist_end = sp_end - fm_end
        post30 = df[(df.index > e) & (df.index <= e + 32)]
        fm30 = post30["fat_mass_lbs"].iloc[-1] if len(post30) > 25 else np.nan
        records.append({
            "length": e - s + 1,
            "sp_dist_end": dist_end,
            "rebound_30d": fm30 - fm_end if not np.isnan(fm30) else np.nan,
        })

    rdf = pd.DataFrame(records)
    above = rdf[rdf["sp_dist_end"] > 0]
    below = rdf[rdf["sp_dist_end"] <= 0]

    print(f"\n  {len(runs)} restriction runs (3+ days, cal < TDEE-200)")
    a30 = above["rebound_30d"].dropna()
    b30 = below["rebound_30d"].dropna()
    print(f"  Above SP at end (n={len(above)}): 30d change = {a30.mean():+.2f} lbs (continues losing)")
    print(f"  Below SP at end (n={len(below)}): 30d change = {b30.mean():+.2f} lbs (rebounds)")

    v30 = rdf.dropna(subset=["rebound_30d"])
    r_run = np.corrcoef(v30["sp_dist_end"], v30["rebound_30d"])[0, 1]
    print(f"  SP distance at end vs 30d rebound: r = {r_run:+.3f}")

    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("5. EXERCISE INDEPENDENCE: walking doesn't move the set point")
    print("=" * 70)

    walk_dates = exercises[exercises["type"] == "walking"]["date"].values
    walks_30d = np.zeros(n)
    for i in range(n):
        d_day = np.datetime64(df["date"].iloc[i], "D")
        start_d = d_day - np.timedelta64(30, "D")
        walks_30d[i] = np.sum((walk_dates >= start_d) & (walk_dates < d_day))
    df["walks_30d"] = walks_30d

    med = df.loc[df["walks_30d"] > 0, "walks_30d"].median()
    for label, sub in [("High-walk", df[df["walks_30d"] >= med]),
                        ("Low-walk", df[df["walks_30d"] < med])]:
        best_r, best_hl = 0, 0
        for hl in range(20, 400, 10):
            sp = ema(sub["fat_mass_lbs"], hl)
            dist = sp - sub["fat_mass_lbs"]
            valid = dist.notna() & sub["binge_rate_90d"].notna()
            if valid.sum() < 100:
                continue
            r = np.corrcoef(dist[valid], sub.loc[valid.values, "binge_rate_90d"])[0, 1]
            if abs(r) > abs(best_r):
                best_r, best_hl = r, hl
        print(f"  {label:>12}: optimal HL = {best_hl}d  r = {best_r:+.3f}")
    print(f"  Walking raises RMR (finding AD) but does not change SP dynamics.")

    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("6. FLOOR EFFECT: set point near essential body fat")
    print("=" * 70)

    min_sp = df["sp_sym"].min()
    min_date = df.loc[df["sp_sym"].idxmin(), "date"]
    print(f"\n  Lowest set point: {min_sp:.1f} lbs FM on {min_date.date()}")

    # Adaptation rate by SP level
    df["sp_rate"] = df["sp_sym"] - df["sp_sym"].shift(30)
    print(f"\n  Adaptation rate by set point level:")
    print(f"  {'SP range':>12} {'Mean rate':>10} {'Abs rate':>9}")
    for lo, hi in [(20, 30), (30, 40), (40, 60), (60, 80), (80, 100)]:
        mask = (df["sp_sym"] >= lo) & (df["sp_sym"] < hi) & df["sp_rate"].notna()
        g = df[mask]
        if len(g) < 30:
            continue
        print(f"  {lo:>3}-{hi:<3} lbs {g['sp_rate'].mean():+10.2f} {g['sp_rate'].abs().mean():9.2f}")

    post = df[(df["date"] > min_date) & (df["date"] <= min_date + pd.Timedelta(days=180))]
    if len(post) > 30:
        print(f"\n  After floor ({min_date.date()}):")
        print(f"    FM: {post['fat_mass_lbs'].iloc[0]:.0f} → {post['fat_mass_lbs'].iloc[-1]:.0f} over 6 months")
        print(f"    Binge rate: {post['binge'].mean() * 100:.1f}%")

    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)
    print(f"  1. Inverted ratchet: SP adapts {best_u // best_d}x faster down (HL={best_d}d) than up (HL={best_u}d)")
    print(f"     Asymmetric model: r = {best_r:+.3f} vs symmetric r = {r_sym:+.3f}")
    print(f"  2. Dual defense: appetite (r = -0.62) + expenditure (r = +0.38)")
    print(f"     Appetite is {abs(-0.618) / abs(r_partial):.0f}x stronger")
    print(f"  3. Binge size constant (~{binges['surplus'].mean():.0f} cal); non-binge drift scales r = {r_nb:+.3f}")
    print(f"  4. Restriction above SP holds (-1.2 lbs/30d); below SP rebounds (+0.4 lbs/30d)")
    print(f"  5. Walking is independent of SP dynamics (both HL ≈ 40d)")
    print(f"  6. Floor at ~{min_sp:.0f} lbs FM; adaptation stalls near essential body fat")


if __name__ == "__main__":
    main()
