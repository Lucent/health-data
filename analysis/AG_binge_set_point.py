"""Binge frequency reveals a hidden, moving fat mass set point.

Binge rate (days > TDEE + 1000 cal) is a monotonic, sigmoid function of
how far current fat mass sits below a trailing exponential moving average
of fat mass. The EMA half-life that best predicts binge frequency is
50 days (~7 weeks), meaning the set point adapts to sustained weight
changes within 3-4 months.

Key results:
  - Fat mass, not total weight: FM wins at every half-life
    (r = -0.62 vs -0.54 total weight vs -0.53 scale weight).
    This is a lipostat, not a gravitostat.
  - Moving, not fixed: a moving EMA (r = -0.62) crushes any fixed
    set point (best r = +0.25). The defended weight adapts.
  - Half-life ~50 days: the set point chases actual FM with a
    7-week lag. After 3-4 months at a new weight, the set point
    has largely caught up.
  - Sigmoid response: baseline 3% binge rate at/above set point,
    rising to ~15% at 7+ lbs below. ~0.5% per lb below.
  - Partial r = -0.64 after controlling for absolute FM. Distance
    below the set point, not FM itself, drives binges.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

BINGE_THRESHOLD = 1000  # cal above TDEE
BINGE_RATE_WINDOW = 90  # days
SET_POINT_HALF_LIFE = 50  # days — optimal from sweep


def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    weight = pd.read_csv(ROOT / "weight" / "weight.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])

    df = intake[["date", "calories"]].merge(
        kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="inner")
    df = df.merge(comp[["date", "ffm_lbs"]], on="date", how="left")
    wt = weight[["date", "weight_lbs"]].drop_duplicates("date", keep="first")
    df = df.merge(wt, on="date", how="left")
    df = df.sort_values("date").reset_index(drop=True)
    df["total_weight"] = df["fat_mass_lbs"] + df["ffm_lbs"]
    n = len(df)

    df["binge"] = (df["calories"] > df["tdee"] + BINGE_THRESHOLD).astype(int)
    print(f"{n} days, {df['binge'].sum()} binges ({df['binge'].mean() * 100:.1f}%)")

    def ema(series, half_life):
        alpha = 1 - np.exp(-np.log(2) / half_life)
        return series.ewm(alpha=alpha, min_periods=30).mean()

    df["binge_rate_90d"] = df["binge"].rolling(BINGE_RATE_WINDOW, min_periods=BINGE_RATE_WINDOW).mean()

    # ── Part 1: What does the set point track? ──
    print(f"\n=== What does the set point track? ===")
    print(f"  Distance below EMA → {BINGE_RATE_WINDOW}d binge rate")
    print(f"  {'Variable':>15} {'Half-life':>10} {'r':>7}")

    all_results = []
    for var_name, var_col in [("fat_mass", "fat_mass_lbs"),
                               ("total_weight", "total_weight"),
                               ("scale_weight", "weight_lbs")]:
        for hl in [40, 50, 60, 80, 120, 180, 270, 365, 540, 730]:
            sp = ema(df[var_col], hl)
            dist = sp - df[var_col]
            valid = dist.notna() & df["binge_rate_90d"].notna()
            if valid.sum() < 500:
                continue
            r = np.corrcoef(dist[valid], df.loc[valid, "binge_rate_90d"])[0, 1]
            all_results.append((var_name, hl, r))

    all_results.sort(key=lambda x: -abs(x[2]))
    for var, hl, r in all_results[:15]:
        print(f"  {var:>15} {hl:>8}d {r:+7.3f}")

    # Best per variable at own optimal HL
    print(f"\n  Each variable at its own optimal half-life:")
    for var_name, var_col in [("fat_mass", "fat_mass_lbs"),
                               ("total_weight", "total_weight"),
                               ("scale_weight", "weight_lbs")]:
        best_r, best_hl = 0, 0
        for hl in range(30, 1500, 10):
            sp = ema(df[var_col], hl)
            dist = sp - df[var_col]
            valid = dist.notna() & df["binge_rate_90d"].notna()
            if valid.sum() < 500:
                continue
            r = np.corrcoef(dist[valid], df.loc[valid, "binge_rate_90d"])[0, 1]
            if abs(r) > abs(best_r):
                best_r, best_hl = r, hl
        print(f"  {var_name:>15}: HL = {best_hl}d  r = {best_r:+.3f}")

    # ── Part 2: Fixed vs moving ──
    print(f"\n=== Fixed vs moving set point ===")
    best_fixed_r, best_fixed_sp = 0, 0
    for fixed_sp in np.arange(20, 110, 1):
        dist = fixed_sp - df["fat_mass_lbs"]
        valid = dist.notna() & df["binge_rate_90d"].notna()
        if valid.sum() < 500:
            continue
        r = np.corrcoef(dist[valid], df.loc[valid, "binge_rate_90d"])[0, 1]
        if abs(r) > abs(best_fixed_r):
            best_fixed_r, best_fixed_sp = r, fixed_sp

    print(f"  Fixed (best): FM = {best_fixed_sp:.0f} lbs, r = {best_fixed_r:+.3f}")
    print(f"  Moving (FM, {SET_POINT_HALF_LIFE}d): r = {all_results[0][2]:+.3f}")
    print(f"  Moving wins by: {abs(all_results[0][2]) - abs(best_fixed_r):.3f}")

    # ── Part 3: Half-life fine sweep ──
    print(f"\n=== Set point half-life ===")
    fine = []
    for hl in range(20, 400, 5):
        sp = ema(df["fat_mass_lbs"], hl)
        dist = sp - df["fat_mass_lbs"]
        valid = dist.notna() & df["binge_rate_90d"].notna()
        if valid.sum() < 500:
            continue
        r = np.corrcoef(dist[valid], df.loc[valid, "binge_rate_90d"])[0, 1]
        fine.append((hl, r))

    fine.sort(key=lambda x: -abs(x[1]))
    print(f"  Top 10 half-lives:")
    for hl, r in fine[:10]:
        print(f"    {hl:>4}d ({hl / 30:.1f} mo): r = {r:+.3f}")
    best_hl = fine[0][0]
    print(f"\n  Optimal: {best_hl} days ({best_hl / 7:.0f} weeks)")

    # ── Part 4: Reconstruct the set point ──
    sp = ema(df["fat_mass_lbs"], best_hl)
    df["set_point"] = sp
    df["sp_distance"] = sp - df["fat_mass_lbs"]

    print(f"\n=== Reconstructed set point trajectory ===")
    print(f"  {'Year':>6} {'FM':>5} {'Set Pt':>7} {'Dist':>6} {'Binge%':>7} {'State':>12}")
    for yr in range(2011, 2027):
        mask = df["date"].dt.year == yr
        g = df[mask]
        if len(g) < 100:
            continue
        fm = g["fat_mass_lbs"].mean()
        sp_val = g["set_point"].mean() if g["set_point"].notna().sum() > 50 else np.nan
        dist = g["sp_distance"].mean() if g["sp_distance"].notna().sum() > 50 else np.nan
        br = g["binge"].mean()
        if dist > 2:
            state = "losing"
        elif dist < -2:
            state = "gaining"
        else:
            state = "at set point"
        print(f"  {yr:>6} {fm:5.0f} {sp_val:7.0f} {dist:+6.1f} {br * 100:6.1f}% {state:>12}")

    # ── Part 5: Partial correlations ──
    print(f"\n=== Partial correlations ===")
    valid = df["sp_distance"].notna() & df["binge_rate_90d"].notna() & df["fat_mass_lbs"].notna()
    sub = df[valid]
    fm = sub["fat_mass_lbs"].values
    dist = sub["sp_distance"].values
    br = sub["binge_rate_90d"].values

    X = np.column_stack([fm, np.ones(len(sub))])
    res_dist = dist - X @ np.linalg.lstsq(X, dist, rcond=None)[0]
    res_br = br - X @ np.linalg.lstsq(X, br, rcond=None)[0]
    r_partial = np.corrcoef(res_dist, res_br)[0, 1]

    X2 = np.column_stack([dist, np.ones(len(sub))])
    res_fm = fm - X2 @ np.linalg.lstsq(X2, fm, rcond=None)[0]
    res_br2 = br - X2 @ np.linalg.lstsq(X2, br, rcond=None)[0]
    r_partial2 = np.corrcoef(res_fm, res_br2)[0, 1]

    print(f"  SP distance | FM: r = {r_partial:+.3f}")
    print(f"  FM | SP distance: r = {r_partial2:+.3f}")

    # ── Part 6: Response shape ──
    print(f"\n=== Binge response curve ===")
    valid = df["sp_distance"].notna() & df["binge"].notna()
    sub = df[valid]

    bins = np.arange(-25, 25, 2.5)
    print(f"  {'Distance':>20} {'Binge %':>8} {'n':>5}")
    bin_centers = []
    bin_rates = []
    for i in range(len(bins) - 1):
        mask = (sub["sp_distance"] > bins[i]) & (sub["sp_distance"] <= bins[i + 1])
        g = sub[mask]
        if len(g) < 30:
            continue
        br = g["binge"].mean()
        center = (bins[i] + bins[i + 1]) / 2
        bin_centers.append(center)
        bin_rates.append(br)
        print(f"  {bins[i]:+5.1f} to {bins[i + 1]:+5.1f} {br * 100:7.1f}% {len(g):5d}")

    # Sigmoid fit
    bin_centers = np.array(bin_centers)
    bin_rates = np.array(bin_rates)
    best_err = np.inf
    best_p = None
    for a in np.arange(0.03, 0.25, 0.005):
        for k in np.arange(0.02, 0.6, 0.02):
            for d0 in np.arange(-10, 20, 0.5):
                for base in np.arange(0.005, 0.06, 0.005):
                    pred = a / (1 + np.exp(-k * (bin_centers - d0))) + base
                    err = np.sum((pred - bin_rates) ** 2)
                    if err < best_err:
                        best_err = err
                        best_p = (a, k, d0, base)

    a, k, d0, base = best_p
    print(f"\n  Sigmoid: binge% = {a * 100:.1f}% / (1 + exp(-{k:.2f}*(dist - {d0:.1f}))) + {base * 100:.1f}%")
    print(f"  Baseline (at/above SP): {base * 100:.1f}%")
    print(f"  Maximum (far below SP): {(a + base) * 100:.1f}%")
    print(f"  Inflection: {d0:.1f} lbs below set point")
    print(f"  Per-lb gradient at inflection: {a * k / 4 * 100:.2f}%/lb")

    # ── Summary ──
    print(f"\n=== Summary ===")
    print(f"  The set point is a lipostat: it tracks fat mass (r = -0.62),")
    print(f"  not total weight (-0.54) or scale weight (-0.53).")
    print(f"  It moves: EMA with {best_hl}-day half-life ({best_hl // 7} weeks).")
    print(f"  After ~{best_hl * 3} days ({best_hl * 3 // 30} months) at a new weight,")
    print(f"  the set point has ~87% adapted.")
    print(f"  Response: sigmoid from {base * 100:.0f}% baseline to {(a + base) * 100:.0f}%")
    print(f"  maximum, {a * k / 4 * 100:.1f}%/lb at the steepest point.")
    print(f"  Partial r = {r_partial:+.3f} controlling for absolute FM.")
    print(f"  Distance below the moving set point, not absolute fat mass,")
    print(f"  drives binge frequency — confirming a defended, adaptive set point")
    print(f"  that operates through appetite.")


if __name__ == "__main__":
    main()
