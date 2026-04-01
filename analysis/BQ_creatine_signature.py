"""
BQ_creatine_signature.py — Detect creatine water-retention signature in body composition.

Creatine monohydrate causes intracellular water retention (~2-5 lbs at 5g/day).
This water is lean mass, invisible to the Kalman filter which derives FM as
smoothed_weight - known_lean. When creatine starts, FM inflates; when it stops,
FM deflates — unless a composition scan captures the change in lean mass.

Documented creatine history:
  Period 1: 2018-10-29 to 2021-04-28, 2g/day
  Period 2: 2025-05-09 to present, 5g/day

Strategy:
  1. Regression discontinuity: sweep candidate transition dates ±30d,
     fit weight = trend + step at each candidate, report step size and p-value.
  2. Innovation analysis: cumulative Kalman innovations pre vs post transition.
  3. Permutation test: is the cessation innovation drop significant?
  4. Composition scans: BIA lean mass ON vs OFF creatine (period 2).
  5. Confounder controls: restrict period 2 to constant tirz dose.
"""

import pandas as pd
import numpy as np
from scipy import stats
from datetime import timedelta

np.random.seed(42)

# ── Load data ─────────────────────────────────────────────────────
kalman = pd.read_csv("analysis/P4_kalman_daily.csv", parse_dates=["date"])
p1     = pd.read_csv("analysis/P1_smoothed_weight.csv", parse_dates=["date"])
p2     = pd.read_csv("analysis/P2_daily_composition.csv", parse_dates=["date"])
comp   = pd.read_csv("composition/composition.csv", parse_dates=["date"])
intake = pd.read_csv("intake/intake_daily.csv", parse_dates=["date"])
tirz   = pd.read_csv("drugs/tirzepatide.csv", parse_dates=["date"])

df = kalman.merge(p1[["date", "smoothed_weight_lbs"]], on="date", how="left")
df = df.merge(p2[["date", "ffm_lbs"]], on="date", how="left")
df = df.merge(intake[["date", "calories"]], on="date", how="left")
df = df.merge(tirz[["date", "dose_mg", "effective_level"]], on="date", how="left")

transitions = [
    ("Creatine ON  2g (2018-10-29)", "2018-10-29", +1),
    ("Creatine OFF    (2021-04-28)", "2021-04-28", -1),
    ("Creatine ON  5g (2025-05-09)", "2025-05-09", +1),
]

# ── 1. Regression discontinuity sweep ────────────────────────────
print("=" * 72)
print("1. REGRESSION DISCONTINUITY: weight = a + b*days + c*step")
print("   Sweep candidate breakpoint ±30d around documented transitions")
print("=" * 72)

for label, center_str, direction in transitions:
    center = pd.Timestamp(center_str)
    print(f"\n--- {label} ---")
    print(f"{'offset':>7} {'step':>7} {'SE':>6} {'t':>6} {'p':>6}")

    results = []
    for offset in range(-30, 31):
        candidate = center + timedelta(days=offset)
        window = df[
            (df.date >= candidate - timedelta(days=45))
            & (df.date <= candidate + timedelta(days=45))
        ].dropna(subset=["smoothed_weight_lbs"]).copy()
        if len(window) < 20:
            continue

        window["days"] = (window.date - candidate).dt.days
        window["post"] = (window.date >= candidate).astype(int)

        X = np.column_stack([np.ones(len(window)), window.days.values, window.post.values])
        y = window.smoothed_weight_lbs.values

        coefs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ coefs
        mse = np.sum(resid ** 2) / (len(y) - 3)
        try:
            cov = mse * np.linalg.inv(X.T @ X)
        except np.linalg.LinAlgError:
            continue
        se_step = np.sqrt(cov[2, 2])
        t_stat = coefs[2] / se_step
        p_val = 2 * stats.t.sf(abs(t_stat), len(y) - 3)
        results.append(dict(offset=offset, step=coefs[2], se=se_step, t=t_stat, p=p_val, n=len(window)))

    if not results:
        print("  (no valid windows — sparse weight data)")
        continue

    rdf = pd.DataFrame(results)
    for _, r in rdf.iterrows():
        if int(r.offset) % 5 == 0:
            sig = " *" if r.p < 0.05 else ""
            print(f"{int(r.offset):>+7d} {r.step:>+7.2f} {r.se:>6.2f} {r.t:>6.2f} {r.p:>6.3f}{sig}")

    if direction > 0:
        best = rdf.loc[rdf.step.idxmax()]
    else:
        best = rdf.loc[rdf.step.idxmin()]
    print(f"  Peak: offset={int(best.offset):+d}, step={best.step:+.2f} lbs, p={best.p:.4f}")

# ── 2. Period 2 restricted to constant tirz dose ─────────────────
print("\n" + "=" * 72)
print("2. PERIOD 2 ONSET: restricted to constant tirz 7.5mg window")
print("   Tirz dose escalated 7.5→10mg on 2025-05-27 (+18d after creatine)")
print("=" * 72)

onset = pd.Timestamp("2025-05-09")
tirz_change = pd.Timestamp("2025-05-27")

print(f"{'offset':>7} {'step':>7} {'SE':>6} {'t':>6} {'p':>6} {'n':>4}")
for offset in range(-15, 18):
    candidate = onset + timedelta(days=offset)
    window = df[
        (df.date >= candidate - timedelta(days=45))
        & (df.date < tirz_change)
    ].dropna(subset=["smoothed_weight_lbs"]).copy()
    if len(window) < 15:
        continue

    window["days"] = (window.date - candidate).dt.days
    window["post"] = (window.date >= candidate).astype(int)
    X = np.column_stack([np.ones(len(window)), window.days.values, window.post.values])
    y = window.smoothed_weight_lbs.values

    coefs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coefs
    mse = np.sum(resid ** 2) / (len(y) - 3)
    try:
        cov = mse * np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        continue
    se_step = np.sqrt(cov[2, 2])
    t_stat = coefs[2] / se_step
    p_val = 2 * stats.t.sf(abs(t_stat), len(y) - 3)
    sig = " *" if p_val < 0.05 else ""
    if offset % 3 == 0:
        print(f"{offset:>+7d} {coefs[2]:>+7.2f} {se_step:>6.2f} {t_stat:>6.2f} {p_val:>6.3f} {len(window):>4d}{sig}")

# ── 3. Innovation analysis with permutation test ─────────────────
print("\n" + "=" * 72)
print("3. KALMAN INNOVATIONS: cumulative 30d pre vs 30d post")
print("   + permutation test for cessation")
print("=" * 72)

for label, center_str, direction in transitions:
    center = pd.Timestamp(center_str)
    inn = df[["date", "innovation"]].dropna()

    pre  = inn[(inn.date >= center - timedelta(days=30)) & (inn.date < center)]["innovation"]
    post = inn[(inn.date >= center) & (inn.date < center + timedelta(days=30))]["innovation"]

    if len(pre) < 3 or len(post) < 3:
        print(f"\n{label}: too few innovations (pre={len(pre)}, post={len(post)})")
        continue

    delta = post.mean() - pre.mean()
    print(f"\n{label}")
    print(f"  Pre  mean innovation: {pre.mean():+.3f} (n={len(pre)})")
    print(f"  Post mean innovation: {post.mean():+.3f} (n={len(post)})")
    print(f"  Delta: {delta:+.3f}")

    # Permutation test
    all_inn = inn["innovation"].values
    n_perm = 10000
    perm_stats = []
    for _ in range(n_perm):
        idx = np.random.randint(30, len(all_inn) - 30)
        fake_pre = all_inn[idx - len(pre):idx]
        fake_post = all_inn[idx:idx + len(post)]
        if len(fake_pre) == len(pre) and len(fake_post) == len(post):
            perm_stats.append(fake_post.mean() - fake_pre.mean())

    perm_stats = np.array(perm_stats)
    if direction < 0:
        p_val = (perm_stats <= delta).mean()
        tail = "drop"
    else:
        p_val = (perm_stats >= delta).mean()
        tail = "rise"
    print(f"  Permutation p ({tail}): {p_val:.4f}")

# ── 4. Composition scans: lean mass ON vs OFF ────────────────────
print("\n" + "=" * 72)
print("4. COMPOSITION SCANS: BIA lean mass ON vs OFF creatine")
print("=" * 72)

# Period 2 scans (Sep 2024 → Mar 2026)
scans = comp[(comp.date >= "2024-09-01") & (comp.date <= "2026-03-01")].copy()
scans["creatine"] = scans.date.apply(lambda d: "ON 5g" if d >= pd.Timestamp("2025-05-09") else "OFF")
scans["lean_pct"] = scans.lean_mass_lbs / scans.weight_lbs * 100

print("\n  Period 2 scans:")
for cr in ["OFF", "ON 5g"]:
    sub = scans[scans.creatine == cr]
    print(f"    {cr:6s}: lean={sub.lean_mass_lbs.mean():.1f}±{sub.lean_mass_lbs.std():.1f} lbs, "
          f"lean%={sub.lean_pct.mean():.1f}±{sub.lean_pct.std():.1f}%, "
          f"weight={sub.weight_lbs.mean():.1f}, n={len(sub)}")

pre  = scans[scans.creatine == "OFF"]
post = scans[scans.creatine == "ON 5g"]
print(f"\n  Lean delta: {post.lean_mass_lbs.mean() - pre.lean_mass_lbs.mean():+.1f} lbs")
print(f"  Weight delta: {post.weight_lbs.mean() - pre.weight_lbs.mean():+.1f} lbs")
print(f"  Fat delta: {post.fat_mass_lbs.mean() - pre.fat_mass_lbs.mean():+.1f} lbs")
print(f"  Lean% delta: {post.lean_pct.mean() - pre.lean_pct.mean():+.1f} pp")
print("  (Lean% confounded by tirz-driven fat loss)")

# Weight-adjusted lean: regress lean on weight, compare residuals
all_scans = comp[comp.lean_mass_lbs.notna() & comp.weight_lbs.notna()].copy()
all_scans["creatine"] = "OFF"
all_scans.loc[(all_scans.date >= "2018-10-29") & (all_scans.date <= "2021-04-28"), "creatine"] = "ON_2g"
all_scans.loc[all_scans.date >= "2025-05-09", "creatine"] = "ON_5g"

slope, intercept = np.polyfit(all_scans.weight_lbs, all_scans.lean_mass_lbs, 1)
all_scans["lean_residual"] = all_scans.lean_mass_lbs - (slope * all_scans.weight_lbs + intercept)

print(f"\n  Weight-adjusted lean residual (all 70 scans, lean ~ {slope:.3f}*weight + {intercept:.1f}):")
for cr in ["OFF", "ON_2g", "ON_5g"]:
    sub = all_scans[all_scans.creatine == cr]["lean_residual"]
    if len(sub):
        print(f"    {cr:6s}: residual = {sub.mean():+.2f} ± {sub.std():.2f} lbs (n={len(sub)})")

on = all_scans[all_scans.creatine != "OFF"]["lean_residual"]
off = all_scans[all_scans.creatine == "OFF"]["lean_residual"]
if len(on) > 2:
    t, p = stats.ttest_ind(on, off, equal_var=False)
    print(f"    ON vs OFF: delta = {on.mean() - off.mean():+.2f} lbs, t={t:.2f}, p={p:.4f}")

# ── 5. TDEE artifact check ──────────────────────────────────────
print("\n" + "=" * 72)
print("5. TDEE ARTIFACT: does creatine water bias Kalman TDEE?")
print("=" * 72)

df2 = df.merge(p2[["date", "expected_rmr"]], on="date", how="left")
df2["tdee_rmr"] = df2["tdee"] / df2["expected_rmr"]

windows = [
    ("Pre-creatine (2016-2018)",    "2016-01-01", "2018-10-28"),
    ("Creatine 2g (2019-2021)",     "2019-01-01", "2021-04-27"),
    ("Post-creatine (2021-2024)",   "2021-04-29", "2024-12-31"),
    ("Creatine 5g (2025-05+)",      "2025-05-10", "2026-03-01"),
]

for label, start, end in windows:
    mask = (df2.date >= start) & (df2.date <= end) & df2.tdee_rmr.notna()
    sub = df2.loc[mask]
    print(f"  {label}: TDEE/RMR={sub.tdee_rmr.mean():.3f}, TDEE={sub.tdee.mean():.0f}, n={len(sub)}")

# ── 6. Intake confounders ────────────────────────────────────────
print("\n" + "=" * 72)
print("6. INTAKE AROUND TRANSITIONS (confounder check)")
print("=" * 72)

for label, center_str, _ in transitions:
    center = pd.Timestamp(center_str)
    pre  = df[(df.date >= center - timedelta(days=30)) & (df.date < center)]["calories"].dropna()
    post = df[(df.date >= center) & (df.date < center + timedelta(days=30))]["calories"].dropna()
    if len(pre) and len(post):
        print(f"  {label}")
        print(f"    Pre  30d: {pre.mean():.0f} cal/day (n={len(pre)})")
        print(f"    Post 30d: {post.mean():.0f} cal/day (n={len(post)})")
        print(f"    Delta: {post.mean() - pre.mean():+.0f} cal/day")

# ── Summary ──────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)
print("""
Period 1 onset (2g, 2018-10-29): NOT TESTABLE
  Only 1 weight reading within ±60 days. No composition scans nearby.

Period 1 cessation (2021-04-28): MARGINAL
  Regression discontinuity: -1.75 lbs at offset -10d (p=0.057)
  Innovation drop: -1.04 (permutation p=0.11)
  Direction correct (weight drops after cessation) but not significant.
  Expected effect at 2g/day: ~1-2 lbs — near detection threshold.

Period 2 onset (5g, 2025-05-09): CONFOUNDED
  Full window RD: +2.92 lbs at offset +13d (p<0.001)
  BUT restricted to constant tirz 7.5mg: +0.85 lbs (p=0.16)
  Tirz dose escalated 7.5→10mg on 2025-05-27 (+18d), dominating the signal.
  Creatine loading kinetics (2-4 weeks to saturate) overlap the dose change.
  Cannot separate creatine water from tirz dose response.

BIA lean mass: +1.1 lbs ON vs OFF (Period 2 scans)
  But weight-adjusted lean residual shows no significant creatine effect.
  Fat loss from tirz drives the apparent lean gain.

Verdict: The data cannot confirm or reject creatine's expected ~2-5 lb
water-retention signature. Period 1 onset has no data. Period 1 cessation
is marginal but underpowered. Period 2 onset is hopelessly confounded by
a concurrent tirz dose change. A clean detection would require a creatine
washout during a period of stable weight and medication.
""")
