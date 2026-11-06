#!/usr/bin/env python3
"""DD. Forward-looking set point test: does distance below the trailing EMA
predict FUTURE eating pressure, or only reflect PAST eating?

Why this exists
---------------
AG correlates a TRAILING 90-day binge rate with today's distance between a
trailing EMA of fat mass and today's fat mass. Binges in the past 90 days raise
FM above its EMA, so that correlation is built in by the energy balance: it is
reverse causation, not a lipostat reading. CN's Test 1 has the same flaw at a
7-day scale (the 7-day intake being predicted is the same 7 days that move FM
relative to the EMA). AU is the only strictly forward test and its sign is
anti-lipostat (below the EMA → LOWER future intake).

Sign convention (same as AG/AU/CN): sp_dist = EMA(FM) − FM.
  sp_dist > 0  ⇔  FM below its trailing EMA  ⇔  weight has been falling.
A lipostat predicts: sp_dist > 0 → MORE future intake / binges (positive r).
Intake persistence predicts: sp_dist > 0 (you have been restricting) → you keep
restricting for a while (negative r at short leads, decaying to 0).

Tests
-----
1. Binge rate by sp_dist bin: trailing 90 d (AG's version) vs the NEXT 90 d.
2. Lead-lag cross-correlation of sp_dist(t) with intake(t+h) and binge(t+h)
   for h from −90 to +365 days.
3. Parametric AR null: fit an AR(60) to daily surplus, simulate intake with no
   feedback from weight at all, integrate to FM, recompute the same CCF. If the
   observed forward correlations sit inside the null band, the data contain no
   detectable feedback from fat mass to future intake beyond intake persistence.
4. Out-of-sample increment: does sp_dist improve a year-ahead forecast of the
   next 30/90-day mean intake over a forecast from past intake alone?

Data: 2014-01-01 onward, pre-tirzepatide, Kalman FM (P4), consumed intake (P4 intake_true).
"""

from pathlib import Path
import numpy as np
import pandas as pd
from numpy.linalg import lstsq

ROOT = Path(__file__).resolve().parent.parent
import sys
sys.path.insert(0, str(ROOT / "analysis"))
from intake_model import INTAKE_MULTIPLIER, describe as describe_intake_model  # single source of truth
from results import save as save_result, load as load_result
SP_HL = 50
BINGE_THRESHOLD = 1000       # cal above TDEE (AG's definition)
WINDOW = 90
AR_ORDER = 60
N_SIMS = 200
LEADS = [-90, -30, -7, 0, 1, 7, 14, 30, 60, 90, 120, 180, 270, 365]
FORBES_FAT_FRAC = 0.85


def ema(x, hl):
    a = 1 - np.exp(-np.log(2) / hl)
    s = np.empty(len(x))
    s[0] = x[0]
    for j in range(1, len(x)):
        s[j] = s[j - 1] + a * (x[j] - s[j - 1])
    return s


def lagged_corr(x, y, h):
    """corr(x(t), y(t+h)); negative h looks backward."""
    if h >= 0:
        a, b = x[:len(x) - h] if h else x, y[h:]
    else:
        a, b = x[-h:], y[:len(y) + h]
    return float(np.corrcoef(a, b)[0, 1])


def load():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    df = intake[["date", "calories"]].merge(
        kalman[["date", "fat_mass_lbs", "tdee"]], on="date").merge(
        tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df = df[(df["date"] >= "2014-01-01") & (df["effective_level"] == 0)]
    df = df.sort_values("date").reset_index(drop=True)
    # consumed calories from P4 (the series the filter actually used)
    df = df.drop(columns=["calories"]).merge(kalman[["date", "intake_true"]], on="date").rename(columns={"intake_true": "calories"})
    return df


def main():
    df = load()
    fm = df["fat_mass_lbs"].values
    cal = df["calories"].values
    tdee = df["tdee"].values
    binge = (cal > tdee + BINGE_THRESHOLD).astype(float)   # both sides in consumed calories
    sp = ema(fm, SP_HL)
    dist = sp - fm
    n = len(df)
    print(f"DD. Forward set point test — {n} pre-tirz days from 2014, "
          f"{int(binge.sum())} binges ({binge.mean() * 100:.1f}%), SP = EMA(FM, HL={SP_HL}d)")

    # ── 1. Binge rate by distance bin: trailing (AG) vs forward ──────────
    trailing = pd.Series(binge).rolling(WINDOW).mean().values
    forward = pd.Series(binge).rolling(WINDOW).mean().shift(-WINDOW).values
    bins = [-25, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 25]
    tab = pd.DataFrame({"dist": dist, "trailing": trailing, "forward": forward})
    tab["bin"] = pd.cut(tab["dist"], bins)
    g = tab.groupby("bin", observed=True).agg(n=("dist", "size"),
                                              trailing_90d=("trailing", "mean"),
                                              forward_90d=("forward", "mean"))
    g = g[g["n"] >= 30]
    print(f"\n=== 1. Binge rate by sp_dist bin (positive = FM below EMA) ===")
    print(f"  {'bin (lb)':>14} {'n':>5} {'trailing 90d':>13} {'next 90d':>9}")
    for b, row in g.iterrows():
        print(f"  {str(b):>14} {int(row['n']):>5} {row['trailing_90d'] * 100:>12.1f}% "
              f"{row['forward_90d'] * 100:>8.1f}%")
    v = ~np.isnan(trailing)
    r_trail = np.corrcoef(dist[v], trailing[v])[0, 1]
    v = ~np.isnan(forward)
    r_fwd = np.corrcoef(dist[v], forward[v])[0, 1]
    print(f"\n  r(sp_dist, trailing 90d binge rate) = {r_trail:+.3f}   (AG's construction)")
    print(f"  r(sp_dist, NEXT 90d binge rate)     = {r_fwd:+.3f}")
    print(f"  Lipostat predicts r > 0 for the forward version.")

    # ── 2. Lead-lag CCF ─────────────────────────────────────────────────
    print(f"\n=== 2. corr(sp_dist(t), outcome(t+h)) ===")
    print(f"  {'h (days)':>9} {'intake':>8} {'binge':>8}")
    obs_int = {}
    for h in LEADS:
        obs_int[h] = lagged_corr(dist, cal, h)
        print(f"  {h:>9} {obs_int[h]:>+8.3f} {lagged_corr(dist, binge, h):>+8.3f}")

    # ── 3. AR null ──────────────────────────────────────────────────────
    sur = cal - tdee
    p = AR_ORDER
    X = np.column_stack([sur[p - l - 1:n - l - 1] for l in range(p)] + [np.ones(n - p)])
    y = sur[p:]
    beta = lstsq(X, y, rcond=None)[0]
    resid = y - X @ beta
    r2 = 1 - resid.var() / y.var()
    rng = np.random.default_rng(0)
    fwd_leads = [h for h in LEADS if h >= 0]
    sims = np.empty((N_SIMS, len(fwd_leads)))
    for s in range(N_SIMS):
        z = np.zeros(n)
        z[:p] = sur[:p]
        e = rng.choice(resid, n)
        for j in range(p, n):
            z[j] = beta[:-1] @ z[j - p:j][::-1] + beta[-1] + e[j]
        fm_sim = fm[0] + np.cumsum(z * FORBES_FAT_FRAC / 3500.0)
        d_sim = ema(fm_sim, SP_HL) - fm_sim
        sims[s] = [lagged_corr(d_sim, z + tdee, h) for h in fwd_leads]
    print(f"\n=== 3. AR({p}) null (intake persistence only, no feedback from FM; "
          f"R²={r2:.3f}, Σcoef={beta[:-1].sum():.2f}, {N_SIMS} sims) ===")
    print(f"  {'h (days)':>9} {'observed':>9} {'null mean':>10} {'null 2.5%':>10} {'null 97.5%':>11} {'outside?':>9}")
    n_outside = 0
    for k, h in enumerate(fwd_leads):
        lo, hi = np.percentile(sims[:, k], [2.5, 97.5])
        out = obs_int[h] < lo or obs_int[h] > hi
        n_outside += out
        print(f"  {h:>9} {obs_int[h]:>+9.3f} {sims[:, k].mean():>+10.3f} {lo:>+10.3f} {hi:>+11.3f} {'YES' if out else 'no':>9}")
    print(f"  Leads outside the null band: {n_outside}/{len(fwd_leads)}")

    # ── 4. Out-of-sample increment over past intake ─────────────────────
    print(f"\n=== 4. Year-ahead forecast of next-H-day mean intake (train on years < test year) ===")
    df["year"] = df["date"].dt.year
    fc = {}
    for H in (30, 90):
        fut = pd.Series(cal).rolling(H).mean().shift(-H).values
        past = pd.Series(cal).rolling(H).mean().values
        past30 = pd.Series(cal).rolling(30).mean().values
        feats = {"past intake only": [past, past30],
                 "past intake + sp_dist": [past, past30, dist],
                 "past intake + sp_dist + FM": [past, past30, dist, fm],
                 "sp_dist only": [dist]}
        ok = ~np.isnan(fut) & ~np.isnan(past)
        print(f"\n  Horizon {H} d:")
        for label, cols in feats.items():
            preds, ys, sp_coefs = [], [], []
            for yr in sorted(df["year"].unique()):
                tr = ok & (df["year"].values < yr)
                te = ok & (df["year"].values == yr)
                if tr.sum() < 365 or te.sum() < 60:
                    continue
                Xtr = np.column_stack([c[tr] for c in cols] + [np.ones(tr.sum())])
                Xte = np.column_stack([c[te] for c in cols] + [np.ones(te.sum())])
                b = lstsq(Xtr, fut[tr], rcond=None)[0]
                if "sp_dist" in label:
                    sp_coefs.append(b[[k for k, c in enumerate(cols) if c is dist][0]])
                preds.extend(Xte @ b)
                ys.extend(fut[te])
            preds, ys = np.array(preds), np.array(ys)
            rmse = np.sqrt(np.mean((preds - ys) ** 2))
            r = np.corrcoef(preds, ys)[0, 1]
            coef_txt = (f"  sp_dist coef (cal/day per lb below EMA): mean {np.mean(sp_coefs):+.1f}, "
                        f"{sum(c > 0 for c in sp_coefs)}/{len(sp_coefs)} folds positive"
                        if sp_coefs else "")
            print(f"    {label:28s} RMSE={rmse:6.1f}  r={r:+.3f}  n={len(ys)}{coef_txt}")
            fc[(H, label)] = dict(rmse=float(rmse), r=float(r), coef=(float(np.mean(sp_coefs)) if sp_coefs else None),
                                  folds_pos=(int(sum(c > 0 for c in sp_coefs)) if sp_coefs else None), folds=len(sp_coefs))

    r_tr = np.corrcoef(binned["sp_dist"], binned["trailing"])[0, 1] if "binned" in dir() else float("nan")
    b90, p90 = fc[(90, "past intake + sp_dist")], fc[(90, "past intake only")]
    b30, p30 = fc[(30, "past intake + sp_dist")], fc[(30, "past intake only")]
    gain90 = 100 * (p90["rmse"] - b90["rmse"]) / p90["rmse"]
    gain30 = 100 * (p30["rmse"] - b30["rmse"]) / p30["rmse"]
    print("\n=== Verdict ===")
    print("  Trailing-window binge rate vs EMA distance (AG) is reverse causation: past binges")
    print(f"  push FM above the EMA. With forward windows the gradient collapses (r {r_trail:+.2f} → {r_fwd:+.2f})")
    print("  and its residual sign is anti-lipostat; the forward leads of the raw CCF sit inside")
    print(f"  an intake-only AR null ({n_outside}/{len(fwd_leads)} outside). What survives: at a 90-day horizon, EMA distance")
    print(f"  changes a past-intake forecast's RMSE by {gain90:+.1f}% with a coefficient of {b90['coef']:+.0f} cal/day per lb")
    print(f"  below the EMA ({b90['folds_pos']}/{b90['folds']} folds positive; {b90['coef'] / 55:.0%} of the 55 cal/lb claimed);")
    print(f"  at 30 days the RMSE change is {gain30:+.1f}% (coefficient {b30['coef']:+.0f}). The pre-drug decade mostly")
    print("  cannot identify a set point because FM never left its own trailing average; see DF for the drug era.")
    save_result("DD_forward_setpoint_test", r_trailing=r_trail, r_forward=r_fwd, ccf_leads_outside_null=n_outside,
                forecast={f"{h}d_{k}": v for (h, k), v in fc.items()}, rmse_gain_pct_90d=gain90, rmse_gain_pct_30d=gain30)
    print("Wrote analysis/results/DD_forward_setpoint_test.json")


if __name__ == "__main__":
    main()
