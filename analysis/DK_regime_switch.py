#!/usr/bin/env python3
"""DK. Behavioral regime switch: date the end of the willpower era objectively, then ask what predicted it.

DH showed no set point model puts its peak gap at the 2014Q3-Q4 binge eruption. The alternative is a
behavioral regime switch: the restriction regime ended. Here intake is modeled as a hidden regime process
and the switch is dated without reference to weight, binges or the subject's own epoch annotations.

Methods (all on log consumed intake = log(P4 intake_true), fasting days floored at 200 kcal):
  split   single best changepoint in the mean (first binary-segmentation split) across windows/transforms
  HMM     Gaussian hidden Markov model, K = 2 and 3, EM (Baum-Welch) + Viterbi, on daily values and on
          non-overlapping 7- and 14-day block means. Era-scale fits constrain the expected dwell to
          ≥ DWELL_ERA days (sticky M-step); the episode-scale fit uses ≥ DWELL_EP days.
          Restriction regime = lowest-mean state.
  PELT    exact penalized changepoint detection on the mean (L2 cost, min segment 30 d, penalty β·σ²·log n).
Switch dates are compared across methods/windows and against intake/diet_epochs.csv.

Switch hazard: on each restriction-regime day, y = 1 if the regime ends within 30 days. Discrete-time
logistic hazard on standardized predictors (univariate, ridge 1e-3; multivariate, ridge 1), year-block
bootstrap CIs. Predictors: FM, LM deficit (running max LM − LM), days in regime, cumulative deficit,
trailing-30 d intake SD, trailing-30 d binge count, Kalman TDEE / expected RMR, season (cos, sin of day
of year), frozen-SP gap (FM₀ − FM) and EMA-100 d SP gap. A "lead" table says how long before the era exit
each predictor had already reached 90 % of its exit value: a pressure that plateaued a year early cannot
date the switch. Secondary: same episode model on 2017-2024 (keto, weekend fasting, potato diets) and a
pooled hazard over every completed pre-drug restriction episode.

Command:  python analysis/DK_regime_switch.py
Artifact: analysis/DK_regime_daily.csv (date, era/episode regime labels per HMM variant, true kcal, binge, FM, LM, TDEE)
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
import sys
sys.path.insert(0, str(ROOT / "analysis"))
from intake_model import INTAKE_MULTIPLIER, describe as describe_intake_model  # single source of truth
from results import save as save_result, load as load_result
MULT = INTAKE_MULTIPLIER
HORIZON = 30          # hazard horizon (days)
DWELL_ERA = 120       # minimum expected dwell (days) for the era-scale HMM
DWELL_EP = 28         # minimum expected dwell (days) for the episode-scale HMM
WIN1 = ("2011-04-20", "2017-12-31")
WIN2 = ("2017-01-01", "2024-09-16")   # pre-drug secondary window
RNG = np.random.default_rng(0)


# ----------------------------------------------------------------------------------------- models
def hmm_fit(x, K, min_dwell, block, iters=300, seed=0):
    """1-D Gaussian HMM by EM with scaled forward-backward; diag(A) floored at 1 − block/min_dwell."""
    n = len(x); rng = np.random.default_rng(seed); floor = 1 - block / min_dwell
    mu = np.quantile(x, (np.arange(K) + 0.5) / K) + rng.normal(0, 0.02, K) * (seed > 0)
    var = np.full(K, x.var()); A = np.full((K, K), 0.02 / (K - 1)); np.fill_diagonal(A, 0.98); pi = np.full(K, 1 / K)
    ll_old = -np.inf
    for _ in range(iters):
        B = np.exp(-0.5 * (x[:, None] - mu) ** 2 / var) / np.sqrt(2 * np.pi * var)
        al = np.zeros((n, K)); c = np.zeros(n)
        al[0] = pi * B[0]; c[0] = al[0].sum(); al[0] /= c[0]
        for t in range(1, n):
            al[t] = (al[t - 1] @ A) * B[t]; c[t] = al[t].sum(); al[t] /= c[t]
        be = np.ones((n, K))
        for t in range(n - 2, -1, -1):
            be[t] = (A @ (B[t + 1] * be[t + 1])) / c[t + 1]
        g = al * be; g /= g.sum(1, keepdims=True)
        xi = (al[:-1].T @ ((B[1:] * be[1:]) / c[1:, None])) * A
        pi = g[0]; A = xi / xi.sum(1, keepdims=True)
        for i in range(K):                      # sticky projection
            if A[i, i] < floor:
                off = A[i].sum() - A[i, i]; A[i] = A[i] * (1 - floor) / max(off, 1e-12); A[i, i] = floor
        mu = (g * x[:, None]).sum(0) / g.sum(0)
        var = np.maximum((g * (x[:, None] - mu) ** 2).sum(0) / g.sum(0), 1e-4)
        ll = np.log(c).sum()
        if ll - ll_old < 1e-6:
            break
        ll_old = ll
    lB = -0.5 * (x[:, None] - mu) ** 2 / var - 0.5 * np.log(2 * np.pi * var); lA = np.log(A)
    V = np.zeros((n, K)); bp = np.zeros((n, K), int); V[0] = np.log(pi + 1e-300) + lB[0]
    for t in range(1, n):
        s = V[t - 1][:, None] + lA; bp[t] = s.argmax(0); V[t] = s.max(0) + lB[t]
    path = np.zeros(n, int); path[-1] = V[-1].argmax()
    for t in range(n - 2, -1, -1):
        path[t] = bp[t + 1, path[t + 1]]
    order = np.argsort(mu); relabel = np.empty(K, int); relabel[order] = np.arange(K)
    return dict(mu=mu[order], sd=np.sqrt(var[order]), A=A[np.ix_(order, order)], ll=ll, path=relabel[path])


def hmm_best(x, K, min_dwell, block, restarts=4):
    return max((hmm_fit(x, K, min_dwell, block, seed=s) for s in range(restarts)), key=lambda f: f["ll"])


def pelt(x, beta, min_len=30):
    """Exact PELT, L2 cost on the mean. Returns changepoint indices (start of each new segment)."""
    n = len(x); S = np.concatenate([[0], np.cumsum(x)]); S2 = np.concatenate([[0], np.cumsum(x * x)])
    F = np.full(n + 1, np.inf); F[0] = -beta; last = np.zeros(n + 1, int); cand = [0]
    for t in range(min_len, n + 1):
        c = np.array([k for k in cand if t - k >= min_len])
        if len(c) == 0:
            continue
        L = t - c; cost = (S2[t] - S2[c]) - (S[t] - S[c]) ** 2 / L
        tot = F[c] + cost + beta; j = tot.argmin(); F[t] = tot[j]; last[t] = c[j]
        keep = [k for k, fk, ck in zip(c, F[c], cost) if fk + ck <= F[t]]
        cand = keep + [k for k in cand if t - k < min_len] + [t]
    cps = []; t = n
    while t > 0:
        t = last[t]
        if t > 0:
            cps.append(t)
    return sorted(cps)


def best_split(x, min_len=90):
    """Index of the single changepoint that most reduces SSE (binary segmentation, first split)."""
    n = len(x); S = np.cumsum(x); S2 = np.cumsum(x * x); k = np.arange(min_len, n - min_len)
    sse = (S2[k - 1] - S[k - 1] ** 2 / k) + ((S2[-1] - S2[k - 1]) - (S[-1] - S[k - 1]) ** 2 / (n - k))
    return int(k[sse.argmin()])


def run_bounds(p):
    cuts = np.flatnonzero(np.diff(p)) + 1
    starts = np.concatenate([[0], cuts]); ends = np.concatenate([cuts, [len(p)]])
    return [(s, e, int(p[s])) for s, e in zip(starts, ends)]


def logistic(X, y, lam=1e-3, iters=50):
    X1 = np.column_stack([np.ones(len(y)), X]); w = np.zeros(X1.shape[1]); R = np.eye(X1.shape[1]) * lam; R[0, 0] = 0
    for _ in range(iters):
        p = 1 / (1 + np.exp(-X1 @ w)); W = p * (1 - p)
        step = np.linalg.solve((X1 * W[:, None]).T @ X1 + R, X1.T @ (y - p) - R @ w); w += step
        if np.abs(step).max() < 1e-8:
            break
    return w


def ema(x, hl):
    a = 1 - np.exp(-np.log(2) / hl); s = np.empty(len(x)); s[0] = x[0]
    for j in range(1, len(x)):
        s[j] = s[j - 1] + a * (x[j] - s[j - 1])
    return s


# ------------------------------------------------------------------------------------------- data
def load():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kal = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    ep = pd.read_csv(ROOT / "intake" / "diet_epochs.csv", parse_dates=["start", "end"])
    d = intake[["date", "calories"]].merge(kal[["date", "fat_mass_lbs", "lm_lbs", "tdee", "intake_true"]], on="date")
    d = d.merge(comp[["date", "expected_rmr"]], on="date", how="left").dropna(subset=["fat_mass_lbs"])
    d = d.sort_values("date").reset_index(drop=True)
    d["cal"] = d["intake_true"]      # consumed calories as used by P4
    d["logcal"] = np.log(np.maximum(d["cal"], 200))
    d["binge"] = (d["calories"] > d["tdee"] + 1000).astype(float)
    d["binge_abs"] = (d["cal"] > 3500).astype(float)
    d["balance"] = d["cal"] - d["tdee"]
    return d, ep


def window(d, lo, hi):
    return d[(d["date"] >= lo) & (d["date"] <= hi)].reset_index(drop=True)


def fit_hmms(w, min_dwell, tag):
    x = w["logcal"].values; n = len(x); out = pd.DataFrame({"date": w["date"]})
    print(f"\n--- HMM fits, {tag} (expected dwell ≥ {min_dwell} d), {n} days ---")
    print(f"  {'model':<14} {'K':>2} {'state means (true kcal/day)':<30} {'sd(log)':<22} {'dwell (d)':<20} loglik")
    for b, bname in ((1, "daily"), (7, "wk7"), (14, "wk14")):
        nb = n // b; xb = x[:nb * b].reshape(nb, b).mean(1)
        for K in (2, 3):
            f = hmm_best(xb, K, min_dwell, b)
            lab = np.repeat(f["path"], b); out[f"hmm_{bname}_K{K}"] = np.concatenate([lab, np.full(n - len(lab), lab[-1])])
            print(f"  {'hmm_' + bname:<14} {K:>2} {str(np.round(np.exp(f['mu'])).astype(int)):<30} "
                  f"{str(np.round(f['sd'], 3)):<22} {str(np.round(b / (1 - np.diag(f['A']))).astype(int)):<20} {f['ll']:.1f}")
    return out


def predictors(w, s, e):
    fm = w["fat_mass_lbs"].values; lm = w["lm_lbs"].values; idx = np.arange(s, e)
    doy = 2 * np.pi * w["date"].dt.day_of_year.values[idx] / 365.25
    return pd.DataFrame({
        "fm": fm[idx], "lm_deficit": (np.maximum.accumulate(lm) - lm)[idx], "days_in_regime": idx - s,
        "cum_deficit_Mcal": -np.cumsum(w["balance"].values[s:e]) / 1000,
        "intake_sd30": w["cal"].rolling(30, min_periods=10).std().values[idx],
        "binge30": w["binge"].rolling(30, min_periods=1).sum().values[idx],
        "tdee_ratio": (w["tdee"] / w["expected_rmr"]).values[idx],
        "season_cos": np.cos(doy), "season_sin": np.sin(doy),
        "gap_frozen": fm[0] - fm[idx], "gap_ema100": (ema(fm, 100) - fm)[idx]}, index=idx)


def episodes(w, lab):
    """Completed restriction episodes (state 0) with predictors and hazard outcome; censored last run dropped."""
    Ps, ys = [], []
    for s, e, l in run_bounds((lab > 0).astype(int)):
        if l == 0 and e < len(w):
            pp = predictors(w, s, e); Ps.append(pp); ys.append(((e - 1 - pp.index.values) < HORIZON).astype(float))
    return Ps, ys


def hazard_table(P, y, years, label, boots=300):
    names = list(P.columns); Z = ((P - P.mean()) / P.std(ddof=0)).fillna(0).values
    mcols = [j for j, k in enumerate(names) if k != "gap_frozen"]          # gap_frozen ≡ −FM within a window
    print(f"\n--- Hazard: {label}: n = {len(y)} restriction days, {int(y.sum())} with exit ≤ {HORIZON} d ---")
    print(f"  {'predictor':<18} {'β (per SD)':>10} {'95% CI (year-block boot)':>26} {'multi β':>9} {'95% CI':>18}")
    uy = np.unique(years); bu = {k: [] for k in names}; bm = {k: [] for k in names}
    w_multi = dict(zip(mcols, logistic(Z[:, mcols], y, lam=1.0)[1:])); w_uni = [logistic(Z[:, [j]], y)[1] for j in range(len(names))]
    for _ in range(boots):
        pick = RNG.choice(uy, len(uy), replace=True); ii = np.concatenate([np.flatnonzero(years == p) for p in pick])
        if y[ii].sum() < 3 or y[ii].sum() > len(ii) - 3:
            continue
        wm = logistic(Z[ii][:, mcols], y[ii], lam=1.0)[1:]
        for jj, j in enumerate(mcols):
            bm[names[j]].append(wm[jj])
        for j, k in enumerate(names):
            bu[k].append(logistic(Z[ii][:, [j]], y[ii])[1])
    rows = []
    for j, k in enumerate(names):
        cu = np.percentile(bu[k], [2.5, 97.5]) if bu[k] else (np.nan, np.nan)
        cm = np.percentile(bm[k], [2.5, 97.5]) if bm[k] else (np.nan, np.nan); wmj = w_multi.get(j, np.nan)
        print(f"  {k:<18} {w_uni[j]:>+10.2f} {f'[{cu[0]:+.2f}, {cu[1]:+.2f}]':>26} {wmj:>+9.2f} {f'[{cm[0]:+.2f}, {cm[1]:+.2f}]':>18}")
        rows.append((label.split(",")[0], k, w_uni[j], cu[0], cu[1], wmj, cm[0], cm[1]))
    print(f"  ({len(bu[names[0]])} usable bootstrap draws; β per SD of predictor; CI excluding 0 = signal; gap_frozen not in multi)")
    if max(abs(b) for b in w_uni) > 10:
        print("  |β| > 10 = perfect separation: the final 30 days are the era's extreme for every monotone predictor, so the")
        print("  one-event hazard identifies calendar order, not which pressure caused the exit.")
    return rows


def regime_table(w, lab, title):
    print(f"\n--- {title} ---")
    print(f"  {'regime':<9} {'start':>10} {'end':>10} {'days':>5} {'kcal/d':>7} {'binge%':>7} {'abs%':>5} {'TDEE':>5} "
          f"{'FM in':>6} {'FM out':>7} {'LM in':>6} {'LM out':>7} {'bal/d':>6}")
    for s, e, l in run_bounds(lab):
        g = w.iloc[s:e]
        print(f"  {('restrict' if l == 0 else 'free' if l == lab.max() else 'mid'):<9} {g['date'].iloc[0].date()!s:>10} "
              f"{g['date'].iloc[-1].date()!s:>10} {e - s:>5} {g['cal'].mean():>7.0f} {g['binge'].mean() * 100:>7.1f} "
              f"{g['binge_abs'].mean() * 100:>5.1f} {g['tdee'].mean():>5.0f} {g['fat_mass_lbs'].iloc[0]:>6.1f} "
              f"{g['fat_mass_lbs'].iloc[-1]:>7.1f} {g['lm_lbs'].iloc[0]:>6.1f} {g['lm_lbs'].iloc[-1]:>7.1f} {g['balance'].mean():>+6.0f}")


def main():
    d, ep = load()
    w1 = window(d, *WIN1); x1 = w1["logcal"].values

    print("=" * 110)
    print("DK. Regime switch: when did the restriction era end, and what predicted it?")
    print("=" * 110)
    print("\n--- Single best changepoint (binary segmentation, first split, min segment 90 d) ---")
    dates = []
    for lo, hi in (("2011-04-20", "2015-12-31"), ("2011-04-20", "2016-12-31"), WIN1, ("2011-04-20", "2018-12-31"), ("2012-01-01", "2017-12-31")):
        ww = window(d, lo, hi)
        for col, nm in (("logcal", "log kcal"), ("cal", "kcal")):
            k = best_split(ww[col].values); dt = ww["date"].iloc[k]; dates.append(dt)
            print(f"  {lo}..{hi} {nm:<8}: {dt.date()}  before {ww['cal'].iloc[:k].mean():.0f} → after {ww['cal'].iloc[k:].mean():.0f} kcal/day")
    sd = pd.Series(dates).sort_values()
    print(f"  median {sd.iloc[len(sd) // 2].date()}, spread {(sd.max() - sd.min()).days} d across {len(sd)} window/transform variants")

    sig2 = np.median(np.abs(np.diff(x1))) ** 2 / (0.6745 ** 2 * 2); n = len(x1)
    print("\n--- PELT on daily log intake, 2011-2017 (first changepoint = restriction era end; others listed) ---")
    for bm in (4, 8, 16):
        cps = pelt(x1, bm * sig2 * np.log(n)); ends = cps[1:] + [n]
        print(f"  β={bm:>2}·σ²·log n: {len(cps):>2} cps: " + ", ".join(f"{w1['date'][c].date()}→{np.exp(x1[c:e].mean()):.0f}" for c, e in zip(cps, ends)))

    era = fit_hmms(w1, DWELL_ERA, "era scale 2011-2017")
    print("\n  Restriction-era end (first exit from the lowest state lasting ≥ 90 d) per era-scale model:")
    for c in [c for c in era.columns if c.startswith("hmm")]:
        rb = run_bounds((era[c].values > 0).astype(int)); ex = [w1["date"].iloc[s] for s, e, l in rb if l == 1 and e - s >= 90]
        print(f"    {c:<14} {ex[0].date() if ex else '-'}   (lowest-state runs: "
              + ", ".join(f"{w1['date'].iloc[s].date()}..{w1['date'].iloc[e - 1].date()}" for s, e, l in rb if l == 0) + ")")
    ref_era = "hmm_wk7_K2"
    regime_table(w1, era[ref_era].values, f"Era-scale regimes ({ref_era}, dwell ≥ {DWELL_ERA} d), 2011-2017")
    epi = fit_hmms(w1, DWELL_EP, "episode scale 2011-2017")
    regime_table(w1, epi[ref_era].values, f"Episode-scale regimes ({ref_era}, dwell ≥ {DWELL_EP} d), 2011-2017")

    print("\n--- Quarterly anatomy of the switch and the eruption ---")
    q = w1.groupby(w1["date"].dt.to_period("Q")).agg(kcal=("cal", "mean"), binge=("binge", "mean"), abs_=("binge_abs", "mean"),
                                                        tdee=("tdee", "mean"), fm=("fat_mass_lbs", "mean"), bal=("balance", "mean"))
    print(f"  {'Q':>7} {'kcal':>6} {'binge%':>7} {'>3500%':>7} {'TDEE':>6} {'FM':>6} {'bal/d':>6}")
    for qq, r in q.loc["2013Q1":"2015Q4"].iterrows():
        print(f"  {str(qq):>7} {r['kcal']:>6.0f} {r['binge'] * 100:>7.1f} {r['abs_'] * 100:>7.1f} {r['tdee']:>6.0f} {r['fm']:>6.1f} {r['bal']:>+6.0f}")
    print("  binge% is TDEE-relative (logged > TDEE + 1000); >3500% is an absolute threshold on true kcal.")

    print("\n--- Regime boundaries vs hand-annotated diet epochs, 2011-2017 ---")
    ep1 = ep[(ep["start"] >= WIN1[0]) & (ep["start"] <= WIN1[1])]
    for nm, lab in (("era", era[ref_era].values), ("episode", epi[ref_era].values)):
        bounds = [w1["date"].iloc[s] for s, e, l in run_bounds(lab)[1:]]
        for _, r in ep1.iterrows():
            near = min(bounds, key=lambda b: abs((b - r["start"]).days))
            print(f"  {nm:<8} epoch {r['label']:<24} starts {r['start'].date()}  nearest boundary {near.date()} ({(near - r['start']).days:+d} d)")
        within = np.mean([min(abs((b - r["start"]).days) for _, r in ep1.iterrows()) <= 30 for b in bounds])
        hit = np.mean([min(abs((b - r["start"]).days) for b in bounds) <= 30 for _, r in ep1.iterrows()])
        print(f"  {nm:<8} {len(bounds)} boundaries: {within * 100:.0f}% within 30 d of an epoch start; {hit * 100:.0f}% of {len(ep1)} epoch starts have a boundary within 30 d")

    # --------------- hazard: era scale (one event) with pressure-lead table, then episode scale pooled
    rows = []
    rb = run_bounds((era[ref_era].values > 0).astype(int))
    e0 = [s for s, e, l in rb if l == 1 and e - s >= 90][0]; s0 = 0; P0 = predictors(w1, s0, e0)
    era_end = w1["date"].iloc[e0 - 1]
    print(f"\n--- Pressure lead in the restriction era {w1['date'].iloc[s0].date()}..{era_end.date()} (start of record to first ≥ 90 d exit) ---")
    print(f"  {'predictor':<18} {'at exit':>9} {'era max':>9} {'exit pctile':>12} {'lead (d)':>9} {'argmax (d before exit)':>23}")
    for k in P0.columns:
        v = P0[k].values; close = np.abs(v - v[-1]) <= 0.1 * (np.nanmax(v) - np.nanmin(v)); first = int(np.argmax(close))
        print(f"  {k:<18} {v[-1]:>9.2f} {np.nanmax(v):>9.2f} {np.nanmean(v < v[-1]) * 100:>11.0f}% {len(v) - 1 - first:>9d} {len(v) - 1 - int(np.nanargmax(v)):>23d}")
    print("  lead = days before exit at which the predictor first came within 10 % (of its era range) of its exit value.")
    y0 = ((e0 - 1 - P0.index.values) < HORIZON).astype(float)
    rows += hazard_table(P0, y0, w1["date"].dt.year.values[P0.index.values],
                         "era-scale 2011-2013, ONE event (the era exit): slow monotone predictors are maximal at exit by construction")
    Ps1, ys1 = episodes(w1, epi[ref_era].values)
    P1 = pd.concat(Ps1); y1 = np.concatenate(ys1); yrs1 = w1["date"].dt.year.values[P1.index.values]
    rows += hazard_table(P1, y1, yrs1, f"episode-scale 2011-2017, {len(Ps1)} completed restriction episodes")

    # --------------- secondary: 2017-2024
    w2 = window(d, *WIN2); epi2 = fit_hmms(w2, DWELL_EP, "episode scale 2017-2024 (pre-drug)")
    lab2 = epi2[ref_era].values
    print(f"\n--- Restriction episodes 2017-2024 ({ref_era}, dwell ≥ {DWELL_EP} d) with annotated epoch covering the exit day ---")
    print(f"  {'start':>10} {'end':>10} {'days':>5} {'kcal/d':>7} {'binge%':>7} {'bal/d':>6} {'FM in':>6} {'FM out':>7} {'LMdef':>6} {'gap100':>7}  epoch at exit")
    ep2 = ep[(ep["end"] >= WIN2[0]) & (ep["start"] <= WIN2[1])]
    for s, e, l in run_bounds((lab2 > 0).astype(int)):
        if l:
            continue
        g = w2.iloc[s:e]; pe = predictors(w2, s, e); dx = g["date"].iloc[-1]
        hit = ep2[(ep2["start"] <= dx) & (ep2["end"] >= dx)]["label"].tolist()
        print(f"  {g['date'].iloc[0].date()!s:>10} {dx.date()!s:>10} {e - s:>5} {g['cal'].mean():>7.0f} {g['binge'].mean() * 100:>7.1f} "
              f"{g['balance'].mean():>+6.0f} {g['fat_mass_lbs'].iloc[0]:>6.1f} {g['fat_mass_lbs'].iloc[-1]:>7.1f} "
              f"{pe['lm_deficit'].iloc[-1]:>6.1f} {pe['gap_ema100'].iloc[-1]:>+7.1f}  {', '.join(hit) or '-'}")
    lab3 = epi2["hmm_wk7_K3"].values
    print("  Annotated restriction epochs: share of days in the K=2 restriction state / in K=3 states 0+1 (below the top state):")
    for _, r in ep2[ep2["label"].str.contains("keto|fasting|potato|restriction|protein")].iterrows():
        m = ((w2["date"] >= r["start"]) & (w2["date"] <= r["end"])).values
        print(f"    {r['label']:<24} {r['start'].date()}..{r['end'].date()} ({m.sum():>3} d, {w2['cal'][m].mean():.0f} kcal): "
              f"K2 restrict {(lab2[m] == 0).mean() * 100:3.0f}%   K3 low+mid {(lab3[m] < 2).mean() * 100:3.0f}%")
    Ps2, ys2 = episodes(w2, lab2)
    P2 = pd.concat(Ps2); y2 = np.concatenate(ys2); yrs2 = w2["date"].dt.year.values[P2.index.values]
    rows += hazard_table(P2, y2, yrs2, f"episode-scale 2017-2024, {len(Ps2)} completed restriction episodes")
    rows += hazard_table(pd.concat([P1, P2]), np.concatenate([y1, y2]), np.concatenate([yrs1, yrs2]),
                         f"pooled episode-scale 2011-2024, {len(Ps1) + len(Ps2)} completed restriction episodes")

    art = pd.concat([era.add_prefix("era_").rename(columns={"era_date": "date"}).merge(epi.add_prefix("ep_").rename(columns={"ep_date": "date"}), on="date"),
                     epi2.add_prefix("ep_").rename(columns={"ep_date": "date"})]).drop_duplicates("date", keep="first")
    art.merge(d[["date", "cal", "binge", "fat_mass_lbs", "lm_lbs", "tdee"]], on="date").to_csv(ROOT / "analysis" / "DK_regime_daily.csv", index=False, float_format="%.2f")

    fm_min = w1["date"].iloc[w1["fat_mass_lbs"].idxmin()]; q_erupt = q.index[(q["binge"] >= 0.10).argmax()]
    post = w1[(w1["date"] > era_end) & (w1["date"] <= era_end + pd.Timedelta(days=180))]; pre_ = w1[w1["date"] <= era_end]
    pooled = [r for r in rows if r[0].startswith("pooled")]
    sig_u = [f"{r[1]} {r[2]:+.2f} [{r[3]:+.2f}, {r[4]:+.2f}]" for r in pooled if r[3] > 0 or r[4] < 0]
    sig_m = [f"{r[1]} {r[5]:+.2f} [{r[6]:+.2f}, {r[7]:+.2f}]" for r in pooled if r[6] > 0 or r[7] < 0]
    slow = {r[1]: r for r in pooled if r[1] in ("gap_frozen", "gap_ema100", "lm_deficit", "days_in_regime", "cum_deficit_Mcal")}
    print("\n=== Verdict ===")
    print(f"  Intake-only dating ends the restriction era at {sd.iloc[len(sd) // 2].date()} (single split; spread {(sd.max() - sd.min()).days} d over windows/transforms),")
    print(f"  HMM era exit {era_end.date()}, PELT {w1['date'][pelt(x1, 8 * sig2 * np.log(n))[0]].date()}. That is {(fm_min - era_end).days} d BEFORE the")
    print(f"  Kalman fat-mass minimum ({fm_min.date()}) and {(q_erupt.start_time - era_end).days} d before the first ≥ 10 % binge quarter ({q_erupt}).")
    print(f"  Sequence: intake stepped {post['cal'].mean() - pre_['cal'].mean():+.0f} kcal/day (era mean → next 180 d) while TDEE was still "
          f"{post['tdee'].mean():.0f}, so FM kept falling; binges (TDEE-relative AND absolute) erupted only after TDEE had dropped to "
          f"{q.loc[q_erupt, 'tdee']:.0f} ({q.loc[q_erupt, 'tdee'] - post['tdee'].mean():+.0f}) at the FM floor.")
    print("  Within the era every slow pressure (frozen gap, LM deficit, cumulative deficit, time) is monotone, so the")
    print("  one-event hazard cannot separate them from calendar time; the lead table shows how long each sat near its exit value.")
    print(f"  Pooled hazard over all completed pre-drug episodes, univariate β per SD [95 % CI]: "
          + "; ".join(f"{k} {r[2]:+.2f} [{r[3]:+.2f}, {r[4]:+.2f}]" for k, r in slow.items()))
    print(f"  Univariate CI excludes 0: {', '.join(sig_u) if sig_u else 'none'}.  Multivariate: {', '.join(sig_m) if sig_m else 'none'}.")
    print("  Artifact: analysis/DK_regime_daily.csv")


if __name__ == "__main__":
    main()
