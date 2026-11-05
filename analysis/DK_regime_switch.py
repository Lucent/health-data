#!/usr/bin/env python3
"""DK. Behavioral regime switch: date the end of the willpower era objectively, then ask what predicted it.

DH showed no set point model puts its peak gap at the 2014Q3-Q4 binge eruption. The alternative is a
behavioral regime switch: the restriction regime ended. Here intake is modeled as a hidden regime process
and the switch is dated without reference to weight, binges or the subject's own epoch annotations.

Methods (all on log true intake = log(1.09 × logged kcal), fasting days floored at 200 kcal):
  HMM     Gaussian hidden Markov model, K = 2 and 3, EM (Baum-Welch) + Viterbi, on daily values and on
          non-overlapping 7- and 14-day block means. Restriction regime = lowest-mean state.
  PELT    exact penalized changepoint detection on the mean (L2 cost, min segment 30 d, penalty β·σ²·log n).
Switch dates are compared across methods/windows and against intake/diet_epochs.csv.

Switch hazard: on each restriction-regime day, y = 1 if the regime ends within 30 days. Discrete-time
logistic hazard on standardized predictors (univariate, ridge λ = 1e-3; multivariate, λ = 1), year-block
bootstrap CIs. Predictors: FM, LM deficit (running max LM − LM), days in regime, cumulative deficit,
trailing-30 d intake SD, trailing-30 d binge count, Kalman TDEE / expected RMR, season (cos, sin of day
of year), frozen-SP gap (FM₀ − FM) and EMA-100 d SP gap. A "lead" statistic says how long before the exit
each predictor had already reached 90 % of its exit value: a pressure that plateaued a year early cannot
date the switch. Secondary: same model on 2017-2024 (keto, weekend fasting, potato diets) and a pooled
hazard over every pre-drug restriction episode.

Command:  python analysis/DK_regime_switch.py
Artifact: analysis/DK_regime_daily.csv (date, true kcal, regime labels per method, hazard predictors)
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
MULT = 1.09
HORIZON = 30          # hazard horizon (days)
MIN_RUN = 28          # Viterbi runs shorter than this are merged into their neighbours
WIN1 = ("2011-04-20", "2017-12-31")
WIN2 = ("2017-01-01", "2024-09-16")   # pre-drug secondary window
RNG = np.random.default_rng(0)


# ----------------------------------------------------------------------------------------- models
def hmm_fit(x, K, iters=300, seed=0):
    """1-D Gaussian HMM by EM with scaled forward-backward. Returns dict with means, sds, A, loglik, path."""
    n = len(x); rng = np.random.default_rng(seed)
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
        mu = (g * x[:, None]).sum(0) / g.sum(0)
        var = np.maximum((g * (x[:, None] - mu) ** 2).sum(0) / g.sum(0), 1e-4)
        ll = np.log(c).sum()
        if ll - ll_old < 1e-6:
            break
        ll_old = ll
    # Viterbi
    lB = -0.5 * (x[:, None] - mu) ** 2 / var - 0.5 * np.log(2 * np.pi * var); lA = np.log(A)
    V = np.zeros((n, K)); bp = np.zeros((n, K), int); V[0] = np.log(pi) + lB[0]
    for t in range(1, n):
        s = V[t - 1][:, None] + lA; bp[t] = s.argmax(0); V[t] = s.max(0) + lB[t]
    path = np.zeros(n, int); path[-1] = V[-1].argmax()
    for t in range(n - 2, -1, -1):
        path[t] = bp[t + 1, path[t + 1]]
    order = np.argsort(mu); relabel = np.empty(K, int); relabel[order] = np.arange(K)
    return dict(mu=mu[order], sd=np.sqrt(var[order]), A=A[np.ix_(order, order)], ll=ll, path=relabel[path])


def hmm_best(x, K, restarts=4):
    fits = [hmm_fit(x, K, seed=s) for s in range(restarts)]
    return max(fits, key=lambda f: f["ll"])


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


def merge_short_runs(path, min_run):
    p = path.copy()
    changed = True
    while changed:
        changed = False
        runs = run_bounds(p)
        for i, (s, e, lab) in enumerate(runs):
            if e - s < min_run and len(runs) > 1:
                left = runs[i - 1][2] if i > 0 else None; right = runs[i + 1][2] if i + 1 < len(runs) else None
                p[s:e] = left if (right is None or (left is not None and runs[i - 1][1] - runs[i - 1][0] >= runs[i + 1][1] - runs[i + 1][0])) else right
                changed = True
                break
    return p


def run_bounds(p):
    cuts = np.flatnonzero(np.diff(p)) + 1
    starts = np.concatenate([[0], cuts]); ends = np.concatenate([cuts, [len(p)]])
    return [(s, e, int(p[s])) for s, e in zip(starts, ends)]


def logistic(X, y, lam=1e-3, iters=50):
    X1 = np.column_stack([np.ones(len(y)), X]); w = np.zeros(X1.shape[1]); R = np.eye(X1.shape[1]) * lam; R[0, 0] = 0
    for _ in range(iters):
        p = 1 / (1 + np.exp(-X1 @ w)); W = p * (1 - p)
        g = X1.T @ (y - p) - R @ w; Hs = (X1 * W[:, None]).T @ X1 + R
        step = np.linalg.solve(Hs, g); w += step
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
    d = intake[["date", "calories"]].merge(kal[["date", "fat_mass_lbs", "lm_lbs", "tdee"]], on="date")
    d = d.merge(comp[["date", "expected_rmr"]], on="date", how="left").dropna(subset=["fat_mass_lbs"])
    d = d.sort_values("date").reset_index(drop=True)
    d["cal"] = d["calories"] * MULT
    d["logcal"] = np.log(np.maximum(d["cal"], 200))
    d["binge"] = (d["calories"] > d["tdee"] + 1000).astype(float)
    d["balance"] = d["cal"] - d["tdee"]
    return d, ep


def block_means(x, b):
    nb = len(x) // b
    return x[:nb * b].reshape(nb, b).mean(1)


def fit_window(d, lo, hi, tag):
    """Fit every regime model on one window; return per-day DataFrame of labels and a table of switch dates."""
    w = d[(d["date"] >= lo) & (d["date"] <= hi)].reset_index(drop=True); x = w["logcal"].values; n = len(x)
    out = pd.DataFrame({"date": w["date"]}); fits = {}
    print(f"\n=== Window {tag}: {lo} .. {hi}, {n} days ===")
    print(f"  {'model':<16} {'K':>2} {'state means (true kcal/day)':<34} {'sd(log)':<20} {'dwell (d)':<22} loglik")
    for b, bname in ((1, "daily"), (7, "wk7"), (14, "wk14")):
        xb = block_means(x, b)
        for K in (2, 3):
            f = hmm_best(xb, K); fits[(bname, K)] = f
            lab = np.repeat(f["path"], b); lab = np.concatenate([lab, np.full(n - len(lab), lab[-1])])
            out[f"hmm_{bname}_K{K}"] = lab
            dwell = b / (1 - np.diag(f["A"]))
            print(f"  {'hmm_' + bname:<16} {K:>2} {str(np.round(np.exp(f['mu'])).astype(int)):<34} "
                  f"{str(np.round(f['sd'], 3)):<20} {str(np.round(dwell).astype(int)):<22} {f['ll']:.1f}")
    sig2 = np.median(np.abs(np.diff(x))) ** 2 / (0.6745 ** 2 * 2)
    for beta_mult in (2, 4, 8):
        cps = pelt(x, beta_mult * sig2 * np.log(n), 30)
        lab = np.zeros(n, int)
        for c in cps:
            lab[c:] += 1
        seg_mean = pd.Series(x).groupby(lab).transform("mean").values
        out[f"pelt_b{beta_mult}"] = (seg_mean > np.median(x)).astype(int)   # high/low vs window median
        out[f"pelt_b{beta_mult}_seg"] = lab
        print(f"  PELT β={beta_mult}·σ²·log n: {len(cps)} changepoints; segment means (kcal): "
              + ", ".join(f"{w['date'][c].date()}→{np.exp(x[c:(cps + [n])[i + 1]].mean()):.0f}" for i, c in enumerate(cps)))
    return out, fits


def main_switch(dates, lab, x, half=180):
    """Boundary with the largest upward step in mean log intake (±half days) among regime boundaries."""
    best = None
    for s, e, _ in run_bounds(lab)[1:]:
        if s < half or s + half > len(x):
            continue
        step = x[s:s + half].mean() - x[s - half:s].mean()
        if best is None or step > best[1]:
            best = (dates[s], step)
    return best


def restriction_segments(lab):
    """Restriction = state 0 after merging runs < MIN_RUN. Returns list of (start_idx, end_idx_exclusive)."""
    p = merge_short_runs((lab > 0).astype(int), MIN_RUN)
    return [(s, e) for s, e, l in run_bounds(p) if l == 0], p


def predictors(w, seg):
    s, e = seg; fm = w["fat_mass_lbs"].values; lm = w["lm_lbs"].values
    idx = np.arange(s, e)
    P = pd.DataFrame({
        "fm": fm[idx],
        "lm_deficit": (np.maximum.accumulate(lm) - lm)[idx],
        "days_in_regime": idx - s,
        "cum_deficit_Mcal": -np.cumsum(w["balance"].values[s:e]) / 1000,
        "intake_sd30": w["cal"].rolling(30, min_periods=10).std().values[idx],
        "binge30": w["binge"].rolling(30, min_periods=1).sum().values[idx],
        "tdee_ratio": (w["tdee"] / w["expected_rmr"]).values[idx],
        "season_cos": np.cos(2 * np.pi * w["date"].dt.dayofyear.values[idx] / 365.25),
        "season_sin": np.sin(2 * np.pi * w["date"].dt.dayofyear.values[idx] / 365.25),
        "gap_frozen": fm[0] - fm[idx],
        "gap_ema100": (ema(fm, 100) - fm)[idx],
    }, index=idx)
    return P


def hazard_table(P, y, years, label, boots=300):
    names = list(P.columns); Z = (P - P.mean()) / P.std(ddof=0); Z = Z.fillna(0).values
    print(f"\n--- Hazard: {label}: n = {len(y)} restriction days, {int(y.sum())} days with exit ≤ {HORIZON} d ---")
    print(f"  {'predictor':<18} {'β (per SD)':>10} {'95% CI (year-block boot)':>26} {'multi β':>9} {'95% CI':>18}")
    uy = np.unique(years)
    boot_u = {k: [] for k in names}; boot_m = {k: [] for k in names}
    w_multi = logistic(Z, y, lam=1.0)[1:]
    w_uni = [logistic(Z[:, [j]], y)[1] for j in range(len(names))]
    for _ in range(boots):
        pick = RNG.choice(uy, len(uy), replace=True); ii = np.concatenate([np.flatnonzero(years == p) for p in pick])
        if y[ii].sum() < 3 or y[ii].sum() > len(ii) - 3:
            continue
        wm = logistic(Z[ii], y[ii], lam=1.0)[1:]
        for j, k in enumerate(names):
            boot_m[k].append(wm[j]); boot_u[k].append(logistic(Z[ii][:, [j]], y[ii])[1])
    rows = []
    for j, k in enumerate(names):
        cu = np.percentile(boot_u[k], [2.5, 97.5]) if boot_u[k] else (np.nan, np.nan)
        cm = np.percentile(boot_m[k], [2.5, 97.5]) if boot_m[k] else (np.nan, np.nan)
        print(f"  {k:<18} {w_uni[j]:>+10.2f} {f'[{cu[0]:+.2f}, {cu[1]:+.2f}]':>26} {w_multi[j]:>+9.2f} {f'[{cm[0]:+.2f}, {cm[1]:+.2f}]':>18}")
        rows.append((k, w_uni[j], cu[0], cu[1], w_multi[j], cm[0], cm[1]))
    print(f"  ({len(boot_u[names[0]])} usable bootstrap draws; univariate ridge 1e-3, multivariate ridge 1 on standardized predictors)")
    return rows


def lead_table(P, label):
    print(f"\n--- Pressure lead, {label}: days before exit at which each predictor first reached 90 % of its exit value ---")
    print(f"  {'predictor':<18} {'at exit':>9} {'seg max':>9} {'exit pctile':>12} {'lead (d)':>9}  {'argmax (d before exit)':>22}")
    n = len(P)
    for k in P.columns:
        v = P[k].values; rng_ = v.max() - v.min()
        thr = v[-1] - 0.1 * rng_ if rng_ > 0 else v[-1]
        first = int(np.argmax(v >= thr)) if (v >= thr).any() else n - 1
        pct = (v < v[-1]).mean() * 100
        print(f"  {k:<18} {v[-1]:>9.2f} {v.max():>9.2f} {pct:>11.0f}% {n - 1 - first:>9d}  {n - 1 - int(np.argmax(v)):>22d}")


def main():
    d, ep = load()
    labs1, _ = fit_window(d, *WIN1, "primary 2011-2017")
    w1 = d[(d["date"] >= WIN1[0]) & (d["date"] <= WIN1[1])].reset_index(drop=True); x1 = w1["logcal"].values

    print("\n--- Main switch (largest upward 180 d step at a regime boundary) by method, 2011-2017 ---")
    method_cols = [c for c in labs1.columns if c.startswith("hmm_") or (c.startswith("pelt_") and not c.endswith("_seg"))]
    switches = {}
    for c in method_cols:
        ms = main_switch(w1["date"].values, labs1[c].values, x1)
        if ms:
            switches[c] = pd.Timestamp(ms[0])
            print(f"  {c:<16} {pd.Timestamp(ms[0]).date()}  step {np.exp(ms[1]):.3f}× (log {ms[1]:+.3f})")
    sw = pd.Series(switches); med = sw.sort_values().iloc[len(sw) // 2]
    print(f"  median switch date {med.date()}; range {sw.min().date()} .. {sw.max().date()} "
          f"(spread {(sw.max() - sw.min()).days} d across {len(sw)} method/window variants)")

    ref = "hmm_wk7_K2"
    segs, p_merged = restriction_segments(labs1[ref].values)
    print(f"\n--- Regimes under {ref} (runs < {MIN_RUN} d merged), 2011-2017 ---")
    print(f"  {'regime':<12} {'start':>10} {'end':>10} {'days':>5} {'kcal/day':>9} {'binge%':>7} {'TDEE':>6} {'FM in':>6} {'FM out':>7} {'LM in':>6} {'LM out':>7} {'bal/d':>6}")
    for s, e, l in run_bounds(p_merged):
        g = w1.iloc[s:e]
        print(f"  {'restrict' if l == 0 else 'free':<12} {g['date'].iloc[0].date()!s:>10} {g['date'].iloc[-1].date()!s:>10} {e - s:>5} "
              f"{g['cal'].mean():>9.0f} {g['binge'].mean() * 100:>7.1f} {g['tdee'].mean():>6.0f} {g['fat_mass_lbs'].iloc[0]:>6.1f} "
              f"{g['fat_mass_lbs'].iloc[-1]:>7.1f} {g['lm_lbs'].iloc[0]:>6.1f} {g['lm_lbs'].iloc[-1]:>7.1f} {g['balance'].mean():>+6.0f}")
    for K, c in ((2, "hmm_wk7_K2"), (3, "hmm_wk7_K3")):
        st = labs1[c].values
        print(f"  {c} state means: " + "; ".join(
            f"state {k}: {w1['cal'][st == k].mean():.0f} kcal, binge {w1['binge'][st == k].mean() * 100:.1f}%, "
            f"TDEE {w1['tdee'][st == k].mean():.0f}, bal {w1['balance'][st == k].mean():+.0f}, {int((st == k).sum())} d" for k in range(K)))

    print("\n--- HMM/PELT boundaries vs hand-annotated diet epochs (2011-2017) ---")
    ep1 = ep[(ep["start"] >= WIN1[0]) & (ep["start"] <= WIN1[1])]
    bounds = [w1["date"].iloc[s] for s, e, l in run_bounds(p_merged)[1:]]
    for _, r in ep1.iterrows():
        near = min(bounds, key=lambda b: abs((b - r["start"]).days)) if bounds else None
        print(f"  epoch {r['label']:<26} starts {r['start'].date()}  nearest {ref} boundary {near.date() if near is not None else '-'} "
              f"({(near - r['start']).days:+d} d)" if near is not None else "")
    for b in bounds:
        near = ep1.iloc[(ep1["start"] - b).abs().dt.days.argmin()]
        print(f"  {ref} boundary {b.date()} → nearest epoch start {near['start'].date()} ({near['label']}, {(near['start'] - b).days:+d} d)")

    # ---------------- hazard, primary window
    rows_all = []
    main_seg = max(segs, key=lambda se: se[1] - se[0])
    P = predictors(w1, main_seg)
    lead_table(P, f"longest restriction segment {w1['date'].iloc[main_seg[0]].date()} .. {w1['date'].iloc[main_seg[1] - 1].date()}")
    Ps, ys = [], []
    for s, e in segs:
        if e >= len(w1):     # censored at window end: drop days whose horizon runs past the data
            e_obs = e - HORIZON; pp = predictors(w1, (s, e)).iloc[:max(e_obs - s, 0)]; yy = np.zeros(len(pp))
        else:
            pp = predictors(w1, (s, e)); yy = ((e - 1 - pp.index.values) < HORIZON).astype(float)
        Ps.append(pp); ys.append(yy)
    P1 = pd.concat(Ps); y1 = np.concatenate(ys); yrs1 = w1["date"].dt.year.values[P1.index.values]
    rows_all += [("2011-2017",) + r for r in hazard_table(P1, y1, yrs1, f"2011-2017, {len(segs)} restriction segment(s)")]

    # ---------------- secondary window and pooled hazard
    labs2, _ = fit_window(d, *WIN2, "secondary 2017-2024 (pre-drug)")
    w2 = d[(d["date"] >= WIN2[0]) & (d["date"] <= WIN2[1])].reset_index(drop=True)
    segs2, p2 = restriction_segments(labs2[ref].values)
    print(f"\n--- Restriction episodes under {ref}, 2017-2024 (runs < {MIN_RUN} d merged) ---")
    print(f"  {'start':>10} {'end':>10} {'days':>5} {'kcal/day':>9} {'binge%':>7} {'bal/d':>6} {'FM in':>6} {'FM out':>7} {'LM def out':>10} {'gap100 out':>10}  annotated epoch at exit")
    ep2 = ep[(ep["end"] >= WIN2[0]) & (ep["start"] <= WIN2[1])]
    for s, e in segs2:
        g = w2.iloc[s:e]; pe = predictors(w2, (s, e)); d_exit = g["date"].iloc[-1]
        hit = ep2[(ep2["start"] <= d_exit) & (ep2["end"] >= d_exit)]["label"].tolist()
        print(f"  {g['date'].iloc[0].date()!s:>10} {d_exit.date()!s:>10} {e - s:>5} {g['cal'].mean():>9.0f} {g['binge'].mean() * 100:>7.1f} "
              f"{g['balance'].mean():>+6.0f} {g['fat_mass_lbs'].iloc[0]:>6.1f} {g['fat_mass_lbs'].iloc[-1]:>7.1f} "
              f"{pe['lm_deficit'].iloc[-1]:>10.1f} {pe['gap_ema100'].iloc[-1]:>+10.1f}  {', '.join(hit) or '-'}")
    for _, r in ep2[ep2["label"].str.contains("keto|fasting|potato|restriction|protein")].iterrows():
        frac = (p2[(w2["date"] >= r["start"]) & (w2["date"] <= r["end"])] == 0).mean()
        print(f"  epoch {r['label']:<24} {r['start'].date()}..{r['end'].date()}: {frac * 100:3.0f}% of days in HMM restriction state")
    Ps2, ys2 = [], []
    for s, e in segs2:
        if e >= len(w2):
            continue
        pp = predictors(w2, (s, e)); Ps2.append(pp); ys2.append(((e - 1 - pp.index.values) < HORIZON).astype(float))
    if Ps2:
        P2 = pd.concat(Ps2); y2 = np.concatenate(ys2); yrs2 = w2["date"].dt.year.values[P2.index.values]
        rows_all += [("2017-2024",) + r for r in hazard_table(P2, y2, yrs2, f"2017-2024, {len(Ps2)} completed restriction episodes")]
        Pp = pd.concat([P1, P2]); yp = np.concatenate([y1, y2]); yrp = np.concatenate([yrs1, yrs2])
        rows_all += [("pooled",) + r for r in hazard_table(Pp, yp, yrp, "pooled 2011-2024 pre-drug episodes")]

    # ---------------- artifact + verdict
    art = pd.concat([labs1.assign(window="2011-2017"), labs2.assign(window="2017-2024")])
    art = art.merge(d[["date", "cal", "binge", "fat_mass_lbs", "lm_lbs", "tdee"]], on="date", how="left")
    art.to_csv(ROOT / "analysis" / "DK_regime_daily.csv", index=False, float_format="%.3f")
    pd.DataFrame(rows_all, columns=["window", "predictor", "beta_uni", "lo_uni", "hi_uni", "beta_multi", "lo_multi", "hi_multi"]) \
        .to_csv(ROOT / "analysis" / "DK_hazard_coefficients.csv", index=False, float_format="%.3f")

    print("\n=== Verdict ===")
    print(f"  The restriction regime ends at {med.date()} (method spread {(sw.max() - sw.min()).days} d). In the primary window every")
    print("  slow biological pressure (frozen gap, LM deficit, cumulative deficit, time in regime) is a monotone function of")
    print("  calendar time inside a single 3-year segment, so a one-event hazard cannot separate them: see the lead table —")
    print("  a predictor that sat within 10 % of its exit value for months before the exit did not time the switch.")
    print("  Fast predictors (gap_ema100, trailing binge count, intake SD, season) are the ones that can rise INTO the")
    print("  switch; the pooled multi-episode hazard says which of them actually did, across every pre-drug episode.")
    print("  Artifacts: analysis/DK_regime_daily.csv, analysis/DK_hazard_coefficients.csv")


if __name__ == "__main__":
    main()
