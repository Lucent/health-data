#!/usr/bin/env python3
"""DJ. Protein-stat (Dulloo) vs cumulative-restriction confound for the lean-mass appetite signal.

DH found that, pre-tirzepatide, lower lean mass at a given fat mass and past intake predicts
HIGHER intake over the next 90 d (partial r ≈ +0.56 Kalman LM, +0.44 scan-interpolated LM).
Note: the running max of LM is the 2011 BOD POD value (159.3 lb) for the whole record, so
"LM deficit = max(LM) − LM" is just a constant minus LM. The signal is LM *level* | FM.

Two readings:  (a) protein-stat — an FFM deficit drives hyperphagia until FFM is restored;
               (b) confound — low LM | FM marks long past restriction, and post-diet rebound
                   raises intake regardless of LM.
Discriminators run here (all outcomes strictly forward of the predictor date; year-block
bootstrap 95% CIs; Kalman LM, scan-interpolated P2 LM and raw scans in parallel):
  1. Restriction-history controls (past-180/365-d intake, past-365-d energy balance, year FE).
  2. Training instrument: trailing workout count → LM (first stage); 2SLS LM slope vs OLS.
     Protein-stat: training-built LM lowers future intake like diet-built LM. Confound: no effect.
  3. Timescale: LM reference (max / EMA 180-365-730 / last scan) × horizon (30/90/180 d).
  4. Thermogenesis: LM deficit vs Kalman TDEE − expected RMR, and vs calorimetry RMR residual.
  5. Protein seeking: next-90-d protein % of calories and g/day vs LM deficit.
  6. Binge rate: next-90-d binge rate vs LM deficit | trailing rate, FM.
  7. Drug era: DF design with LM deficit vs frozen-FM gap; pre-drug slope × observed ΔLM.
  8. Scan-level: 49 pre-drug composition scans, next-90-d intake ~ scan LM + scan FM + past intake.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.linalg import lstsq
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent.parent
MULT = 1.09
DRUG_START = pd.Timestamp("2024-09-17")
WEEKLY_END = "2026-03-31"
TRAIN_END = "2025-08-19"
NBOOT = 400
rng = np.random.default_rng(7)


def ema(x, hl):
    a = 1 - np.exp(-np.log(2) / hl)
    s = np.empty(len(x)); s[0] = x[0]
    for j in range(1, len(x)):
        s[j] = s[j - 1] + a * (x[j] - s[j - 1])
    return s


def ols(y, cols):
    X = np.column_stack([np.ones(len(y))] + list(cols))
    return lstsq(X, y, rcond=None)[0]


def pcorr(y, x, ctrls):
    X = np.column_stack([np.ones(len(y))] + list(ctrls))
    ry = y - X @ lstsq(X, y, rcond=None)[0]
    rx = x - X @ lstsq(X, x, rcond=None)[0]
    return float(np.corrcoef(ry, rx)[0, 1])


def boot(stat, years, nb=NBOOT):
    """Year-block bootstrap: resample calendar years with replacement, recompute stat(idx)."""
    uy = np.unique(years); out = []
    for _ in range(nb):
        pick = rng.choice(uy, len(uy), replace=True)
        idx = np.concatenate([np.where(years == u)[0] for u in pick])
        try:
            out.append(stat(idx))
        except Exception:
            pass
    return np.nanpercentile(out, [2.5, 97.5])


def fmt(v, ci, d=3):
    return f"{v:+.{d}f} [{ci[0]:+.{d}f}, {ci[1]:+.{d}f}]"


def fwd(x, h): return pd.Series(x).rolling(h).mean().shift(-h).values   # t+1..t+h
def bwd(x, h): return pd.Series(x).rolling(h).mean().values             # t-h+1..t


def load():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kal = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    p2 = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "composition" / "composition.csv", parse_dates=["date"])
    comp["lm_scan"] = comp["weight_lbs"] - comp["fat_mass_lbs"]
    comp = comp.groupby("date", as_index=False)[["lm_scan", "fat_mass_lbs"]].mean()
    strength = pd.read_csv(ROOT / "workout" / "strength.csv", parse_dates=["date"])
    d = intake.merge(kal[["date", "fat_mass_lbs", "lm_lbs", "tdee"]], on="date")
    d = d.merge(p2[["date", "ffm_lbs", "expected_rmr"]], on="date", how="left")
    d = d.merge(tirz[["date", "blood_level", "days_since_injection", "dose_mg"]], on="date", how="left")
    d = d.merge(comp[["date", "lm_scan"]], on="date", how="left")
    d = d.dropna(subset=["fat_mass_lbs", "lm_lbs"]).sort_values("date").reset_index(drop=True)
    d["lm_scan_ff"] = d["lm_scan"].ffill()           # last scan carried forward: no look-ahead, no Kalman
    d["blood_level"] = d["blood_level"].fillna(0)
    d["workout"] = d["date"].isin(set(strength["date"].dt.normalize())).astype(float)
    pre = intake[(intake["date"] >= "2022-01-01") & (intake["date"] < DRUG_START)].copy()
    dow_adj = pre.groupby(pre["date"].dt.day_of_week)["calories"].mean() - pre["calories"].mean()
    d["cal_adj"] = d["calories"] - d["date"].dt.day_of_week.map(dow_adj)
    d["cal"] = d["calories"] * MULT
    d["binge"] = (d["calories"] > d["tdee"] + 1000).astype(float)
    return d, comp


def main():
    d, comp = load()
    fm, tdee, cal, yr = d["fat_mass_lbs"].values, d["tdee"].values, d["cal"].values, d["date"].dt.year.values
    LM = {"kalman": d["lm_lbs"].values, "p2": d["ffm_lbs"].values, "scan_ff": d["lm_scan_ff"].values}
    wk = {n: pd.Series(d["workout"]).rolling(n, min_periods=1).sum().values for n in (30, 60, 90, 180)}
    F, Bk = {h: fwd(cal, h) for h in (30, 90, 180, 365)}, {h: bwd(cal, h) for h in (30, 90, 180, 365)}
    eb365 = bwd(cal - tdee, 365)
    premask = lambda h: (d["date"] < DRUG_START - pd.Timedelta(days=h)).values

    def ok_of(h, *arrs):
        m = premask(h) & ~np.isnan(F[h]) & ~np.isnan(Bk[h])
        for a in arrs: m &= ~np.isnan(a)
        return m

    print("=" * 100)
    print("DJ. Protein-stat vs cumulative-restriction confound — pre-tirzepatide record, forward windows")
    print("=" * 100)
    print(f"  running max of Kalman LM = {np.maximum.accumulate(LM['kalman']).max():.1f} lb on 100% of days → "
          "'LM deficit' ≡ −LM; all tests below use LM level (slope sign flipped: protein-stat ⇒ negative).")

    # ---- 1. Baseline + restriction-history controls -------------------------------------------
    print("\n--- 1. Next-90-d intake ~ LM + FM + past-90 intake, then adding restriction-history controls ---")
    print(f"  {'LM series':<9} {'controls':<34} {'n':>5} {'partial r [95% CI]':>26} {'cal/day per lb LM [95% CI]':>30}")
    base = {}
    for name, lm in LM.items():
        for lab, extra in [("FM, past90", []), ("+past180, past365", [Bk[180], Bk[365]]),
                           ("+past365 energy balance", [eb365]), ("+past365, EB365, wk90, TDEE", [Bk[365], eb365, wk[90], tdee]),
                           ("+year FE", "FE")]:
            ok = ok_of(90, lm, Bk[365], eb365)
            y, x, c = F[90][ok], lm[ok], [fm[ok], Bk[90][ok]]
            if extra == "FE":
                c = c + [(yr[ok] == u).astype(float) for u in np.unique(yr[ok])[1:]]
            else:
                c = c + [e[ok] for e in extra]
            pr = lambda i, y=y, x=x, c=c: pcorr(y[i], x[i], [a[i] for a in c])
            sl = lambda i, y=y, x=x, c=c: ols(y[i], [x[i]] + [a[i] for a in c])[1]
            I = np.arange(len(y))
            r, b = pr(I), sl(I)
            rci, bci = boot(pr, yr[ok]), boot(sl, yr[ok])
            if lab == "FM, past90": base[name] = (b, bci)
            print(f"  {name:<9} {lab:<34} {ok.sum():>5} {fmt(r, rci):>26} {fmt(b, bci, 1):>30}")
    print("  Confound predicts the LM slope collapses once past restriction is controlled; protein-stat predicts it survives.")
    print("\n  1b. Split by LM direction over the past 90 d (rebound confound ⇒ effect concentrated while LM is falling):")
    for name in ("kalman", "p2"):
        lm = LM[name]; dlm = pd.Series(lm).diff(90).values
        for lab, m in [("LM falling", dlm < 0), ("LM rising", dlm >= 0)]:
            ok = ok_of(90, lm, dlm) & m
            y, x, c = F[90][ok], lm[ok], [fm[ok], Bk[90][ok]]
            pr = lambda i, y=y, x=x, c=c: pcorr(y[i], x[i], [a[i] for a in c])
            print(f"     {name:<7} {lab:<11} n={ok.sum():>4}  partial r {fmt(pr(np.arange(len(y))), boot(pr, yr[ok]))}")

    # ---- 2. Training instrument ---------------------------------------------------------------
    print("\n--- 2. Training instrument: trailing workouts (and workouts × deficit, CB's lean-sparing channel) → LM ---")
    print("  first stage LM ~ z + deficitN + FM + past90; strength = partial r(LM_hat, LM | controls); 2SLS cal/day per lb")
    print(f"  {'LM':<7} {'win':>4} {'instrument':<12} {'fs strength':>11} {'fs lb/unit [CI]':>26} {'reduced cal/d/unit [CI]':>26} "
          f"{'2SLS cal/d/lb [CI]':>24} {'OLS':>6} {'OLS+wk':>7}")
    for name in ("kalman", "p2"):
        lm = LM[name]
        for n in (60, 90, 180):
            dfc = np.clip(bwd(tdee - cal, n), 0, None) / 1000.0        # mean deficit over window, kcal/day
            for zlab, z_ in [("wk", wk[n]), ("wk×deficit", wk[n] * dfc)]:
                ok = ok_of(90, lm, dfc)
                y, x, z, c = F[90][ok], lm[ok], z_[ok], [dfc[ok], fm[ok], Bk[90][ok]]
                fs = lambda i: ols(x[i], [z[i]] + [a[i] for a in c])[1]
                rf = lambda i: ols(y[i], [z[i]] + [a[i] for a in c])[1]
                def tsls(i):
                    X = [z[i]] + [a[i] for a in c]
                    xhat = np.column_stack([np.ones(len(i))] + X) @ ols(x[i], X)
                    return ols(y[i], [xhat] + [a[i] for a in c])[1]
                I = np.arange(len(y))
                X = [z] + c; xhat = np.column_stack([np.ones(len(y))] + X) @ ols(x, X)
                print(f"  {name:<7} {n:>4} {zlab:<12} {pcorr(x, xhat, c):>+11.3f} {fmt(fs(I), boot(fs, yr[ok]), 3):>26} "
                      f"{fmt(rf(I), boot(rf, yr[ok]), 1):>26} {fmt(tsls(I), boot(tsls, yr[ok]), 0):>24} "
                      f"{ols(y, [x] + c)[1]:>+6.0f} {ols(y, [x, z] + c)[1]:>+7.0f}")
    print("  Protein-stat: 2SLS negative, comparable to OLS. Confound: 2SLS ≈ 0. Neither is readable if fs strength ≲ 0.2:")
    print("  the instrument is then too weak (CB/CC: workouts spare lean for ~1 day; there is no chronic training-built stock).")

    # ---- 3. Timescale sweep -------------------------------------------------------------------
    print("\n--- 3. Timescale sweep: partial r(ref − LM, next-H intake | FM, past-H intake) ---")
    rows = []
    for name, lm in LM.items():
        refs = {"max": np.full(len(lm), np.nanmax(lm)), "ema180": None, "ema365": None, "ema730": None}
        lmf = pd.Series(lm).ffill().bfill().values
        for hl in (180, 365, 730): refs[f"ema{hl}"] = ema(lmf, hl)
        for rname, ref in refs.items():
            gap = ref - lm
            row = {"lm": name, "ref": rname}
            for h in (30, 90, 180, 365):
                ok = ok_of(h, gap)
                row[f"r{h}"] = pcorr(F[h][ok], gap[ok], [fm[ok], Bk[h][ok]])
            rows.append(row)
    sw = pd.DataFrame(rows)
    print(sw.to_string(index=False, float_format=lambda v: f"{v:+.3f}"))
    st = sw.set_index(["lm", "ref"]).stack()
    print(f"  peak: {st.idxmax()}  r = {st.max():+.3f}.  The fixed reference (LM level) beats every adaptive EMA reference;")
    print("  the signal grows with horizon, i.e. LM tracks a slow component of appetite, not a 30-day rebound.")
    sw.to_csv(ROOT / "analysis" / "DJ_proteinstat_sweep.csv", index=False, float_format="%.4f")

    # ---- 4. Thermogenesis ---------------------------------------------------------------------
    print("\n--- 4. Does low LM also suppress expenditure? (Dulloo: FFM deficit lowers thermogenesis) ---")
    res = (d["tdee"] - d["expected_rmr"]).values
    for name in ("kalman", "p2"):
        lm = LM[name]; ok = premask(0) & ~np.isnan(res) & ~np.isnan(lm) & ~np.isnan(Bk[90])
        y, x, c = res[ok], lm[ok], [fm[ok], Bk[90][ok]]
        pr = lambda i: pcorr(y[i], x[i], [a[i] for a in c]); sl = lambda i: ols(y[i], [x[i]] + [a[i] for a in c])[1]
        I = np.arange(len(y))
        print(f"  Kalman TDEE − expected RMR ~ LM_{name} | FM, past90: partial r {fmt(pr(I), boot(pr, yr[ok]))}, "
              f"slope {fmt(sl(I), boot(sl, yr[ok]), 1)} cal/day per lb (Dulloo ⇒ positive: more LM, more heat)")
    rmr = pd.read_csv(ROOT / "RMR" / "rmr.csv", parse_dates=["date"])
    rm = rmr.merge(d[["date", "fat_mass_lbs", "lm_lbs", "ffm_lbs", "expected_rmr"]], on="date")
    rm = rm[rm["date"] < DRUG_START]
    rres = (rm["rmr_kcal"] - rm["expected_rmr"]).values
    for col in ("lm_lbs", "ffm_lbs"):
        r_ = pcorr(rres, rm[col].values, [rm["fat_mass_lbs"].values])
        b_ = ols(rres, [rm[col].values, rm["fat_mass_lbs"].values])[1]
        print(f"  calorimetry RMR − expected RMR (n={len(rm)} pre-drug dates) ~ {col} | FM: partial r {r_:+.3f}, slope {b_:+.1f} cal/day/lb")
    print("  (expected_rmr is itself a×FFM + b×FM + c from P2, so the residual is what LM explains beyond its fitted RMR share.)")

    # ---- 5. Protein seeking -------------------------------------------------------------------
    print("\n--- 5. Protein seeking: next-90-d protein share and g/day vs LM | FM, past intake, past protein ---")
    pg = d["protein_g"].values.astype(float); cr = d["calories"].values.astype(float)
    pf_f = pd.Series(pg * 4).rolling(90).sum().shift(-90).values / pd.Series(cr).rolling(90).sum().shift(-90).values
    pf_b = pd.Series(pg * 4).rolling(90).sum().values / pd.Series(cr).rolling(90).sum().values
    pg_f, pg_b = fwd(pg, 90), bwd(pg, 90)
    for name in ("kalman", "p2"):
        lm = LM[name]; ok = ok_of(90, lm, pf_f, pf_b)
        for lab, y, c in [("protein % of cal", pf_f, [pf_b]), ("protein g/day", pg_f, [pg_b])]:
            yy, x, cc = y[ok], lm[ok], [fm[ok], Bk[90][ok]] + [a[ok] for a in c]
            pr = lambda i, yy=yy, x=x, cc=cc: pcorr(yy[i], x[i], [a[i] for a in cc])
            print(f"  {name:<7} {lab:<18} partial r {fmt(pr(np.arange(len(yy))), boot(pr, yr[ok]))}  (protein-stat ⇒ negative: low LM → seek protein)")

    # ---- 6. Binge -----------------------------------------------------------------------------
    print("\n--- 6. Next-90-d binge rate (cal > TDEE + 1000) vs LM | trailing-90 binge rate, FM ---")
    bf, bb = fwd(d["binge"].values, 90), bwd(d["binge"].values, 90)
    for name in ("kalman", "p2"):
        lm = LM[name]; ok = ok_of(90, lm, bf, bb)
        y, x, c = bf[ok], lm[ok], [bb[ok], fm[ok]]
        pr = lambda i: pcorr(y[i], x[i], [a[i] for a in c]); sl = lambda i: ols(y[i], [x[i]] + [a[i] for a in c])[1] * 100
        I = np.arange(len(y))
        print(f"  {name:<7} partial r {fmt(pr(I), boot(pr, yr[ok]))}, slope {fmt(sl(I), boot(sl, yr[ok]), 2)} binge-%-points per lb LM")

    # ---- 7. Drug era --------------------------------------------------------------------------
    print("\n--- 7. Drug era (DF design, weekday-adjusted intake on weekly-dosing days): LM deficit vs frozen-FM gap ---")
    g = d[(d["date"] >= DRUG_START) & (d["date"] < WEEKLY_END) & (d["days_since_injection"] <= 6)].reset_index(drop=True)
    y, B = g["cal_adj"].values, g["blood_level"].values
    gaps = {"none": np.zeros(len(g)), "FM0 − FM": g["fat_mass_lbs"].iloc[0] - g["fat_mass_lbs"].values,
            "LM0 − LM (Kalman)": g["lm_lbs"].iloc[0] - g["lm_lbs"].values, "LM0 − LM (P2)": g["ffm_lbs"].iloc[0] - g["ffm_lbs"].values}
    tr, te = np.where(g["date"] < TRAIN_END)[0], np.where(g["date"] >= TRAIN_END)[0]
    print(f"  n={len(y)}; FM {g['fat_mass_lbs'].iloc[0]:.1f}→{g['fat_mass_lbs'].iloc[-1]:.1f}, LM Kalman {g['lm_lbs'].iloc[0]:.1f}→{g['lm_lbs'].iloc[-1]:.1f}, "
          f"P2 {g['ffm_lbs'].iloc[0]:.1f}→{g['ffm_lbs'].iloc[-1]:.1f} lb")
    print(f"  {'gap term':<20} {'P cal/day/lb':>13} {'s':>7} {'AIC':>8} {'ΔAIC vs none':>13} {'OOS RMSE 12.5mg':>16} {'OOS bias':>9}")
    for lab, gp in gaps.items():
        k = 2 if lab == "none" else 3
        f = (lambda p, gp=gp: p[0] + p[1] * B + (p[2] * gp if len(p) > 2 else 0))
        x0, lb, ub = ([2300, -80, 20][:k], [1000, -600, -600][:k], [3500, 0, 600][:k])
        fit = lambda idx: least_squares(lambda p: y[idx] - f(p)[idx], x0, bounds=(lb, ub)).x
        pa = fit(np.arange(len(y))); sse = np.sum((y - f(pa)) ** 2); aic = len(y) * np.log(sse / len(y)) + 2 * k
        if lab == "none": aic0 = aic
        pt = fit(tr); pred = f(pt)
        print(f"  {lab:<20} {(pa[2] if k == 3 else 0):>+13.1f} {pa[1]:>+7.1f} {aic:>8.1f} {aic - aic0:>+13.1f} "
              f"{np.sqrt(np.mean((y[te] - pred[te]) ** 2)):>16.0f} {np.mean(y[te] - pred[te]):>+9.0f}")
    b_lm = base["kalman"][0]
    first8, last8 = g.iloc[:56]["cal_adj"].mean(), g.iloc[-56:]["cal_adj"].mean()
    print(f"  Pre-drug LM slope {b_lm:+.1f} cal/day/lb × drug-era ΔLM {g['lm_lbs'].iloc[-1] - g['lm_lbs'].iloc[0]:+.1f} lb ⇒ "
          f"{b_lm * (g['lm_lbs'].iloc[-1] - g['lm_lbs'].iloc[0]):+.0f} cal/day predicted; observed first-8-wk → last-8-wk {first8:.0f} → {last8:.0f} ({last8 - first8:+.0f})")
    print("  The LM deficit has no room to move on-drug; whatever drives the on-drug intake rise is not lean mass.")

    # ---- 8. Scan-level ------------------------------------------------------------------------
    print("\n--- 8. Raw composition scans: next-90-d intake ~ scan LM + scan FM + past-90 intake (no Kalman, no interpolation) ---")
    sc = comp.merge(d[["date"]].assign(fwd90=F[90], bwd90=Bk[90], wk90=wk[90]), on="date")
    sc = sc[(sc["date"] < DRUG_START - pd.Timedelta(days=90)) & sc["fwd90"].notna() & sc["bwd90"].notna()].reset_index(drop=True)
    sy = sc["date"].dt.year.values
    y, x, f_, b_, w_ = sc["fwd90"].values, sc["lm_scan"].values, sc["fat_mass_lbs"].values, sc["bwd90"].values, sc["wk90"].values
    for lab, c in [("| FM, past90", [f_, b_]), ("| FM, past90, wk90", [f_, b_, w_])]:
        pr = lambda i, c=c: pcorr(y[i], x[i], [a[i] for a in c]); sl = lambda i, c=c: ols(y[i], [x[i]] + [a[i] for a in c])[1]
        I = np.arange(len(y))
        print(f"  n={len(sc)} scans ({sy.min()}–{sy.max()}) {lab:<22} partial r {fmt(pr(I), boot(pr, sy))}, slope {fmt(sl(I), boot(sl, sy), 1)} cal/day/lb")
    thin = sc[sc["date"].diff().dt.days.fillna(999) >= 30].reset_index(drop=True)
    r_t = pcorr(thin["fwd90"].values, thin["lm_scan"].values, [thin["fat_mass_lbs"].values, thin["bwd90"].values])
    print(f"  thinned to ≥30 d apart: n={len(thin)}, partial r {r_t:+.3f}")

    print("\n=== Verdict ===")
    print("  Read sections 1, 2 and 8 together: if the LM slope survives restriction-history controls, the 2SLS")
    print("  (training-built LM) slope is negative and comparable to OLS, and raw scans reproduce it, the protein-stat")
    print("  reading stands. If the slope dies under past-365 controls or 2SLS ≈ 0, it is cumulative restriction.")


if __name__ == "__main__":
    main()
