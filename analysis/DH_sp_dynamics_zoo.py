#!/usr/bin/env python3
"""DH. Set point dynamics model zoo over the full 2011-2026 record.

Question: does the set point move down and chase current fat mass, or stay put?
The willpower era (2011-14: FM 100 → 15 lb on 1200-1900 cal/day) is where SP models
disagree most. If the SP froze or chased slowly, the subject spent 2012-13 under
1,000-3,000 cal/day of eating pressure held off by willpower, and the leak — binges —
should have built with the gap. If the SP chased fast (EMA, HL ≤ 100 d), pressure was
never more than a few hundred cal/day and the 2014-15 eruption needs another cause.

Every model produces a daily SP from the Kalman FM (or LM) path; gap = SP − FM.
Scores, all forward-looking where an outcome is involved:
  F-binge  r(gap_t, binge rate over t+1..t+90) and partial | trailing-90 binge rate   [pre-drug]
  F-intake partial r(gap_t, mean intake t+1..t+90 | past-90 intake, FM)              [pre-drug]
  Metab    partial r(gap_t, Kalman TDEE − expected RMR | FM)                          [pre-drug; CY caveat]
  Drug     AIC and P of  y = a + s·B + P·gap  on the 528 weekly-dosing days (DF design)
  Pressure mean gap in 2012-13 × 43 cal/lb = implied cal/day the subject overrode
  Timing   quarter of maximum gap (the eruption is 2014Q4 / 2015Q3; FM minimum 2014Q2)

Models
  frozen_init      SP = FM(2011-04)                         (never moves)
  frozen_max       SP = running max of FM                   (Lowe's weight-suppression reference)
  ema_H            SP = EMA(FM, HL = H d),  H ∈ {30, 50, 100, 200, 365, 730}
  asym_dA_uB       EMA with HL_down = A (when FM < SP), HL_up = B
  latch_subj       SmoothLatch (3 lb, 14 d, 0.01/d)         (BD/CF subject-fit)
  latch_trial      SmoothLatch (1 lb, 30 d, 0.0025/d)       (CF/BG trial-fit; DF)
  attractor_H      dSP = α(FM − SP) + β(SP_gen − SP), α: HL H d toward FM, β: HL 2 y toward SP_gen = FM(2011)
  lm_ema365        Dulloo protein-stat: gap = EMA(LM, 365 d) − LM   (FFM deficit drives appetite)
  lm_max           gap = running max of LM − LM
  floor_*          pressure gated by FM/(FM + 18 lb): leptin-floor reading of AC+BK

Appendix: two theories the project has data for but has not tested — Flatt's
glycogenostat (does yesterday's glycogen state predict today's intake?) and
seasonality of binges.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from numpy.linalg import lstsq
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent.parent
import sys
sys.path.insert(0, str(ROOT / "analysis"))
from intake_model import INTAKE_MULTIPLIER, describe as describe_intake_model  # single source of truth
from results import save as save_result, load as load_result
DRUG_START = "2024-09-17"
WEEKLY_END = "2026-03-31"
P_REF = 43.0          # cal/day per lb (DG, true cal) for the implied-pressure column
H = 90


def ema(x, hl):
    a = 1 - np.exp(-np.log(2) / hl)
    s = np.empty(len(x)); s[0] = x[0]
    for j in range(1, len(x)):
        s[j] = s[j - 1] + a * (x[j] - s[j - 1])
    return s


def asym_ema(x, hl_down, hl_up):
    ad, au = (1 - np.exp(-np.log(2) / hl_down)), (1 - np.exp(-np.log(2) / hl_up))
    s = np.empty(len(x)); s[0] = x[0]
    for j in range(1, len(x)):
        a = ad if x[j] < s[j - 1] else au
        s[j] = s[j - 1] + a * (x[j] - s[j - 1])
    return s


def latch(x, tol, hold, rate):
    s = np.empty(len(x)); s[0] = x[0]
    for i in range(1, len(x)):
        s[i] = s[i - 1]
        if i + 1 >= hold:
            w = x[i + 1 - hold:i + 1]; c = w.mean()
            if np.max(np.abs(w - c)) <= tol:
                s[i] = s[i - 1] + rate * (c - s[i - 1])
    return s


def attractor(x, hl_fm, hl_gen=730, sp_gen=None):
    a, b = (1 - np.exp(-np.log(2) / hl_fm)), (1 - np.exp(-np.log(2) / hl_gen))
    g = x[0] if sp_gen is None else sp_gen
    s = np.empty(len(x)); s[0] = x[0]
    for j in range(1, len(x)):
        s[j] = s[j - 1] + a * (x[j] - s[j - 1]) + b * (g - s[j - 1])
    return s


def partial_corr(y, x, controls):
    X = np.column_stack([np.ones(len(y))] + list(controls))
    ry = y - X @ lstsq(X, y, rcond=None)[0]
    rx = x - X @ lstsq(X, x, rcond=None)[0]
    return float(np.corrcoef(ry, rx)[0, 1])


def load():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    p1 = pd.read_csv(ROOT / "analysis" / "P1_smoothed_weight.csv", parse_dates=["date"])
    d = intake.merge(kalman[["date", "fat_mass_lbs", "lm_lbs", "tdee", "intake_true"]], on="date")
    d = d.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    d = d.merge(tirz[["date", "blood_level", "days_since_injection"]], on="date", how="left")
    d = d.dropna(subset=["fat_mass_lbs"]).sort_values("date").reset_index(drop=True)
    d["blood_level"] = d["blood_level"].fillna(0)
    d["binge"] = (d["calories"] > d["tdee"] + 1000).astype(float)
    pre = intake[(intake["date"] >= "2022-01-01") & (intake["date"] < DRUG_START)].copy()
    pre["dow"] = pre["date"].dt.day_of_week
    dow_adj = pre.groupby("dow")["calories"].mean() - pre["calories"].mean()
    d["cal_adj"] = d["calories"] - d["date"].dt.day_of_week.map(dow_adj)
    return d, p1


def build_models(fm, lm):
    fm_floor_gate = fm / (fm + 18.0)
    M = {
        "frozen_init": np.full(len(fm), fm[0]),
        "frozen_max": np.maximum.accumulate(fm),
    }
    for h in (30, 50, 100, 200, 365, 730):
        M[f"ema_{h}"] = ema(fm, h)
    M["asym_d365_u50"] = asym_ema(fm, 365, 50)
    M["asym_d50_u365"] = asym_ema(fm, 50, 365)
    M["latch_subj"] = latch(fm, 3.0, 14, 0.01)
    M["latch_trial"] = latch(fm, 1.0, 30, 0.0025)
    M["attractor_100"] = attractor(fm, 100)
    M["attractor_365"] = attractor(fm, 365)
    gaps = {k: v - fm for k, v in M.items()}
    gaps["lm_ema365"] = ema(lm, 365) - lm
    gaps["lm_max"] = np.maximum.accumulate(lm) - lm
    gaps["floor_frozen_max"] = gaps["frozen_max"] * fm_floor_gate
    gaps["floor_ema_100"] = gaps["ema_100"] * fm_floor_gate
    gaps["floor_latch_trial"] = gaps["latch_trial"] * fm_floor_gate
    return gaps


def main():
    d, p1 = load()
    fm, lm = d["fat_mass_lbs"].values, d["lm_lbs"].values
    gaps = build_models(fm, lm)
    pre = (d["date"] < DRUG_START).values
    cal = d["intake_true"].values     # consumed calories as used by P4
    binge = d["binge"].values
    fwd_binge = pd.Series(binge).rolling(H).mean().shift(-H).values
    bwd_binge = pd.Series(binge).rolling(H).mean().values
    fwd_cal = pd.Series(cal).rolling(H).mean().shift(-H).values
    bwd_cal = pd.Series(cal).rolling(H).mean().values
    tdee_res = (d["tdee"] - d["expected_rmr"]).values
    q = d["date"].dt.to_period("Q")
    drug = (d["date"] >= DRUG_START) & (d["date"] < WEEKLY_END) & (d["days_since_injection"] <= 6)
    drug = drug.values
    yd, Bd = d.loc[drug, "cal_adj"].values, d.loc[drug, "blood_level"].values
    era1213 = ((d["date"] >= "2012-01-01") & (d["date"] < "2014-01-01")).values

    print("=" * 118)
    print("DH. Set point dynamics zoo — full record, forward-looking scores (gap = SP − FM, + = FM below SP)")
    print("=" * 118)
    print(f"  pre-drug days {pre.sum()}, drug-era days {drug.sum()}, binge = logged cal > Kalman TDEE + 1000")
    hdr = (f"  {'model':<18} {'F-binge r':>9} {'|past':>7} {'F-intake|':>9} {'Metab|FM':>9} "
           f"{'drug AIC':>9} {'drug P':>7} {'drug s':>7} {'gap12-13':>9} {'cal/day':>8} {'max-gap Q':>10}")
    print(hdr)
    rows = []
    for name, gap in gaps.items():
        ok = pre & ~np.isnan(fwd_binge) & ~np.isnan(bwd_binge) & ~np.isnan(fwd_cal) & ~np.isnan(bwd_cal)
        r_fb = np.corrcoef(gap[ok], fwd_binge[ok])[0, 1]
        pr_fb = partial_corr(fwd_binge[ok], gap[ok], [bwd_binge[ok]])
        pr_fi = partial_corr(fwd_cal[ok], gap[ok], [bwd_cal[ok], fm[ok]])
        okm = pre & ~np.isnan(tdee_res)
        pr_m = partial_corr(tdee_res[okm], gap[okm], [fm[okm]])
        g = gap[drug]
        f = lambda p: p[0] + p[1] * Bd + p[2] * g
        r = least_squares(lambda p: yd - f(p), [2300, -80, 20], bounds=([1000, -600, -300], [3500, 0, 300]))
        sse = np.sum((yd - f(r.x)) ** 2)
        aic = len(yd) * np.log(sse / len(yd)) + 6
        gap_1213 = gap[era1213].mean()
        qmax = str(q[int(np.argmax(np.where(pre, gap, -np.inf)))])
        rows.append((name, r_fb, pr_fb, pr_fi, pr_m, aic, r.x[2], r.x[1], gap_1213, qmax))
        print(f"  {name:<18} {r_fb:>+9.3f} {pr_fb:>+7.3f} {pr_fi:>+9.3f} {pr_m:>+9.3f} "
              f"{aic:>9.1f} {r.x[2]:>+7.1f} {r.x[1]:>+7.1f} {gap_1213:>+9.1f} {gap_1213 * P_REF:>8.0f} {qmax:>10}")
    print("\n  F-binge r   : corr of gap with the NEXT 90-day binge rate; '|past' partials out the trailing 90-day rate.")
    print("  F-intake|   : partial corr with next-90-day intake, controlling past-90-day intake and FM. Lipostat: > 0.")
    print("  Metab|FM    : partial corr with Kalman TDEE − expected RMR controlling FM (CY attribution caveat).")
    print("  drug AIC/P  : DF design on the drug era with this model's gap. gap12-13: mean gap in 2012-13 and the")
    print("                eating pressure it implies at 43 cal/lb — the willpower the model requires.")

    # Quarterly picture: gap under three representative models vs binge rate, through the willpower era
    print("\n--- Quarterly: where did each model put the pressure, and when did binges come? ---")
    qd = pd.DataFrame({"q": q, "fm": fm, "binge": binge,
                       "frozen_max": gaps["frozen_max"], "ema_100": gaps["ema_100"], "ema_365": gaps["ema_365"],
                       "latch_trial": gaps["latch_trial"], "lm_ema365": gaps["lm_ema365"]})
    qa = qd.groupby("q").mean(numeric_only=True)
    qa["binge_next"] = qa["binge"].shift(-1)
    sub = qa.loc["2011Q2":"2016Q4"]
    print(f"  {'quarter':>8} {'FM':>6} {'binge%':>7} {'next-Q%':>8} {'frozen_max':>11} {'ema_100':>8} {'ema_365':>8} {'latch_trial':>12} {'lm_ema365':>10}")
    for qq, r_ in sub.iterrows():
        print(f"  {str(qq):>8} {r_['fm']:>6.1f} {r_['binge'] * 100:>7.1f} {r_['binge_next'] * 100:>8.1f} {r_['frozen_max']:>+11.1f} "
              f"{r_['ema_100']:>+8.1f} {r_['ema_365']:>+8.1f} {r_['latch_trial']:>+12.1f} {r_['lm_ema365']:>+10.1f}")
    pre_q = qa.loc[:"2024Q2"].dropna()
    print("\n  Quarterly corr(gap, NEXT quarter binge rate), 2011Q2-2024Q2:")
    for m_ in ("frozen_max", "ema_100", "ema_365", "latch_trial", "lm_ema365"):
        print(f"    {m_:<12} r = {np.corrcoef(pre_q[m_], pre_q['binge_next'])[0, 1]:+.3f}")

    # Stress tests on the two signals that stand out
    print("\n--- Stress tests: lean-mass deficit (Dulloo) and frozen FM gap ---")
    print("  NOTE: a frozen gap (SP₀ − FM) is collinear with −FM, so any partial 'controlling FM' is 0 by construction;")
    print("  its F-binge r is the statement 'lower FM → more binges', indistinguishable from a frozen set point.")
    strength = pd.read_csv(ROOT / "workout" / "strength.csv", parse_dates=["date"])
    ws = set(strength["date"].dt.normalize())
    wk30 = pd.Series([1.0 if x in ws else 0.0 for x in d["date"]]).rolling(30, min_periods=1).sum().values
    tdee = d["tdee"].values
    yr = d["date"].dt.year.values
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    lm_p2 = d[["date"]].merge(comp[["date", "ffm_lbs"]], on="date", how="left")["ffm_lbs"].values
    gap_lm_p2 = np.maximum.accumulate(np.nan_to_num(lm_p2, nan=np.nanmax(lm_p2))) - lm_p2
    tests = [("lm_max (Kalman LM)", gaps["lm_max"]), ("lm_max (scan-interpolated LM, no Kalman)", gap_lm_p2),
             ("lm_ema365", gaps["lm_ema365"]), ("ema_100 (FM)", gaps["ema_100"])]
    print(f"  {'model / window':<44} {'n':>5} {'F-int|past,FM':>13} {'+wk,TDEE':>9} {'+year FE':>9} {'F-binge|past':>12} {'+FM,yr':>8}")
    for name, gap in tests:
        for label, mask in [("pre-drug", pre), ("2017+", pre & (d["date"] >= "2017-01-01").values),
                            ("2014-16", (d["date"] >= "2014-01-01").values & (d["date"] < "2017-01-01").values)]:
            ok = mask & ~np.isnan(fwd_cal) & ~np.isnan(bwd_cal) & ~np.isnan(gap) & ~np.isnan(fwd_binge) & ~np.isnan(bwd_binge)
            yfe = [(yr[ok] == y_).astype(float) for y_ in np.unique(yr[ok])[1:]]
            base = [bwd_cal[ok], fm[ok]]
            r1 = partial_corr(fwd_cal[ok], gap[ok], base)
            r2 = partial_corr(fwd_cal[ok], gap[ok], base + [wk30[ok], tdee[ok]])
            r3 = partial_corr(fwd_cal[ok], gap[ok], base + [wk30[ok], tdee[ok]] + yfe)
            b1 = partial_corr(fwd_binge[ok], gap[ok], [bwd_binge[ok]])
            b3 = partial_corr(fwd_binge[ok], gap[ok], [bwd_binge[ok], fm[ok]] + yfe)
            print(f"  {name + ' / ' + label:<44} {ok.sum():>5} {r1:>+13.3f} {r2:>+9.3f} {r3:>+9.3f} {b1:>+12.3f} {b3:>+8.3f}")
    print("  Lower lean mass at a given fat mass and past intake predicts MORE intake over the next 90 days, in every")
    print("  window, with workout density, TDEE and year controlled, and with scan-only LM. This is the strongest")
    print("  forward appetite signal in the pre-drug record (vs ≤ +0.1 for any FM-based gap). It is what Dulloo's")
    print("  protein-stat / collateral-fattening model predicts; it is also what cumulative restriction history would")
    print("  look like, since LM-at-given-FM integrates past deficits. On-drug the LM gap is ~0 and the FM gap drives intake.")

    # Lean vs fat vs total weight as the level variable (DJ's caveat: LM | FM ≡ weight | FM)
    print("\n--- Which level carries the forward signal: lean, fat, or total weight? (next-90d intake ~ past-90d intake + level) ---")
    comp_fm = d[["date"]].merge(comp[["date", "fm_lbs"]], on="date", how="left")["fm_lbs"].values
    rng = np.random.default_rng(0)
    for src, fm_v, lm_v in (("Kalman", fm, lm), ("scan-interpolated (P2)", comp_fm, lm_p2)):
        ok = pre & ~np.isnan(fwd_cal) & ~np.isnan(bwd_cal) & ~np.isnan(fm_v) & ~np.isnan(lm_v)
        y, P_ = fwd_cal[ok], bwd_cal[ok]
        w_v = fm_v + lm_v
        yrs_ok = yr[ok]
        print(f"  {src} (n={ok.sum()}):")
        for label, cols in (("FM + LM, free", [fm_v[ok], lm_v[ok]]), ("weight only", [w_v[ok]]),
                            ("FM only", [fm_v[ok]]), ("LM only", [lm_v[ok]])):
            X = np.column_stack([np.ones(ok.sum()), P_] + cols)
            b = lstsq(X, y, rcond=None)[0]
            r2 = 1 - np.var(y - X @ b) / np.var(y)
            boots = []
            for _ in range(200):
                pick = rng.choice(np.unique(yrs_ok), len(np.unique(yrs_ok)))
                idx = np.concatenate([np.where(yrs_ok == p_)[0] for p_ in pick])
                boots.append(lstsq(X[idx], y[idx], rcond=None)[0][2:])
            lo, hi = np.percentile(boots, [2.5, 97.5], axis=0)
            coefs = "  ".join(f"{v:+.1f} [{l:+.0f},{h:+.0f}]" for v, l, h in zip(b[2:], lo, hi))
            print(f"    {label:16s} R² = {r2:.3f}   cal/day per lb: {coefs}")
    print("  Lean mass alone nearly doubles the R² of past intake; weight and fat mass add nothing; in the free fit")
    print("  the fat-mass coefficient is ~0 and the lean-mass coefficient is −20 to −36 cal/day per lb. It is lean mass.")

    # Appendix A: glycogenostat
    print("\n--- Appendix A. Flatt glycogenostat: yesterday's glycogen state → today's intake ---")
    gly = p1[["date", "glycogen_g"]].copy()
    dd = d[["date", "calories"]].merge(gly, on="date", how="inner")
    dd["cal_prev"] = dd["calories"].shift(1)
    dd["cal_prev7"] = dd["calories"].rolling(7).mean().shift(1)
    okg = dd.dropna()
    r_raw = np.corrcoef(okg["glycogen_g"], okg["calories"])[0, 1]
    r_part = partial_corr(okg["calories"].values, okg["glycogen_g"].values, [okg["cal_prev"].values, okg["cal_prev7"].values])
    print(f"  weigh-in days with P1 glycogen state: n = {len(okg)}; r(glycogen_lagged, today's intake) raw {r_raw:+.3f},")
    print(f"  partial | yesterday's and past-7-day intake {r_part:+.3f}. Glycogenostat predicts < 0 (depleted stores → eat more).")

    # Appendix B: seasonality
    print("\n--- Appendix B. Binge rate by month (pre-drug, 2014+) ---")
    s = d[pre & (d["date"] >= "2014-01-01").values]
    bym = s.groupby(s["date"].dt.month)["binge"].mean() * 100
    print("  " + "  ".join(f"{m:>2}:{v:4.1f}%" for m, v in bym.items()))
    print(f"  max/min month ratio {bym.max() / bym.min():.2f}; walks/sunlight peak in summer (AD/AK).")


if __name__ == "__main__":
    main()
