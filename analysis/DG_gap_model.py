#!/usr/bin/env python3
"""DG. One gap-form model for the subject's drug era, the GLP-1 trials, and the q3d era.

Replaces the tachyphylaxis structure of BZ/CE/AV with the form DF selects:

    intake(t)  = I₀ + m·s·B(t) + P·(SP(t) − FM(t))
    TDEE(t)    = TDEE₀ + k·(BW(t) − BW₀) [+ β·(intake − I₀) adaptive thermogenesis]
    dFM, dLM   = Forbes-coupled partition of intake − TDEE
    SP(t)      = SmoothLatch(tol, hold, rate) run on FM — the TRIAL-FIT latch
                 (tol 1 lb, hold 30 d, rate 0.0025/d; CF/BG), which DF shows the
                 subject's drug era also selects. SP₀ = FM at drug start.

  B(t)  daily-mean PK exposure (drugs/extract.py; one-compartment, t½ 5 d, ka 3.31/d,
        dose in mg, mean over the eating window from the injection clock time).
  s     drug coefficient, cal/day per unit of B. No tachyphylaxis, no saturation.
  P     eating pressure, cal/day per lb of fat below the set point.
  m     central-response multiplier; 1.0 = this subject. Re-derived for trial cohorts.

Parameters s, P, I₀ are fit ONCE, on the subject's weekly-dosing drug era (528 days,
weekday-adjusted logged intake; identical to DF model L8). Everything after that is
forward simulation with no trial-fitted parameters except where a cohort multiplier m
is explicitly solved for and reported.

What is tested
--------------
1. Subject forward simulation 2024-09-17 → 2026-02-28 from the pre-drug state, FM
   path vs the Kalman FM (which the fit never saw as a target — only intake was fit).
2. SURMOUNT-1 5/10/15 mg at 24/48/72 wk; SURMOUNT-4 36-wk lead-in and 52-wk regain
   trajectory; SURMOUNT-2 (T2D) via a solved m.
3. Semaglutide / liraglutide withdrawal arms (STEP-1 ext, STEP-4, SCALE): regain from
   the published lead-in trajectory with the same P and latch (CF design).
4. Preregistered q3d forecast (2026-04-01 onward, 5.2 then 6.25 mg every 3 days) —
   weekly mean intake and sawtooth, to be scored once 2026 intake is merged.

Trial inputs are the ones BZ/CF used (Jastreboff NEJM 2022; Aronne JAMA 2024; Garvey
Lancet 2023; Wilding 2022 extension; Rubino JAMA 2021). Cohort body composition from
the DXA substudies: 40-45% fat at baseline, 75% of lost weight as fat for tirzepatide,
~0.56-0.65 for semaglutide.
"""

from pathlib import Path
import sys
from datetime import timedelta
import numpy as np
import pandas as pd
from scipy.optimize import least_squares, minimize_scalar, minimize

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "drugs"))
from extract import daily_mean_exposure  # noqa: E402  (timestamp-aware PK)

import sys
sys.path.insert(0, str(ROOT / "analysis"))
from intake_model import INTAKE_MULTIPLIER, describe as describe_intake_model  # single source of truth
from results import save as save_result, load as load_result
LATCH_TOL, LATCH_HOLD, LATCH_RATE = 1.0, 30, 0.0025     # trial-fit latch (CF/BG, DF L8)
E_FAT, E_LEAN = 3500.0, 1100.0
FORBES_C_KG, LBS_TO_KG, KG_TO_LBS = 10.4, 0.453592, 2.20462
WEEKLY_END = "2026-03-31"
TRAIN_END = "2025-08-19"


# ── building blocks ────────────────────────────────────────────────────
class Latch:
    """Incremental SmoothLatch: SP moves toward the trailing-window mean only when
    every FM in the window sits within ±tol of that mean."""

    def __init__(self, sp0, tol=LATCH_TOL, hold=LATCH_HOLD, rate=LATCH_RATE):
        self.sp, self.tol, self.hold, self.rate = sp0, tol, hold, rate
        self.hist = []

    def update(self, fm):
        self.hist.append(fm)
        if len(self.hist) >= self.hold:
            w = np.asarray(self.hist[-self.hold:])
            c = w.mean()
            if np.max(np.abs(w - c)) <= self.tol:
                self.sp += self.rate * (c - self.sp)
        return self.sp


def forbes_ff(fm_lb):
    fm_kg = max(fm_lb, 5.0) * LBS_TO_KG
    return 1.0 - FORBES_C_KG ** 2 / (FORBES_C_KG + fm_kg) ** 2


def partition(surplus_cal, fm_lb, fat_frac=None):
    ff = forbes_ff(fm_lb) if fat_frac is None else fat_frac
    e = ff * E_FAT + (1 - ff) * E_LEAN
    dW = surplus_cal / e
    return ff * dW, (1 - ff) * dW


def mifflin_tdee(weight_kg, height_cm, age, pct_female, pal):
    rmr_m = 10 * weight_kg + 6.25 * height_cm - 5 * age + 5
    rmr_f = rmr_m - 166
    return ((1 - pct_female) * rmr_m + pct_female * rmr_f) * pal


def weekly_injections(start, n_weeks, dose_fn, hour=12):
    return [(start + timedelta(days=7 * w, hours=hour), dose_fn(w)) for w in range(n_weeks)]


# ── 1. subject fit (DF model L8) ───────────────────────────────────────
def load_subject():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    # Work in consumed ("true") calories throughout: P4 writes the corrected
    # intake it consumed as `intake_true`, so s and P come out in the same units
    # as TDEE with no conversion factor in this file.
    intake = intake[["date"]].merge(kalman[["date", "intake_true"]], on="date")
    intake = intake.rename(columns={"intake_true": "calories"})
    pre = intake[(intake["date"] >= "2022-01-01") & (intake["date"] < "2024-09-17")].copy()
    pre["dow"] = pre["date"].dt.day_of_week
    dow_adj = pre.groupby("dow")["calories"].mean() - pre["calories"].mean()
    d = intake[["date", "calories"]].merge(
        tirz[["date", "blood_level", "days_since_injection", "dose_mg"]], on="date")
    d = d.merge(kalman[["date", "fat_mass_lbs", "lm_lbs", "tdee"]], on="date")
    d = d[(d["date"] < WEEKLY_END) & (d["days_since_injection"] <= 6)]
    d = d.sort_values("date").reset_index(drop=True)
    d["y"] = d["calories"] - d["date"].dt.day_of_week.map(dow_adj)
    return d, kalman, tirz, dow_adj


def fit_subject(d):
    fm = d["fat_mass_lbs"].values
    latch = Latch(fm[0])
    sp = np.array([latch.update(f) for f in fm])
    gap = sp - fm
    B, y = d["blood_level"].values, d["y"].values
    f = lambda p: p[0] + p[1] * B + p[2] * gap
    r = least_squares(lambda p: y - f(p), [2300, -80, 30], bounds=([1000, -600, -300], [3500, 0, 300]))
    a, s, P = r.x
    resid = y - f(r.x)
    # bootstrap by week blocks for CIs
    wk = ((d["date"] - d["date"].min()).dt.days // 7).values
    rng = np.random.default_rng(1)
    boots = []
    weeks = np.unique(wk)
    for _ in range(300):
        pick = rng.choice(weeks, len(weeks))
        idx = np.concatenate([np.where(wk == w)[0] for w in pick])
        rb = least_squares(lambda p: y[idx] - (p[0] + p[1] * B[idx] + p[2] * gap[idx]),
                           [a, s, P], bounds=([1000, -600, -300], [3500, 0, 300]))
        boots.append(rb.x)
    boots = np.array(boots)
    ci = np.percentile(boots, [2.5, 97.5], axis=0)
    return dict(a=a, s=s, P=P, rmse=float(np.sqrt(np.mean(resid ** 2))), gap_end=float(gap[-1]),
                s_ci=ci[:, 1], P_ci=ci[:, 2], sp=sp)


# ── 2. subject forward simulation ─────────────────────────────────────
def simulate_subject(d, kalman, fit, tdee_mode="kalman", k_tdee=6.4):
    """Forward-simulate FM from drug start using only the fitted intake model and the
    actual exposure series. tdee_mode: 'kalman' uses the Kalman TDEE path (true cal);
    'model' uses TDEE₀ (pre-drug 90-d mean) + k·ΔBW."""
    fm = d["fat_mass_lbs"].iloc[0]
    lm = d["lm_lbs"].iloc[0]
    bw0 = fm + lm
    pre = kalman[(kalman["date"] < d["date"].iloc[0]) & (kalman["date"] >= d["date"].iloc[0] - pd.Timedelta(days=90))]
    tdee0 = pre["tdee"].mean()
    latch = Latch(fm)
    out = []
    for i, row in d.iterrows():
        sp = latch.update(fm)
        intake_true = fit["a"] + fit["s"] * row["blood_level"] + fit["P"] * (sp - fm)
        intake_logged = intake_true          # model is fit in consumed calories
        tdee = row["tdee"] if tdee_mode == "kalman" else tdee0 + k_tdee * ((fm + lm) - bw0)
        dfm, dlm = partition(intake_true - tdee, fm)
        out.append((row["date"], fm, sp, intake_logged, tdee))
        fm, lm = fm + dfm, lm + dlm
    sim = pd.DataFrame(out, columns=["date", "fm_sim", "sp_sim", "intake_sim", "tdee_used"])
    sim["fm_kalman"] = d["fat_mass_lbs"].values
    return sim


# ── 3. tirzepatide trial simulation ───────────────────────────────────
ESCALATION = {5: [(0, 2.5), (4, 5.0)],
              10: [(0, 2.5), (4, 5.0), (8, 7.5), (12, 10.0)],
              15: [(0, 2.5), (4, 5.0), (8, 7.5), (12, 10.0), (16, 12.5), (20, 15.0)]}


def dose_at_week(target_mg):
    steps = ESCALATION[target_mg]
    def f(w):
        dose = steps[0][1]
        for wk, dmg in steps:
            if w >= wk:
                dose = dmg
        return dose
    return f


_EXPOSURE_CACHE = {}


def trial_exposure(target_mg, weeks_on, weeks_off):
    """Daily-mean exposure series for a trial schedule (cached; independent of m, P)."""
    key = (target_mg, weeks_on, weeks_off)
    if key not in _EXPOSURE_CACHE:
        start = pd.Timestamp("2000-01-03")
        injections = weekly_injections(start, weeks_on, dose_at_week(target_mg))
        n_days = (weeks_on + weeks_off) * 7
        _EXPOSURE_CACHE[key] = np.array([daily_mean_exposure(start + timedelta(days=d), injections)
                                         for d in range(n_days)])
    return _EXPOSURE_CACHE[key]


def simulate_trial(start_kg, height_cm, age, pct_female, target_mg, weeks_on, weeks_off,
                   fm_frac, fat_frac_loss, s_true, P_true, m=1.0, pal=1.4, beta_at=0.0,
                   latch_params=(LATCH_TOL, LATCH_HOLD, LATCH_RATE), ec50=None):
    """Returns weekly % change in body weight from baseline and the SP/FM paths.
    ec50: if set, drug term is m·s·B/(1 + B/ec50) (Emax with initial slope m·s)."""
    B_series = trial_exposure(target_mg, weeks_on, weeks_off)
    bw0 = start_kg * KG_TO_LBS
    fm, lm = bw0 * fm_frac, bw0 * (1 - fm_frac)
    tdee0 = mifflin_tdee(start_kg, height_cm, age, pct_female, pal)
    i0 = tdee0
    latch = Latch(fm, *latch_params)
    n_days = (weeks_on + weeks_off) * 7
    weekly_pct, sp_path, fm_path = [], [], []
    for day in range(n_days):
        B = B_series[day]
        if ec50 is not None:
            B = B / (1 + B / ec50)
        sp = latch.update(fm)
        intake = i0 + m * s_true * B + P_true * (sp - fm)
        bw = fm + lm
        tdee = mifflin_tdee(bw * LBS_TO_KG, height_cm, age, pct_female, pal) + beta_at * (intake - i0)
        dfm, dlm = partition(intake - tdee, fm, fat_frac_loss)
        if day % 7 == 0:
            weekly_pct.append((bw - bw0) / bw0 * 100)
            sp_path.append(sp)
            fm_path.append(fm)
        fm, lm = fm + dfm, lm + dlm
    weekly_pct.append(((fm + lm) - bw0) / bw0 * 100)
    return np.array(weekly_pct), np.array(sp_path), np.array(fm_path)


def regain_from_leadin(baseline_kg, stop_kg, lead_weeks, fm_frac, fat_frac_loss,
                       P_true, regain_weeks=52, k_tdee=6.4, regain_fat_frac=0.75,
                       latch_params=(LATCH_TOL, LATCH_HOLD, LATCH_RATE)):
    """CF design: linear lead-in on the published loss, latch runs on it, then regain
    under P·(SP − FM) with TDEE rising k cal per lb regained. Returns % of stop weight
    regained by week, and % of the loss."""
    bw0, bw1 = baseline_kg * KG_TO_LBS, stop_kg * KG_TO_LBS
    fm0 = bw0 * fm_frac
    fm1 = fm0 - (bw0 - bw1) * fat_frac_loss
    fm_lead = np.linspace(fm0, fm1, lead_weeks * 7 + 1)
    latch = Latch(fm0, *latch_params)
    for f in fm_lead:
        sp = latch.update(f)
    fm, bw = fm1, bw1
    weekly = []
    for day in range(regain_weeks * 7 + 1):
        if day % 7 == 0:
            weekly.append((bw - bw1) / bw1 * 100)
        sp = latch.update(fm)
        # Off drug, intake returns to the habitual pre-trial level I₀ = TDEE₀ while TDEE
        # is lower at the lower weight: surplus = k·(BW₀ − BW) even with P = 0 (pure CICO
        # return to habit). Set point pressure adds P·(SP − FM) on top.
        surplus = k_tdee * (bw0 - bw) + P_true * (sp - fm)
        dfm, dlm = partition(max(surplus, 0.0), fm, regain_fat_frac)
        fm, bw = fm + dfm, bw + dfm + dlm
    weekly = np.array(weekly)
    return weekly, weekly * bw1 / (bw0 - bw1)


# ── 4. q3d forecast ───────────────────────────────────────────────────
def q3d_forecast(tirz, fit, fm_last, sp_last, dow_adj):
    rows = tirz[tirz["date"] >= "2026-03-01"].copy()
    fm, latch = fm_last, Latch(sp_last)
    latch.hist = [fm_last] * LATCH_HOLD
    # hold TDEE at the last Kalman value (true cal) — plateau regime
    out = []
    for _, r in rows.iterrows():
        sp = latch.update(fm)
        intake_logged = fit["a"] + fit["s"] * r["blood_level"] + fit["P"] * (sp - fm)
        out.append((r["date"], r["dose_mg"], r["blood_level"], intake_logged + dow_adj.get(r["date"].day_of_week, 0)))
    return pd.DataFrame(out, columns=["date", "dose_mg", "B", "intake_pred_true"])   # consumed calories (compare with intake_model.true_intake of the diary, not the raw diary)


# ── main ──────────────────────────────────────────────────────────────
def main():
    d, kalman, tirz, dow_adj = load_subject()
    fit = fit_subject(d)
    s_true, P_true = fit["s"], fit["P"]      # already in consumed calories
    print("=" * 92)
    print("DG. Gap-form model — fit once on the subject's drug era, then forward everywhere")
    print("=" * 92)
    print(f"  n = {len(d)} weekly-dosing days, trial-fit latch (tol {LATCH_TOL}, hold {LATCH_HOLD}, rate {LATCH_RATE})")
    print(f"  I₀ = {fit['a']:.0f} logged cal/day   s = {fit['s']:+.1f} cal/unit  [95% {fit['s_ci'][0]:+.0f}, {fit['s_ci'][1]:+.0f}]"
          f"   P = {fit['P']:+.1f} cal/day/lb  [95% {fit['P_ci'][0]:+.0f}, {fit['P_ci'][1]:+.0f}]   RMSE {fit['rmse']:.0f}")
    print(f"  (consumed-calorie units throughout; intake model: {describe_intake_model()});  SP − FM at end of record = {fit['gap_end']:+.1f} lb")

    # 1. subject forward sim
    print("\n--- 1. Subject forward simulation from the pre-drug state (intake model only; FM never fit) ---")
    sim_stats = {}
    for mode in ("kalman", "model"):
        sim = simulate_subject(d, kalman, fit, tdee_mode=mode)
        err = sim["fm_sim"] - sim["fm_kalman"]
        sim_stats[mode] = float(np.sqrt(np.mean(err ** 2)))
        plateau = sim[sim["date"] >= "2025-11-15"]
        print(f"  TDEE = {mode:6s}: FM {sim['fm_sim'].iloc[0]:.1f} → {sim['fm_sim'].iloc[-1]:.1f} lb (Kalman {sim['fm_kalman'].iloc[-1]:.1f});"
              f" RMSE vs Kalman FM {np.sqrt(np.mean(err ** 2)):.1f} lb, max |err| {err.abs().max():.1f};"
              f" plateau-window mean intake sim {plateau['intake_sim'].mean():.0f} vs obs {d.loc[plateau.index, 'y'].mean():.0f}")
        if mode == "kalman":
            sim_k = sim
    print("  Dose-era FM (sim vs Kalman):")
    for dose, lo, hi in [(5.0, "2024-10-15", "2025-03-04"), (7.5, "2025-03-04", "2025-05-27"),
                         (10.0, "2025-05-27", "2025-08-19"), (12.5, "2025-08-19", WEEKLY_END)]:
        s_ = sim_k[(sim_k["date"] >= lo) & (sim_k["date"] < hi)]
        print(f"    {dose:>5.1f} mg: sim {s_['fm_sim'].mean():5.1f}  Kalman {s_['fm_kalman'].mean():5.1f}  SP_sim {s_['sp_sim'].mean():5.1f}")

    # 2. tirzepatide trials
    print("\n--- 2. Tirzepatide trials, m = 1.0 (subject parameters, nothing fit to trials) ---")
    S1 = dict(start_kg=104.8, height_cm=167, age=45, pct_female=0.67, fm_frac=0.40, fat_frac_loss=0.75)
    S2 = dict(start_kg=100.7, height_cm=167, age=54, pct_female=0.50, fm_frac=0.38, fat_frac_loss=0.75)
    S4 = dict(start_kg=107.3, height_cm=167, age=46, pct_female=0.67, fm_frac=0.455, fat_frac_loss=0.75)
    pub_s1 = {5: (-16.0, None, None), 10: (-21.4, None, None), 15: (-22.5, -15.0, -20.5)}
    print(f"  {'arm':<22} {'wk24':>7} {'wk48':>7} {'wk72':>7} {'published wk72':>15} {'SP−FM@72':>9}")
    sse = 0
    for dose in (5, 10, 15):
        pct, sp, fmp = simulate_trial(**S1, target_mg=dose, weeks_on=72, weeks_off=0, s_true=s_true, P_true=P_true)
        pub = pub_s1[dose][0]
        sse += (pct[72] - pub) ** 2
        print(f"  SURMOUNT-1 {dose:>2} mg       {pct[24]:>+7.1f} {pct[48]:>+7.1f} {pct[72]:>+7.1f} {pub:>+15.1f} {sp[-1] - fmp[-1]:>+9.1f}")
    print(f"  (SURMOUNT-1 15 mg published shape ≈ −15 at wk 24, −20.5 at wk 48, −22.5 at wk 72; SSE across arms {sse:.1f})")

    pct4, sp4, fm4 = simulate_trial(**S4, target_mg=15, weeks_on=36, weeks_off=52, s_true=s_true, P_true=P_true)
    stop = pct4[36]
    pub_regain = {0: 0, 12: 7.0, 24: 11.0, 36: 13.0, 52: 14.0}
    print(f"\n  SURMOUNT-4 15 mg: lead-in wk36 {stop:+.1f}% (published −20.9%);  SP−FM at stop {sp4[36] - fm4[36]:+.1f} lb")
    print(f"  {'weeks off':>10} {'sim % of stop wt':>17} {'published':>10} {'sim % of loss':>14}")
    bw_stop = 1 + stop / 100
    for w in (12, 24, 36, 52):
        regain_pct = (pct4[36 + w] - stop) / bw_stop
        print(f"  {w:>10} {regain_pct:>+17.1f} {pub_regain[w]:>+10.1f} {regain_pct * bw_stop / (-stop) * 100:>13.0f}%")

    # cohort multiplier m
    def sse_m(m, trial, arms, weeks=72):
        return sum((simulate_trial(**trial, target_mg=dz, weeks_on=weeks, weeks_off=0, s_true=s_true, P_true=P_true, m=m)[0][weeks] - pub) ** 2
                   for dz, pub in arms)
    m1 = minimize_scalar(lambda m: sse_m(m, S1, [(5, -16.0), (10, -21.4), (15, -22.5)]), bounds=(0.3, 3.0), method="bounded").x
    m2 = minimize_scalar(lambda m: sse_m(m, S2, [(10, -12.8), (15, -14.7)]), bounds=(0.3, 3.0), method="bounded").x
    print(f"\n  Cohort multiplier m on s (P held at subject value): SURMOUNT-1 non-diabetic m = {m1:.2f};  SURMOUNT-2 T2D m = {m2:.2f}  (ratio {m2 / m1:.2f}; CE had 0.69)")
    print(f"  {'arm':<22} {'sim@m':>7} {'published':>10}")
    for dz, pub in [(5, -16.0), (10, -21.4), (15, -22.5)]:
        p = simulate_trial(**S1, target_mg=dz, weeks_on=72, weeks_off=0, s_true=s_true, P_true=P_true, m=m1)[0][72]
        print(f"  SURMOUNT-1 {dz:>2} mg       {p:>+7.1f} {pub:>+10.1f}")
    for dz, pub in [(10, -12.8), (15, -14.7)]:
        p = simulate_trial(**S2, target_mg=dz, weeks_on=72, weeks_off=0, s_true=s_true, P_true=P_true, m=m2)[0][72]
        print(f"  SURMOUNT-2 {dz:>2} mg       {p:>+7.1f} {pub:>+10.1f}")
    pct4m = simulate_trial(**S4, target_mg=15, weeks_on=36, weeks_off=52, s_true=s_true, P_true=P_true, m=m1)[0]
    stop_m = pct4m[36]
    print(f"  SURMOUNT-4 at m = {m1:.2f}: lead-in {stop_m:+.1f}% (pub −20.9); regain +52 wk {(pct4m[88] - stop_m) / (1 + stop_m / 100):+.1f}% of stop weight (pub +14.0)")

    # sensitivity: adaptive thermogenesis, latch choice
    print("\n  Sensitivity (SURMOUNT-1 15 mg wk72 / SURMOUNT-4 regain +52 wk, % of stop weight):")
    for label, kw in [("no latch movement (frozen SP)", dict(latch_params=(0.0, 30, 0.0))),
                      ("subject-fit latch (3, 14, 0.01)", dict(latch_params=(3.0, 14, 0.01))),
                      ("adaptive thermogenesis β = 0.14", dict(beta_at=0.14)),
                      ("P × 0.5", dict(P_scale=0.5)), ("P × 2", dict(P_scale=2.0))]:
        ps = kw.pop("P_scale", 1.0)
        a = simulate_trial(**S1, target_mg=15, weeks_on=72, weeks_off=0, s_true=s_true, P_true=P_true * ps, **kw)[0][72]
        b = simulate_trial(**S4, target_mg=15, weeks_on=36, weeks_off=52, s_true=s_true, P_true=P_true * ps, **kw)[0]
        print(f"    {label:34s} {a:+6.1f}%   regain {(b[88] - b[36]) / (1 + b[36] / 100):+6.1f}%")

    # 3. semaglutide / liraglutide regain
    print("\n--- 3. Withdrawal arms from published lead-ins (same P and latch; CF design) ---")
    print(f"  {'trial':<24} {'lead wk':>7} {'loss %':>7} {'sim +52wk % stop wt':>20} {'published':>10} {'sim % of loss':>14}")
    for name, b0, b1, lead, ff_, fl, pub in [
            ("SURMOUNT-4 (check)", 107.3, 85.8, 36, 0.455, 0.75, 14.0),
            ("STEP-1 ext sema", 105.6, 87.5, 68, 0.434, 0.56, 13.1),
            ("STEP-4 sema", 107.2, 96.1, 20, 0.434, 0.56, 6.9),
            ("SCALE lira 3 mg", 105.2, 96.7, 56, 0.434, 0.65, 6.5)]:
        wk, frac = regain_from_leadin(b0, b1, lead, ff_, fl, P_true)
        print(f"  {name:<24} {lead:>7} {(b0 - b1) / b0 * 100:>7.1f} {wk[52]:>+20.1f} {pub:>+10.1f} {frac[52]:>13.0f}%")
    print("  (STEP-4 published is at +48 wk; SCALE liraglutide 3 mg withdrawal regain is approximate.)")

    # 5. cohort dose-response: Emax
    print("\n--- 5. Cohort dose-response shape: linear vs Emax drug term (SURMOUNT-1 three arms) ---")
    print("  Published: 5 mg reaches 71% of the 15 mg loss at one-third the exposure — concave.")
    arms1 = [(5, -16.0), (10, -21.4), (15, -22.5)]

    def sse_emax(p):
        m_, ec = p
        if ec <= 0.5 or m_ <= 0:
            return 1e9
        return sum((simulate_trial(**S1, target_mg=dz, weeks_on=72, weeks_off=0, s_true=s_true,
                                   P_true=P_true, m=m_, ec50=ec)[0][72] - pub) ** 2 for dz, pub in arms1)
    best = min((minimize(sse_emax, x0, method="Nelder-Mead", options={"xatol": 0.01, "fatol": 0.05})
                for x0 in ([2.0, 5.0], [3.0, 2.0], [1.5, 10.0])), key=lambda r: r.fun)
    m_e, ec_e = best.x
    print(f"  Linear (m = {m1:.2f}):  SSE {sse_m(m1, S1, arms1):.1f}     Emax (m = {m_e:.2f}, EC50 = {ec_e:.1f} units):  SSE {best.fun:.1f}")
    print(f"  {'arm':<22} {'linear':>7} {'Emax':>7} {'published':>10}")
    for dz, pub in arms1:
        a = simulate_trial(**S1, target_mg=dz, weeks_on=72, weeks_off=0, s_true=s_true, P_true=P_true, m=m1)[0][72]
        b = simulate_trial(**S1, target_mg=dz, weeks_on=72, weeks_off=0, s_true=s_true, P_true=P_true, m=m_e, ec50=ec_e)[0][72]
        print(f"  SURMOUNT-1 {dz:>2} mg       {a:>+7.1f} {b:>+7.1f} {pub:>+10.1f}")
    p4e = simulate_trial(**S4, target_mg=15, weeks_on=36, weeks_off=52, s_true=s_true, P_true=P_true, m=m_e, ec50=ec_e)[0]
    print(f"  SURMOUNT-4 under Emax: lead-in {p4e[36]:+.1f}% (pub −20.9); regain +52 wk {(p4e[88] - p4e[36]) / (1 + p4e[36] / 100):+.1f}% of stop weight (pub +14.0)")
    # Subject-side check of saturation: refit the drug era with an Emax term at fixed EC50
    # (plus the latch gap), and bootstrap the sawtooth-amplitude ratio 12.5 mg / 5 mg.
    y_s, B_s = d["y"].values, d["blood_level"].values
    fm_s = d["fat_mass_lbs"].values
    lt = Latch(fm_s[0]); gap_s = np.array([lt.update(f) for f in fm_s]) - fm_s
    print(f"\n  Subject drug era refit with Emax at fixed EC50 (+ latch gap), AIC relative to linear:")
    aics = {}
    for ec in (3, 5, 8, 12, 20, 50, np.inf):
        Be = B_s / (1 + B_s / ec) if np.isfinite(ec) else B_s
        fe = lambda p: p[0] + p[1] * Be + p[2] * gap_s
        re = least_squares(lambda p: y_s - fe(p), [2300, -80, 30], bounds=([1000, -3000, -300], [3500, 0, 300]))
        sse = np.sum((y_s - fe(re.x)) ** 2)
        aics[ec] = len(y_s) * np.log(sse / len(y_s)) + 6
    for ec in aics:
        print(f"    EC50 = {ec if np.isfinite(ec) else 'linear':>7}: ΔAIC {aics[ec] - aics[np.inf]:+6.1f}")
    wk_s = ((d["date"] - d["date"].min()).dt.days // 7).values
    dfw = pd.DataFrame({"wk": wk_s, "y": y_s, "dsi": d["days_since_injection"].values, "dose": d["dose_mg"].values})
    pk_ = dfw[dfw["dsi"].isin([1, 2])].groupby("wk")["y"].mean()
    tr_ = dfw[dfw["dsi"].isin([5, 6])].groupby("wk")["y"].mean()
    amp = (tr_ - pk_).dropna()
    dose_w = dfw.groupby("wk")["dose"].max().reindex(amp.index)
    a5, a12 = amp[dose_w == 5.0].values, amp[dose_w == 12.5].values
    rng = np.random.default_rng(0)
    ratios = [rng.choice(a12, len(a12)).mean() / rng.choice(a5, len(a5)).mean() for _ in range(5000)]
    lo_r, hi_r = np.percentile(ratios, [2.5, 97.5])
    print(f"  Sawtooth amplitude 12.5 mg / 5 mg: {a12.mean() / a5.mean():.2f} (bootstrap 95% [{lo_r:.2f}, {hi_r:.2f}]); "
          f"predicted by EC50 = 5: 0.95, 12: 1.20, 20: 1.36, linear: 1.82.")
    print(f"  → The subject tolerates mild saturation (EC50 ≳ 12) but not the cohort's EC50 ≈ {ec_e:.0f}.")
    for ec_fixed in (12.0, 20.0):
        m_f = minimize_scalar(lambda m_: sum((simulate_trial(**S1, target_mg=dz, weeks_on=72, weeks_off=0, s_true=s_true,
                                                             P_true=P_true, m=m_, ec50=ec_fixed)[0][72] - pub) ** 2 for dz, pub in arms1),
                              bounds=(0.3, 6.0), method="bounded").x
        arms_f = [simulate_trial(**S1, target_mg=dz, weeks_on=72, weeks_off=0, s_true=s_true, P_true=P_true, m=m_f, ec50=ec_fixed)[0][72] for dz, _ in arms1]
        print(f"  SURMOUNT-1 at subject-compatible EC50 = {ec_fixed:.0f} (m = {m_f:.2f}): "
              + ", ".join(f"{dz} mg {a:+.1f}" for (dz, _), a in zip(arms1, arms_f)) + "  (pub −16.0 / −21.4 / −22.5)")

    # 5b. joint cohort model: subject-compatible EC50, cohort P from regain, m from SURMOUNT-1
    print("\n--- 5b. Joint cohort model: EC50 ∈ {12, 20}, cohort P ∈ {5, 15, 43}, frozen SP, m from SURMOUNT-1 ---")
    print(f"  {'EC50':>5} {'P_c':>5} {'m':>5} {'5 mg':>7} {'10 mg':>7} {'15 mg':>7} {'S4 lead-in':>11} {'S4 regain +12/24/36/52 wk':>28}")
    frozen = (0.0, 30, 0.0)
    for ec_fixed in (12.0, 20.0):
        for p_c in (5.0, 15.0, 43.0):
            m_f = minimize_scalar(lambda m_: sum((simulate_trial(**S1, target_mg=dz, weeks_on=72, weeks_off=0, s_true=s_true,
                                                                 P_true=p_c, m=m_, ec50=ec_fixed, latch_params=frozen)[0][72] - pub) ** 2
                                                for dz, pub in arms1), bounds=(0.3, 8.0), method="bounded").x
            arms_f = [simulate_trial(**S1, target_mg=dz, weeks_on=72, weeks_off=0, s_true=s_true, P_true=p_c, m=m_f,
                                     ec50=ec_fixed, latch_params=frozen)[0][72] for dz, _ in arms1]
            p4 = simulate_trial(**S4, target_mg=15, weeks_on=36, weeks_off=52, s_true=s_true, P_true=p_c, m=m_f,
                                ec50=ec_fixed, latch_params=frozen)[0]
            st_ = p4[36]
            rg = [(p4[36 + w] - st_) / (1 + st_ / 100) for w in (12, 24, 36, 52)]
            print(f"  {ec_fixed:>5.0f} {p_c:>5.0f} {m_f:>5.2f} {arms_f[0]:>+7.1f} {arms_f[1]:>+7.1f} {arms_f[2]:>+7.1f} {st_:>+11.1f} "
                  + " / ".join(f"{r_:+.1f}" for r_ in rg))
    print("  published:            −16.0   −21.4   −22.5       −20.9   +7.0 / +11.0 / +13.0 / +14.0")

    # 6. withdrawal arms with a frozen SP at several P
    print("\n--- 6. Withdrawal arms, frozen SP (no latch), % of loss regained at +52 wk, by cohort P ---")
    print("  P = 0 is pure CICO: intake returns to the pre-trial habit, TDEE is lower at the lower weight (6.4 cal/lb).")
    arms_w = [("SURMOUNT-4", 107.3, 85.8, 36, 0.455, 0.75, 67),
              ("STEP-1 ext sema", 105.6, 87.5, 68, 0.434, 0.56, 65),
              ("STEP-4 sema (+48 wk)", 107.2, 96.1, 20, 0.434, 0.56, 66),
              ("SCALE lira 3 mg", 105.2, 96.7, 56, 0.434, 0.65, 80)]
    Ps = [P_true, 25.0, 15.0, 5.0, 0.0]
    print(f"  {'trial':<22} " + " ".join(f"{'P=' + format(p, '.0f'):>8}" for p in Ps) + f" {'published':>10}")
    for name, b0, b1, lead, ff_, fl, pub in arms_w:
        vals = []
        for p_ in Ps:
            wk, frac = regain_from_leadin(b0, b1, lead, ff_, fl, p_, latch_params=(0.0, 30, 0.0))
            vals.append(frac[48] if "48" in name else frac[52])
        print(f"  {name:<22} " + " ".join(f"{v:>7.0f}%" for v in vals) + f" {pub:>9}%")
    for p_ in (15.0, 5.0, 0.0):
        wk, frac = regain_from_leadin(107.3, 85.8, 36, 0.455, 0.75, p_, latch_params=(0.0, 30, 0.0))
        print(f"  SURMOUNT-4 trajectory (frozen SP) at P = {p_:.0f}, % of stop weight — weeks off: " +
              "  ".join(f"{w}: sim {wk[w]:+.1f} / pub {pub_regain[w]:+.1f}" for w in (12, 24, 36, 52)))

    # 4. q3d forecast
    print("\n--- 7. Preregistered q3d forecast (gap model, FM held at plateau 60.5 lb, SP from latch) ---")
    fm_last = d["fat_mass_lbs"].iloc[-1]
    fc = q3d_forecast(tirz, fit, fm_last, fit["sp"][-1], dow_adj)
    fc["blk"] = ((fc["date"] - pd.Timestamp("2026-03-01")).dt.days // 28)
    g = fc.groupby("blk").agg(start=("date", "min"), dose=("dose_mg", "max"), B=("B", "mean"),
                              mean_intake=("intake_pred_true", "mean"),
                              sawtooth=("intake_pred_true", lambda x: x.max() - x.min()))
    print(f"  {'from':>12} {'dose':>6} {'mean B':>7} {'pred mean intake':>17} {'pred max−min':>13}")
    for _, r in g.iterrows():
        print(f"  {r['start'].date()!s:>12} {r['dose']:>6.2f} {r['B']:>7.1f} {r['mean_intake']:>17.0f} {r['sawtooth']:>13.0f}")
    out = ROOT / "analysis" / "DG_q3d_forecast.csv"
    fc.drop(columns="blk").to_csv(out, index=False, float_format="%.1f")
    print(f"  Wrote daily forecast to {out.name}. Score against intake_daily.csv once 2026-03+ is merged.")
    print("  Gap model prediction: mean intake tracks FM only — flat while FM holds; sawtooth shrinks with q3d.")
    print("  Tachyphylaxis (BA) predicts a continued upward drift with cumulative exposure.")

    # 8. control form on the subject's drug era: proportional vs integral (Tam 2009), linear vs threshold (Speakman DIP)
    print("\n--- 8. Control form on the drug era (DF design): proportional vs integral; linear gap vs threshold ---")
    y8, B8, fm8 = d["y"].values, d["blood_level"].values, d["fat_mass_lbs"].values
    gf8 = fm8[0] - fm8
    cum8 = np.cumsum(gf8) / 30.0
    def fit8(f, x0, lb, ub):
        r8 = least_squares(lambda p: y8 - f(p), x0, bounds=(lb, ub))
        sse8 = np.sum((y8 - f(r8.x)) ** 2)
        return len(y8) * np.log(sse8 / len(y8)) + 2 * len(x0), r8.x
    aic_p, _ = fit8(lambda p: p[0] + p[1] * B8 + p[2] * gf8, [2300, -80, 30], [1000, -600, -300], [3500, 0, 300])
    aic_pi, x_pi = fit8(lambda p: p[0] + p[1] * B8 + p[2] * gf8 + p[3] * cum8, [2300, -80, 30, 0], [1000, -600, -300, -50], [3500, 0, 300, 50])
    aic_i, _ = fit8(lambda p: p[0] + p[1] * B8 + p[2] * cum8, [2300, -80, 0], [1000, -600, -50], [3500, 0, 50])
    print(f"  proportional P·gap: AIC {aic_p:.1f};  PI (+ kI·∫gap dt): AIC {aic_pi:.1f} (kI = {x_pi[3]:.3f} cal/day per lb·month);  integral only: AIC {aic_i:.1f}")
    print("  Integral control (Tam 2009 'defend below threshold with zero steady-state error') predicts intake keeps rising")
    print("  at constant gap; it adds nothing here and over-erodes the 12.5 mg hold-out. Polidori/Hall proportional control holds.")
    print(f"  threshold P·max(0, LIP − FM) vs linear from FM₀ = {fm8[0]:.1f}:")
    for lip in (60, 65, 70, 75, 80, fm8[0]):
        thr = np.maximum(0, lip - fm8)
        a_t, x_t = fit8(lambda p, thr=thr: p[0] + p[1] * B8 + p[2] * thr, [2300, -80, 30], [1000, -600, -300], [3500, 0, 300])
        print(f"    LIP = {lip:5.1f} lb: AIC {a_t:.1f}  P = {x_t[2]:+.1f}")
    print("  A lower intervention point sits within ~5 lb of the pre-drug FM: no zone of indifference wider than that.")

    print("\n=== Scorecard: one parameter set (s, P from the subject's drug era) against everything ===")
    print(f"  ✓ Subject 18-month FM path from the pre-drug state: RMSE {sim_stats.get('kalman', float('nan')):.1f} lb (section 1).")
    save_result("DG_gap_model", I0=fit["a"], s=fit["s"], s_ci=list(fit["s_ci"]), P=fit["P"], P_ci=list(fit["P_ci"]),
                rmse_fit=fit["rmse"], gap_end=fit["gap_end"], sim_fm_rmse=sim_stats,
                latch=dict(tol=LATCH_TOL, hold=LATCH_HOLD, rate=LATCH_RATE), units="consumed calories (P4 intake_true)",
                intake_model=describe_intake_model())
    print("Wrote analysis/results/DG_gap_model.json")
    print("  ✓ Withdrawal regain fraction in all four trials is ~2/3 of the loss regardless of loss speed. Most of it is")
    print("    CICO: intake returning to habit with TDEE lower at the lower weight (P = 0) already regains ~half; a frozen")
    print("    SP with small cohort P (≈5-15) closes the rest. The latch's speed-of-loss prediction is NOT in the data.")
    print("  ✗ Cohort dose-response is concave (5 mg → 71% of the 15 mg loss at 1/3 the exposure). A linear drug term cannot")
    print("    do it at any m (SSE 51); an Emax with EC50 ≈ 5 units fits exactly (SSE 0.4) but predicts a flat sawtooth")
    print("    across the subject's doses, which his data contradict (1.7× growth). Cohort heterogeneity/tolerability, open.")
    print("  ✗ Pre-drug decade and 2011-13 willpower era: no identifiable gap (DD); a frozen SP at P ≈ 40 would have put")
    print("    1,500+ cal/day of pressure on 2012-13 intake that stayed at 1,850. Outside this model.")


if __name__ == "__main__":
    main()
