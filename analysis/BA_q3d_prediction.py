#!/usr/bin/env python3
"""BA. Preregistered prediction for q3d dosing switch (2026-04-01).

Model G: two-component drug effect on intake.
  calories = base + acute × (blood - EMA₁ᵈ(blood)) + tonic × blood × exp(-decay × weeks)

Fitted to 529 days of q7d data. Predicts forward under q3d (25 clicks,
5.2 mg every 3 days) and q7d counterfactual (12.5 mg weekly).

The multiplicative model (AX: -74.5 × blood × tachy) fails to predict
the 2026 sawtooth because tachyphylaxis crushes the within-week slope
to -19 cal/unit at 80 weeks, while observed is -70. Model G separates:
  - Acute response to blood level *changes* (does not decay)
  - Tonic suppression at steady-state blood level (decays, HL=55 weeks)
"""

from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy.optimize import minimize
import openpyxl

ROOT = Path(__file__).resolve().parent.parent

# ═══════════════════════════════════════════════════════════════
# PK parameters (FDA prescribing information)
# ═══════════════════════════════════════════════════════════════
HALF_LIFE = 5.0
KE = math.log(2) / HALF_LIFE
KA = 3.3122


def pk_curve(dose, t):
    if t <= 0:
        return 0.0
    return dose * KA / (KA - KE) * (math.exp(-KE * t) - math.exp(-KA * t))


def compute_ema(values, hl_days):
    alpha = 1 - math.exp(-math.log(2) / max(hl_days, 0.1))
    ema = np.empty(len(values))
    ema[0] = values[0]
    for i in range(1, len(values)):
        ema[i] = ema[i - 1] + alpha * (values[i] - ema[i - 1])
    return ema


def blood_on_day(day, injections):
    level = 0.0
    for inj_dt, inj_dose in injections:
        t = (day - inj_dt).total_seconds() / 86400
        if t < 0:
            break
        if t > 30:
            continue
        level += pk_curve(inj_dose, t)
    return level


def load_injections():
    """Load injection history from medicine.xlsx."""
    xlsx = ROOT / "drugs" / "medicine.xlsx"
    wb = openpyxl.load_workbook(xlsx, data_only=True)
    ws = wb.active
    injections = []
    for row in ws.iter_rows(min_row=2, values_only=True):
        date_val, drug, dose_val, clicks_val = row[0], row[1], row[2], row[3]
        if drug and ("tirzepatide" in str(drug).lower() or "zepbound" in str(drug).lower()):
            has_clicks = clicks_val is not None and str(clicks_val).strip() not in ("", "None")
            if has_clicks:
                dose_mg = int(float(clicks_val)) * float(dose_val) / 60
            else:
                dose_mg = float(dose_val)
            injections.append((pd.Timestamp(date_val), dose_mg))
    return injections


def main():
    # ═══════════════════════════════════════════════════════════
    # STEP 1: Fit Model G to historical q7d data
    # ═══════════════════════════════════════════════════════════
    print("=" * 70)
    print("1. FIT: Model G on historical data (raw calories ~ blood level)")
    print("=" * 70)

    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    df = intake[["date", "calories"]].merge(
        tirz[["date", "blood_level", "days_since_injection"]], on="date")
    df = df[df["blood_level"] > 0].sort_values("date").reset_index(drop=True)
    drug_start = df["date"].min()
    df["cum_days"] = (df["date"] - drug_start).dt.days

    blood = df["blood_level"].values
    cal = df["calories"].values
    cum_days = df["cum_days"].values
    n = len(df)

    # Grid search adapt_hl, optimize remaining params
    ADAPT_HL = 1  # days (best from grid search)
    ema = compute_ema(blood, ADAPT_HL)
    dev = blood - ema

    def objective(params):
        base, slope, tonic, decay = params
        if decay < 0:
            return 1e12
        pred = base + slope * dev + tonic * blood * np.exp(-decay * cum_days / 7)
        return np.mean((cal - pred) ** 2)

    res = minimize(objective, [2200, -100, -50, 0.013], method="Nelder-Mead",
                   options={"maxiter": 20000, "xatol": 0.01, "fatol": 0.1})
    BASE, ACUTE, TONIC, DECAY = res.x
    TACHY_HL = math.log(2) / DECAY if DECAY > 0 else float("inf")

    pred_fit = BASE + ACUTE * dev + TONIC * blood * np.exp(-DECAY * cum_days / 7)
    rmse = np.sqrt(np.mean((cal - pred_fit) ** 2))
    r = np.corrcoef(cal, pred_fit)[0, 1]

    print(f"\n  Model G parameters (n={n} days):")
    print(f"    base     = {BASE:.0f} cal (intake without drug)")
    print(f"    acute    = {ACUTE:+.1f} cal/unit (response to blood − EMA₁ᵈ)")
    print(f"    tonic    = {TONIC:+.1f} cal/unit blood (steady-state suppression)")
    print(f"    tachy HL = {TACHY_HL:.0f} weeks (tonic decay)")
    print(f"    RMSE     = {rmse:.0f} cal, r = {r:.3f}")

    # Validate sawtooth by era
    df["pred"] = pred_fit
    df["cum_weeks"] = df["cum_days"] / 7
    print(f"\n  Sawtooth validation (day 0-1 vs day 4-5):")
    for label, w0, w1 in [("Wk 0-20", 0, 20), ("Wk 20-40", 20, 40),
                           ("Wk 40-60", 40, 60), ("Wk 60-80", 60, 100)]:
        sub = df[(df["cum_weeks"] >= w0) & (df["cum_weeks"] < w1)]
        if len(sub) < 30:
            continue
        obs_sw = (sub[sub["days_since_injection"].isin([4, 5])]["calories"].mean()
                  - sub[sub["days_since_injection"].isin([0, 1])]["calories"].mean())
        pred_sw = (sub[sub["days_since_injection"].isin([4, 5])]["pred"].mean()
                   - sub[sub["days_since_injection"].isin([0, 1])]["pred"].mean())
        obs_mean = sub["calories"].mean()
        pred_mean = sub["pred"].mean()
        print(f"    {label}: sawtooth obs={obs_sw:+.0f} pred={pred_sw:+.0f}  "
              f"mean obs={obs_mean:.0f} pred={pred_mean:.0f}")

    # ═══════════════════════════════════════════════════════════
    # STEP 2: Forward prediction — q3d vs q7d
    # ═══════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("2. FORWARD PREDICTION: q3d (5.2 mg) vs q7d (12.5 mg)")
    print("=" * 70)

    past_injections = load_injections()
    last_inj = past_injections[-1]

    sim_start = pd.Timestamp("2026-04-01")
    sim_end = pd.Timestamp("2026-05-04")
    days = pd.date_range(sim_start, sim_end)
    warmup = pd.date_range(sim_start - pd.Timedelta(days=7),
                           sim_start - pd.Timedelta(days=1))

    # Build future schedules
    q3d_dose = 25 * 12.5 / 60
    q3d_future = []
    nxt = last_inj[0] + pd.Timedelta(days=3)
    while nxt <= sim_end:
        q3d_future.append((nxt, q3d_dose))
        nxt += pd.Timedelta(days=3)

    q7d_future = []
    nxt = last_inj[0] + pd.Timedelta(days=7)
    while nxt <= sim_end:
        q7d_future.append((nxt, 12.5))
        nxt += pd.Timedelta(days=7)

    def simulate(injections, future_injs):
        all_days = list(warmup) + list(days)
        bloods = [blood_on_day(d, injections) for d in all_days]
        ema_vals = compute_ema(np.array(bloods), ADAPT_HL)

        warmup_n = len(warmup)
        results = []
        for i, day in enumerate(days):
            idx = warmup_n + i
            bl = bloods[idx]
            em = ema_vals[idx]
            deviation = bl - em
            weeks = (day - drug_start).days / 7.0
            tachy = math.exp(-DECAY * weeks)

            acute_eff = ACUTE * deviation
            tonic_eff = TONIC * bl * tachy
            predicted = BASE + acute_eff + tonic_eff

            is_inj = any(d.date() == day.date() for d, _ in future_injs)
            results.append(dict(
                date=day, blood=bl, ema=em, deviation=deviation,
                tachy=tachy, acute=acute_eff, tonic=tonic_eff,
                predicted=predicted, is_inj=is_inj,
            ))
        return pd.DataFrame(results)

    r_q3d = simulate(past_injections + q3d_future, q3d_future)
    r_q7d = simulate(past_injections + q7d_future, q7d_future)

    for label, rdf in [("q3d (25 clicks, 5.2 mg q3d)", r_q3d),
                        ("q7d (12.5 mg weekly)", r_q7d)]:
        print(f"\n  === {label} ===\n")
        print(f"  {'Date':>8} {'Blood':>6} {'EMA':>6} {'Dev':>6} {'Acute':>7} "
              f"{'Tonic':>7} {'Intake':>7}")
        print(f"  {'-' * 56}")
        for _, r in rdf.iterrows():
            mark = " ←inj" if r["is_inj"] else ""
            print(f"  {r['date'].strftime('%b %d'):>8} {r['blood']:6.1f} "
                  f"{r['ema']:6.1f} {r['deviation']:+6.1f} {r['acute']:>+7.0f} "
                  f"{r['tonic']:>+7.0f} {r['predicted']:>7.0f}{mark}")

    # ═══════════════════════════════════════════════════════════
    # STEP 3: Summary
    # ═══════════════════════════════════════════════════════════
    TDEE = 2063  # 2026 mean

    print(f"\n{'=' * 70}")
    print("3. PREREGISTERED PREDICTIONS (2026-04-01)")
    print("=" * 70)
    print(f"\n  {'':>24} {'q3d 5.2mg':>12} {'q7d 12.5mg':>12} {'Diff':>8}")
    print(f"  {'Mean intake':>24} {r_q3d['predicted'].mean():>12.0f} "
          f"{r_q7d['predicted'].mean():>12.0f} "
          f"{r_q3d['predicted'].mean() - r_q7d['predicted'].mean():>+8.0f}")
    print(f"  {'Intake std':>24} {r_q3d['predicted'].std():>12.0f} "
          f"{r_q7d['predicted'].std():>12.0f} "
          f"{r_q3d['predicted'].std() - r_q7d['predicted'].std():>+8.0f}")
    print(f"  {'Min (deepest day)':>24} {r_q3d['predicted'].min():>12.0f} "
          f"{r_q7d['predicted'].min():>12.0f} "
          f"{r_q3d['predicted'].min() - r_q7d['predicted'].min():>+8.0f}")
    print(f"  {'Max (weakest day)':>24} {r_q3d['predicted'].max():>12.0f} "
          f"{r_q7d['predicted'].max():>12.0f} "
          f"{r_q3d['predicted'].max() - r_q7d['predicted'].max():>+8.0f}")
    print(f"  {'Mean acute effect':>24} {r_q3d['acute'].mean():>+12.0f} "
          f"{r_q7d['acute'].mean():>+12.0f}")
    print(f"  {'Mean tonic effect':>24} {r_q3d['tonic'].mean():>+12.0f} "
          f"{r_q7d['tonic'].mean():>+12.0f}")
    print(f"\n  TDEE ≈ {TDEE} cal (2026 mean)")
    print(f"  q3d deficit: {r_q3d['predicted'].mean() - TDEE:+.0f} cal/day")
    print(f"  q7d deficit: {r_q7d['predicted'].mean() - TDEE:+.0f} cal/day")

    print(f"\n  Falsification criteria:")
    print(f"    Model is wrong if q3d mean intake deviates from "
          f"{r_q3d['predicted'].mean():.0f} by >100 cal")
    print(f"    Model is wrong if q3d sawtooth (3-day range) exceeds "
          f"{r_q3d['predicted'].max() - r_q3d['predicted'].min():.0f} cal")


if __name__ == "__main__":
    main()
