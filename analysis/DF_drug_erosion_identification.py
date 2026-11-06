#!/usr/bin/env python3
"""DF. What erodes tirzepatide's intake suppression: receptor tachyphylaxis,
set point pressure, or neither?

AX/BZ fit the drug era as  intake = base + s·B(t)·exp(−K·AUC(t))  and read the
decay as receptor tachyphylaxis. But every term that rises monotonically through
the drug era — cumulative AUC, calendar time, the gap between a frozen or
lagging set point and falling fat mass — produces the same within-sample fit.
The erosion is not identified by the fit; it is identified by what happens
when the candidate drivers DIVERGE:

  * At constant dose, AUC and time keep growing even when FM plateaus. A set
    point gap stops growing when FM stops falling.
  * Dose steps raise B abruptly. Receptor saturation (an Emax curve) would make
    each step buy less than the last; tachyphylaxis says the per-unit effect is
    the same at every dose but decays with exposure.

Design
------
Exposure B(t) is the daily-mean concentration over the eating window (drugs/
extract.py, timestamp-aware). Intake is weekday-adjusted using the pre-drug
2022+ weekday profile (injections are always Tuesdays, so day-of-cycle is
collinear with day-of-week). Candidate models for weekday-adjusted intake y:

  L0   y = a + s·B                                   (no erosion)
  L1   y = a + s·B·exp(−K·AUC)                       (AX/BZ tachyphylaxis)
  L2   y = a + s·B + c·t                              (linear time drift)
  L3   y = a + s·B + P·(SP_ema − FM)                  (lagging EMA set point, HL 45 d)
  L4   y = a + s·B + P·(FM_0 − FM)                    (frozen set point = latch open)
  L5   y = a + s·B/(1 + B/B50)                         (Emax saturation, no decay)
  L6   y = a + s·B + P·(FM_0 − FM) + K-form           (frozen SP + tachyphylaxis)

Each is fit on the first 48 weeks (through the 10 mg era) and scored on the
12.5 mg era (28 weeks), and also fit on the full record for AIC. The FM
plateau at ~61 lb from late 2025 is the discriminating window: AUC and t keep
rising there, FM_0 − FM does not.

Also reported: within constant-dose eras, first-half vs second-half mean
intake and sawtooth amplitude. Tachyphylaxis predicts both the mean and the
peak-day suppression erode; an additive pressure shifts the mean only.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent.parent
TRAIN_END = "2025-08-19"     # start of 12.5 mg era
WEEKLY_END = "2026-03-31"    # q3d dosing starts 2026-04-01 (no intake data yet)
SP_HL = 45.0


def ema(x, hl):
    a = 1 - np.exp(-np.log(2) / hl)
    s = np.empty(len(x))
    s[0] = x[0]
    for j in range(1, len(x)):
        s[j] = s[j - 1] + a * (x[j] - s[j - 1])
    return s


def smooth_latch(fm, tol, hold, rate):
    """SmoothLatch SP: adapts toward the window mean only when the trailing
    `hold` days all sit within ±tol of that mean (CF's smoothlatch_step)."""
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        sp[i] = sp[i - 1]
        if i + 1 >= hold:
            window = fm[i + 1 - hold:i + 1]
            center = window.mean()
            if np.max(np.abs(window - center)) <= tol:
                sp[i] = sp[i - 1] + rate * (center - sp[i - 1])
    return sp


def load():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    # Pre-drug weekday profile (2022+)
    pre = intake[(intake["date"] >= "2022-01-01") & (intake["date"] < "2024-09-17")].copy()
    pre["dow"] = pre["date"].dt.day_of_week
    dow_adj = pre.groupby("dow")["calories"].mean() - pre["calories"].mean()

    # Set point from pre-drug FM, propagated with its own dynamics (as AX does)
    fm_all = kalman["fat_mass_lbs"].ffill().bfill().values
    sp_all = ema(fm_all, SP_HL)
    kalman = kalman.assign(sp_ema=sp_all)

    d = intake[["date", "calories"]].merge(
        tirz[["date", "blood_level", "days_since_injection", "dose_mg"]], on="date")
    d = d.merge(kalman[["date", "fat_mass_lbs", "tdee", "sp_ema"]], on="date")
    d = d[(d["date"] < WEEKLY_END) & (d["days_since_injection"] <= 6)]
    d = d.sort_values("date").reset_index(drop=True)
    d["y"] = d["calories"] - d["date"].dt.day_of_week.map(dow_adj)
    d["t_wk"] = (d["date"] - d["date"].min()).dt.days / 7.0
    d["auc"] = d["blood_level"].cumsum() / 7.0
    d["gap_ema"] = d["sp_ema"] - d["fat_mass_lbs"]
    d["gap_frozen"] = d["fat_mass_lbs"].iloc[0] - d["fat_mass_lbs"]
    # SmoothLatch set points run on the drug-era FM path from the pre-drug FM:
    # subject-fit (tol 3, hold 14, rate 0.01; BD/CF) and trial-fit (tol 1, hold 30, rate 0.0025; CF/BG)
    fm_path = d["fat_mass_lbs"].values
    d["gap_latch_subj"] = smooth_latch(fm_path, 3.0, 14, 0.01) - fm_path
    d["gap_latch_trial"] = smooth_latch(fm_path, 1.0, 30, 0.0025) - fm_path
    d["wk"] = (d["date"] - d["date"].min()).dt.days // 7
    return d, dow_adj


def fit_models(d):
    y = d["y"].values
    B = d["blood_level"].values
    t = d["t_wk"].values
    auc = d["auc"].values
    g_ema = d["gap_ema"].values
    g_frz = d["gap_frozen"].values
    g_ls = d["gap_latch_subj"].values
    g_lt = d["gap_latch_trial"].values
    print(f"  SP − FM at end of record: EMA {g_ema[-1]:+.1f}, frozen {g_frz[-1]:+.1f}, "
          f"SmoothLatch subject-fit {g_ls[-1]:+.1f}, SmoothLatch trial-fit {g_lt[-1]:+.1f} lb")

    models = {
        "L0 linear, no erosion":        (lambda p: p[0] + p[1] * B,
                                         [2300, -50], [1000, -600], [3500, 0]),
        "L1 linear × exp(−K·AUC)":      (lambda p: p[0] + p[1] * B * np.exp(-p[2] * auc),
                                         [2300, -90, 0.002], [1000, -600, 0], [3500, 0, 0.05]),
        "L2 linear + c·t":              (lambda p: p[0] + p[1] * B + p[2] * t,
                                         [2300, -80, 5], [1000, -600, -100], [3500, 0, 100]),
        "L3 linear + P·(SP_ema − FM)":  (lambda p: p[0] + p[1] * B + p[2] * g_ema,
                                         [2300, -80, 30], [1000, -600, -300], [3500, 0, 300]),
        "L4 linear + P·(FM0 − FM)":     (lambda p: p[0] + p[1] * B + p[2] * g_frz,
                                         [2300, -80, 20], [1000, -600, -300], [3500, 0, 300]),
        "L5 Emax s·B/(1+B/B50)":        (lambda p: p[0] + p[1] * B / (1 + B / p[2]),
                                         [2300, -100, 10], [1000, -600, 0.5], [3500, 0, 500]),
        "L6 frozen SP + exp(−K·AUC)":   (lambda p: p[0] + p[1] * B * np.exp(-p[3] * auc) + p[2] * g_frz,
                                         [2300, -90, 20, 0.001], [1000, -600, -300, 0], [3500, 0, 300, 0.05]),
        "L7 linear + P·(Latch_subj − FM)": (lambda p: p[0] + p[1] * B + p[2] * g_ls,
                                            [2300, -80, 30], [1000, -600, -300], [3500, 0, 300]),
        "L8 linear + P·(Latch_trial − FM)": (lambda p: p[0] + p[1] * B + p[2] * g_lt,
                                             [2300, -80, 30], [1000, -600, -300], [3500, 0, 300]),
    }
    train = np.where(d["date"] < TRAIN_END)[0]
    test = np.where(d["date"] >= TRAIN_END)[0]
    allidx = np.arange(len(y))

    def fit(f, x0, lb, ub, idx):
        r = least_squares(lambda p: y[idx] - f(p)[idx], x0, bounds=(lb, ub))
        pred = f(r.x)
        sse = float(np.sum((y[idx] - pred[idx]) ** 2))
        aic = len(idx) * np.log(sse / len(idx)) + 2 * len(x0)
        return r.x, pred, np.sqrt(sse / len(idx)), aic

    print("=" * 96)
    print("A. Full-record fits (weekday-adjusted intake, daily-mean exposure, n=%d)" % len(y))
    print("=" * 96)
    print(f"  {"model":34s} {"params":48s} {'RMSE':>6} {'AIC':>8}")
    full = {}
    for name, (f, x0, lb, ub) in models.items():
        p, pred, rmse, aic = fit(f, x0, lb, ub, allidx)
        full[name] = (p, aic)
        print(f"  {name:34s} {np.array2string(np.round(p, 4), max_line_width=200):48s} {rmse:6.0f} {aic:8.1f}")
    best = min(full.items(), key=lambda kv: kv[1][1])
    print(f"\n  Lowest AIC: {best[0]}   (ΔAIC vs L1 tachyphylaxis = {best[1][1] - full['L1 linear × exp(−K·AUC)'][1]:+.1f})")
    import sys as _sys; _sys.path.insert(0, str(ROOT / "analysis"))
    from results import save as _save
    _save("DF_drug_erosion_identification", best_model=best[0],
          d_aic_vs_tachyphylaxis=float(best[1][1] - full['L1 linear × exp(−K·AUC)'][1]),
          full_record={k: dict(params=[float(x) for x in np.atleast_1d(v[0])], aic=float(v[1])) for k, v in full.items()})
    print("  Wrote analysis/results/DF_drug_erosion_identification.json")

    print("\n" + "=" * 96)
    print(f"B. Fit through the 10 mg era (n={len(train)}), score on the 12.5 mg era (n={len(test)})")
    print("=" * 96)
    print(f"  {"model":34s} {"params (train)":48s} {'RMSE_in':>7} {'RMSE_out':>8} {'bias_out':>9}")
    oos = {}
    for name, (f, x0, lb, ub) in models.items():
        p, pred, rmse, aic = fit(f, x0, lb, ub, train)
        rmse_out = float(np.sqrt(np.mean((y[test] - pred[test]) ** 2)))
        bias = float(np.mean(y[test] - pred[test]))
        oos[name] = (rmse_out, bias)
        print(f"  {name:34s} {np.array2string(np.round(p, 4), max_line_width=200):48s} {rmse:7.0f} {rmse_out:8.0f} {bias:+9.0f}")
    print("\n  bias_out > 0: model predicts LESS intake than observed (under-erodes);")
    print("  bias_out < 0: model predicts MORE intake than observed (over-erodes).")
    return full, oos


def constant_dose_eras(d):
    print("\n" + "=" * 96)
    print("C. Constant-dose eras: first half vs second half (weekday-adjusted intake)")
    print("=" * 96)
    g = d.groupby("wk")
    w = pd.DataFrame({
        "date": g["date"].min(), "dose": g["dose_mg"].max(), "mean": g["y"].mean(),
        "fm": g["fat_mass_lbs"].mean(), "B": g["blood_level"].mean(),
        "peak": d[d["days_since_injection"].isin([1, 2])].groupby("wk")["y"].mean(),
        "trough": d[d["days_since_injection"].isin([5, 6])].groupby("wk")["y"].mean(),
    })
    w["amp"] = w["trough"] - w["peak"]
    print(f"  {'dose':>5} {'weeks':>5} | {'mean1':>5} {'mean2':>5} {'Δmean':>6} | {'peak1':>5} {'peak2':>5} {'Δpeak':>6} | "
          f"{'amp1':>5} {'amp2':>5} {'Δamp':>5} | {'FM1':>5} {'FM2':>5}")
    for dose, lo, hi in [(5.0, "2024-10-15", "2025-03-04"), (7.5, "2025-03-04", "2025-05-27"),
                         (10.0, "2025-05-27", "2025-08-19"), (12.5, "2025-08-19", WEEKLY_END)]:
        s = w[(w["date"] >= lo) & (w["date"] < hi)]
        h = len(s) // 2
        a, b = s.iloc[:h], s.iloc[h:]
        print(f"  {dose:>5.1f} {len(s):>5d} | {a['mean'].mean():5.0f} {b['mean'].mean():5.0f} {b['mean'].mean() - a['mean'].mean():+6.0f} | "
              f"{a['peak'].mean():5.0f} {b['peak'].mean():5.0f} {b['peak'].mean() - a['peak'].mean():+6.0f} | "
              f"{a['amp'].mean():5.0f} {b['amp'].mean():5.0f} {b['amp'].mean() - a['amp'].mean():+5.0f} | "
              f"{a['fm'].mean():5.1f} {b['fm'].mean():5.1f}")
    print("  Tachyphylaxis: Δmean > 0 AND Δamp < 0 in every era (per-unit effect shrinks).")
    print("  Additive pressure: Δmean > 0 with Δamp ≈ 0.  Nothing: Δmean ≈ 0.")

    # Dose-era means
    print("\n  Dose-era means (exposure → intake):")
    e = w.groupby("dose").agg(weeks=("mean", "size"), B=("B", "mean"), intake=("mean", "mean"),
                              amp=("amp", "mean"), fm=("fm", "mean")).round(1)
    print(e.to_string())


def plateau_window(d):
    print("\n" + "=" * 96)
    print("D. FM plateau window (FM within ±1.5 lb of its final 12.5 mg value) — do AUC and intake keep rising?")
    print("=" * 96)
    last_fm = d["fat_mass_lbs"].iloc[-1]
    plateau = d[(d["fat_mass_lbs"] - last_fm).abs() <= 1.5]
    start = plateau["date"].min()
    p = d[d["date"] >= start]
    half = len(p) // 2
    a, b = p.iloc[:half], p.iloc[half:]
    print(f"  Plateau starts {start.date()} ({len(p)} days, FM {p['fat_mass_lbs'].iloc[0]:.1f} → {p['fat_mass_lbs'].iloc[-1]:.1f} lb)")
    print(f"  first half: intake {a['y'].mean():.0f}  AUC {a['auc'].mean():.0f}  B {a['blood_level'].mean():.1f}")
    print(f"  second half: intake {b['y'].mean():.0f}  AUC {b['auc'].mean():.0f}  B {b['blood_level'].mean():.1f}")
    print(f"  Δintake = {b['y'].mean() - a['y'].mean():+.0f} cal/day while AUC rose {b['auc'].mean() - a['auc'].mean():.0f} unit-wk")


def main():
    d, dow_adj = load()
    print("DF. Drug-era erosion identification")
    print(f"  on-drug weekly-dosing days: {len(d)}  ({d['date'].min().date()} → {d['date'].max().date()})")
    print(f"  pre-drug weekday adjustment (cal): " +
          ", ".join(f"{k}:{v:+.0f}" for k, v in dow_adj.round(0).items()))
    full, oos = fit_models(d)
    constant_dose_eras(d)
    plateau_window(d)
    print("\n=== Reading ===")
    print("  Any monotone-in-time driver fits the erosion in-sample; the AIC differences between")
    print("  tachyphylaxis (L1), time (L2) and set point gap (L3/L4) are what the data can say.")
    print("  The 12.5 mg hold-out and the FM plateau are the only windows where the drivers diverge.")


if __name__ == "__main__":
    main()
