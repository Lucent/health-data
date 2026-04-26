"""CB. What strength workouts actually do, using the Kalman's daily LM state.

Two questions:

  1. Do workouts add lean mass per session, and at what rate?
  2. Do they preserve lean during cuts at a different ratio than cutting without workouts?

AA fit +41 g per session with a 275-day half-life from the 70 composition scans.
With daily LM available as a state, we can test these questions much more directly.

Method (Q1):
  Regress dLM/dt on (today's workout indicator, deficit, FM) over the full series.
  The coefficient on the workout indicator is the per-session lean response.

Method (Q2):
  Filter to days inside loss runs (rolling-7d intake < TDEE − 200 cal).
  Compute "lean-fraction-of-loss" = -dLM / -dW per loss run.
  Compare runs with high vs zero workout density.
  Forbes baseline says (1 - ff(FM)) of weight loss is lean. Compare observed against Forbes.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
LBS_TO_KG = 0.453592
FORBES_C_KG = 10.4


def forbes_lean_frac(fm_lbs):
    fm_kg = max(fm_lbs, 5.0) * LBS_TO_KG
    ff = 1.0 - (FORBES_C_KG ** 2) / (FORBES_C_KG + fm_kg) ** 2
    return 1.0 - ff


def load():
    kf = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    strength = pd.read_csv(ROOT / "workout" / "strength.csv", parse_dates=["date"])
    return kf, intake, strength


def main():
    kf, intake, strength = load()
    df = kf.rename(columns={"fat_mass_lbs": "fm_lbs", "tdee": "tdee_cal"}).merge(
        intake[["date", "calories"]], on="date", how="left"
    )
    df = df.dropna(subset=["fm_lbs", "lm_lbs", "tdee_cal"]).sort_values("date").reset_index(drop=True)

    # Workout indicators
    workout_dates = set(pd.to_datetime(strength["date"]).dt.normalize())
    df["workout"] = df["date"].dt.normalize().isin(workout_dates).astype(int)
    df["workout_7d"] = df["workout"].rolling(7, min_periods=1).sum()
    df["workout_30d"] = df["workout"].rolling(30, min_periods=1).sum()

    # Daily deltas
    df["dLM"] = df["lm_lbs"].diff()
    df["dFM"] = df["fm_lbs"].diff()
    df["dW"] = df["dLM"] + df["dFM"]
    df["deficit"] = df["calories"] - df["tdee_cal"]
    df["forbes_lean_frac"] = df["fm_lbs"].apply(forbes_lean_frac)

    print(f"n = {len(df)} days  workout days = {df.workout.sum()}")
    print(f"daily LM std = {df.dLM.std():.4f} lb/day")

    # =================================================================
    # Q1: Per-session lean response
    # =================================================================
    print("\n" + "=" * 70)
    print("Q1: Per-session lean increase")
    print("=" * 70)

    print("\n--- Same-day workout regression: dLM ~ workout + deficit + FM ---")
    sub = df.dropna(subset=["dLM", "deficit", "fm_lbs", "workout"]).copy()
    X = np.column_stack([sub.workout, sub.deficit / 1000.0, sub.fm_lbs, np.ones(len(sub))])
    y = sub.dLM.values
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    r = np.corrcoef(y, yhat)[0, 1]
    print(f"  workout coef:    {beta[0]:+.4f} lb/session  ({beta[0]*453.6:+.0f} g/session)")
    print(f"  deficit coef:    {beta[1]:+.4f} lb per 1000 cal/day deficit")
    print(f"  FM coef:         {beta[2]:+.5f} lb LM per lb FM")
    print(f"  intercept:       {beta[3]:+.4f}")
    print(f"  fit r:           {r:.3f}")

    # Bootstrap SE
    rng = np.random.default_rng(7)
    boot_ws = []
    n = len(sub)
    for _ in range(500):
        idx = rng.integers(0, n, n)
        Xi = X[idx]; yi = y[idx]
        b, *_ = np.linalg.lstsq(Xi, yi, rcond=None)
        boot_ws.append(b[0])
    print(f"  workout coef bootstrap: median {np.median(boot_ws)*453.6:+.0f} g, "
          f"95% CI [{np.percentile(boot_ws, 2.5)*453.6:+.0f}, {np.percentile(boot_ws, 97.5)*453.6:+.0f}] g")

    print("\n--- Stratified by workout-day vs not (matched FM bands) ---")
    print(f"  {'FM band':<14} {'n_w':>5} {'dLM_w':>10} {'n_no':>5} {'dLM_no':>10} {'Δ (g)':>10}")
    for lo, hi in [(20, 40), (40, 60), (60, 80), (80, 110)]:
        m = (sub.fm_lbs >= lo) & (sub.fm_lbs < hi)
        ww = sub[m & (sub.workout == 1)]
        nn = sub[m & (sub.workout == 0)]
        dw = ww.dLM.mean() if len(ww) else float("nan")
        dn = nn.dLM.mean() if len(nn) else float("nan")
        delta = (dw - dn) * 453.6
        print(f"  FM {lo}-{hi}      {len(ww):>5} {dw:>10.4f} {len(nn):>5} {dn:>10.4f} {delta:>+10.0f}")

    print("\n--- Lag profile: dLM vs workout at lags 0..14 ---")
    for lag in [0, 1, 3, 7, 14]:
        s = sub.copy()
        s["w_lag"] = s.workout.shift(lag).fillna(0)
        b, *_ = np.linalg.lstsq(np.column_stack([s.w_lag, s.deficit/1000, s.fm_lbs, np.ones(len(s))]),
                                 s.dLM.values, rcond=None)
        print(f"  lag {lag:>2}d: workout coef = {b[0]*453.6:+.0f} g/session")

    # =================================================================
    # Q2: Lean-fraction-of-loss in cuts with vs without workouts
    # =================================================================
    print("\n" + "=" * 70)
    print("Q2: Lean preservation during cuts")
    print("=" * 70)

    # Identify loss runs (rolling 7-day deficit < -200 cal/day; intake below TDEE)
    df["deficit_7d"] = df.deficit.rolling(7, min_periods=4).mean()
    df["in_cut"] = (df.deficit_7d < -200).astype(int)

    # Find contiguous cut runs of length >= 14 days
    runs = []
    i = 0
    while i < len(df):
        if df.in_cut.iloc[i] == 1:
            j = i
            while j < len(df) and df.in_cut.iloc[j] == 1:
                j += 1
            if j - i >= 14:
                runs.append((i, j))
            i = j
        else:
            i += 1

    print(f"\nLoss runs (≥14 d, rolling 7d deficit < -200 cal): n = {len(runs)}")
    print(f"  {'start':<12} {'end':<12} {'days':>4} {'workouts':>8} "
          f"{'dW':>7} {'dFM':>7} {'dLM':>7} {'lean%':>7} {'forbes%':>8} {'Δfrom':>7}")

    rows = []
    for s, e in runs:
        days = e - s
        sub = df.iloc[s:e]
        dW = sub.dW.sum()
        dFM = sub.dFM.sum()
        dLM = sub.dLM.sum()
        n_workouts = sub.workout.sum()
        if dW < -1.0:  # actually lost weight
            lean_frac = dLM / dW * 100
            mean_fm = sub.fm_lbs.mean()
            forbes_frac = forbes_lean_frac(mean_fm) * 100
            delta = lean_frac - forbes_frac
            rows.append({
                "start": str(sub.date.iloc[0])[:10],
                "end": str(sub.date.iloc[-1])[:10],
                "days": days,
                "workouts": n_workouts,
                "workout_per_wk": n_workouts / days * 7,
                "dW": dW,
                "dFM": dFM,
                "dLM": dLM,
                "mean_fm": mean_fm,
                "lean_frac_pct": lean_frac,
                "forbes_pct": forbes_frac,
                "delta_pct": delta,
            })

    runs_df = pd.DataFrame(rows)
    if len(runs_df):
        for _, r in runs_df.iterrows():
            print(f"  {r.start:<12} {r.end:<12} {r.days:>4} {r.workouts:>8.0f} "
                  f"{r.dW:>+7.1f} {r.dFM:>+7.1f} {r.dLM:>+7.1f} "
                  f"{r.lean_frac_pct:>+6.1f}% {r.forbes_pct:>+7.1f}% {r.delta_pct:>+6.1f}%")

        # Split by workout density
        runs_df["workout_density"] = pd.cut(
            runs_df.workout_per_wk, bins=[-0.01, 0.01, 1.0, 99],
            labels=["none", "low (<1/wk)", "high (≥1/wk)"],
        )
        print(f"\n--- By workout density ---")
        print(f"  {'density':<15} {'n_runs':>6} {'mean_lean%':>11} {'forbes%':>8} "
              f"{'Δlean%':>8} {'mean dW':>9}")
        for cat in ["none", "low (<1/wk)", "high (≥1/wk)"]:
            g = runs_df[runs_df.workout_density == cat]
            if len(g):
                print(f"  {cat:<15} {len(g):>6} {g.lean_frac_pct.mean():>+10.1f}% "
                      f"{g.forbes_pct.mean():>+7.1f}% {g.delta_pct.mean():>+7.1f}% "
                      f"{g.dW.mean():>+9.1f}")

        # Regression: delta-from-forbes ~ workouts/wk + total dW + mean FM
        if len(runs_df) >= 8:
            X = np.column_stack([runs_df.workout_per_wk, runs_df.dW, runs_df.mean_fm,
                                 np.ones(len(runs_df))])
            y = runs_df.delta_pct.values
            b, *_ = np.linalg.lstsq(X, y, rcond=None)
            print(f"\n--- Regression: lean_frac_minus_forbes ~ workouts_per_wk + dW + FM ---")
            print(f"  workouts/wk:  {b[0]:+.2f} percentage points (positive = workouts preserve lean vs Forbes)")
            print(f"  total dW:     {b[1]:+.3f}")
            print(f"  mean FM:      {b[2]:+.3f}")
            print(f"  intercept:    {b[3]:+.2f}")

            # Bootstrap SE on workouts coef
            rng = np.random.default_rng(11)
            boot_w = []
            n_runs = len(runs_df)
            for _ in range(2000):
                idx = rng.integers(0, n_runs, n_runs)
                Xi = X[idx]; yi = y[idx]
                bi, *_ = np.linalg.lstsq(Xi, yi, rcond=None)
                boot_w.append(bi[0])
            print(f"  workouts/wk bootstrap 95% CI: "
                  f"[{np.percentile(boot_w, 2.5):+.2f}, {np.percentile(boot_w, 97.5):+.2f}]")

    runs_df.to_csv(ROOT / "analysis" / "CB_workouts_lean_runs.csv", index=False)
    print(f"\nWrote run-level data to analysis/CB_workouts_lean_runs.csv")


if __name__ == "__main__":
    main()
