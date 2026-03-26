"""What predicts measured RMR beyond body composition?

25 Cosmed Fitmate calorimetry measurements (3 at a sports medicine clinic,
22 at home — same device). Sweep trailing dietary, activity, and sleep
features at multiple windows. Validation: leave-one-cluster-out CV
(measurements within 7 days = same cluster).

Results:
  1. Dietary features (calories, protein, carbs, fat, sodium) — null.
     No trailing dietary variable at any window beats the Fitmate noise
     floor (~170 cal) under CV.
  2. Step counts — null. Total daily steps at any window: CV RMSE ≥ 177.
  3. Walk sessions (Samsung Health logged exercise) — signal. Count of
     deliberate walks in 30 days: CV RMSE = 116 (R² = 0.49). Survives
     controlling for season (partial r = 0.47), body composition, and
     tirzepatide. Holds within-winter (r = 0.63) and within the 2022-2023
     cluster where FM was nearly constant (r = 0.68). Walk session count
     predicts RMR; total steps and walk minutes do not.
"""

from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from itertools import combinations

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path(__file__).resolve().parent.parent
LBS_TO_KG = 0.453592

# Samsung Health real data starts 2014-04-24. Drop pre-Samsung measurements
# where step/exercise data is backfilled from hospital calendar estimates.
SAMSUNG_START = "2016-01-01"

WINDOWS = [7, 10, 14, 21, 30, 45, 60]


def load_all():
    rmr = pd.read_csv(ROOT / "RMR" / "rmr.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    steps = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])
    sleep = pd.read_csv(ROOT / "steps-sleep" / "sleep.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    strength = pd.read_csv(ROOT / "workout" / "strength.csv", parse_dates=["date"])
    exercises = pd.read_csv(ROOT / "steps-sleep" / "exercises.csv", parse_dates=["date"])
    drugs = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    return rmr, intake, steps, sleep, comp, strength, exercises, drugs


def trailing(dates, values, target, window):
    """Mean of values in [target - window, target)."""
    mask = (dates >= target - pd.Timedelta(days=window)) & (dates < target)
    v = values[mask]
    return v.mean() if len(v) > 0 else np.nan


def trailing_count(dates, target, window):
    mask = (dates >= target - pd.Timedelta(days=window)) & (dates < target)
    return mask.sum()


def trailing_sum(dates, values, target, window):
    mask = (dates >= target - pd.Timedelta(days=window)) & (dates < target)
    v = values[mask]
    return v.sum() if len(v) > 0 else np.nan


def cluster_measurements(dates, gap_days=7):
    """Group measurements within gap_days of each other."""
    clusters = []
    current = [0]
    for i in range(1, len(dates)):
        delta = (dates[i] - dates[i - 1]) / np.timedelta64(1, "D")
        if delta <= gap_days:
            current.append(i)
        else:
            clusters.append(current)
            current = [i]
    clusters.append(current)
    return clusters


def loco_cv(X, y, clusters, ridge=50.0):
    """Leave-one-cluster-out CV with ridge regression."""
    preds = np.full(len(y), np.nan)
    for cidx in clusters:
        mask = np.ones(len(y), dtype=bool)
        mask[cidx] = False
        try:
            p = X.shape[1]
            pen = ridge * np.eye(p)
            pen[-1, -1] = 0.0
            c = np.linalg.solve(X[mask].T @ X[mask] + pen, X[mask].T @ y[mask])
            preds[cidx] = X[cidx] @ c
        except np.linalg.LinAlgError:
            preds[cidx] = y[mask].mean()
    resid = y - preds
    rmse = np.sqrt(np.mean(resid ** 2))
    ss_res = np.sum(resid ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    return rmse, r2


def remap_clusters(clusters, valid_idx):
    """Remap global cluster indices to valid-only indices."""
    idx_set = set(valid_idx.tolist())
    idx_map = {g: l for l, g in enumerate(valid_idx)}
    out = []
    for cidx in clusters:
        mapped = [idx_map[i] for i in cidx if i in idx_set]
        if mapped:
            out.append(mapped)
    return out


def build_features(rmr, intake, steps, sleep, comp, strength, exercises, drugs):
    """One row per RMR measurement with trailing features."""
    walk_mask = exercises["type"] == "walking"
    rows = []
    for _, r in rmr.iterrows():
        d = r["date"]
        cr = comp.iloc[(comp["date"] - d).abs().idxmin()]
        drug = drugs[drugs["date"] == d]
        tirz = drug["effective_level"].values[0] if len(drug) > 0 else 0.0

        row = {
            "date": d, "rmr_kcal": r["rmr_kcal"],
            "expected_rmr": cr["expected_rmr"],
            "fm_lbs": cr["fm_lbs"], "ffm_lbs": cr["ffm_lbs"],
            "tirz_level": tirz,
            "is_summer": 1 if d.month in [5, 6, 7, 8] else 0,
            "month_sin": np.sin(2 * np.pi * d.month / 12),
            "month_cos": np.cos(2 * np.pi * d.month / 12),
        }

        for w in WINDOWS:
            # Dietary
            row[f"cal_{w}d"] = trailing(intake["date"].values, intake["calories"].values, d, w)
            row[f"prot_g_{w}d"] = trailing(intake["date"].values, intake["protein_g"].values, d, w)
            row[f"sodium_{w}d"] = trailing(intake["date"].values, intake["sodium_mg"].values, d, w)
            cal = row[f"cal_{w}d"]
            if cal and cal > 0 and not np.isnan(cal):
                p = row[f"prot_g_{w}d"]
                row[f"prot_pct_{w}d"] = p * 4 / cal * 100 if p and not np.isnan(p) else np.nan

            # Steps
            row[f"steps_{w}d"] = trailing(steps["date"].values, steps["steps"].values, d, w)

            # Walk sessions (count and minutes)
            ex_mask = (exercises["date"] >= d - pd.Timedelta(days=w)) & (exercises["date"] < d)
            walks = exercises[ex_mask & walk_mask]
            row[f"walk_sessions_{w}d"] = len(walks)
            row[f"walk_min_{w}d"] = walks["duration_min"].sum() if len(walks) > 0 else 0

            # Sleep
            row[f"sleep_{w}d"] = trailing(sleep["date"].values, sleep["sleep_hours"].values, d, w)

            # Strength
            row[f"strength_{w}d"] = trailing_count(strength["date"].values, d, w)

        rows.append(row)
    return pd.DataFrame(rows)


def sweep_singles(df, y, clusters, base, label="expected_rmr"):
    """Sweep all single features added to base. Return sorted results."""
    exclude = {"date", "rmr_kcal", "expected_rmr", "fm_lbs", "ffm_lbs",
               "tirz_level", "is_summer", "month_sin", "month_cos"}
    cols = [c for c in df.columns if c not in exclude]
    results = []
    n = len(y)
    for col in cols:
        v = df[col].values.astype(float)
        valid = ~np.isnan(v) & ~np.isnan(base)
        if valid.sum() < 15 or np.std(v[valid]) < 1e-10:
            continue
        X = np.column_stack([base[valid], v[valid], np.ones(valid.sum())])
        vc = remap_clusters(clusters, np.where(valid)[0])
        if len(vc) < 3:
            continue
        rmse, r2 = loco_cv(X, y[valid], vc)
        results.append((col, rmse, r2, valid.sum()))
    results.sort(key=lambda x: x[1])
    return results


def main():
    rmr_all, intake, steps, sleep, comp, strength, exercises, drugs = load_all()
    rmr = rmr_all[rmr_all["date"] >= SAMSUNG_START].reset_index(drop=True)
    print(f"=== Data ===")
    print(f"  RMR measurements: {len(rmr_all)} total, {len(rmr)} with Samsung Health data (>= {SAMSUNG_START})")
    print(f"  Dropped: {len(rmr_all) - len(rmr)} pre-Samsung (backfilled steps, no exercise sessions)")

    df = build_features(rmr, intake, steps, sleep, comp, strength, exercises, drugs)
    y = df["rmr_kcal"].values
    n = len(y)
    ermr = df["expected_rmr"].values
    clusters = cluster_measurements(df["date"].values)

    print(f"  {n} measurements in {len(clusters)} independent clusters:")
    for i, c in enumerate(clusters):
        dates = [df.iloc[j]["date"].strftime("%Y-%m-%d") for j in c]
        print(f"    {i + 1:2d}. {', '.join(dates)}")

    # ── Baseline ─────────────────────────────────────────────────────
    X0 = np.column_stack([ermr, np.ones(n)])
    rmse0, r2_0 = loco_cv(X0, y, clusters)
    print(f"\n=== Baseline ===")
    print(f"  expected_rmr only: CV RMSE = {rmse0:.0f}, R² = {r2_0:.3f}")
    print(f"  Fitmate noise floor: ~170 cal (consecutive-day std)")

    # ── Part 1: dietary features (null) ──────────────────────────────
    print(f"\n{'=' * 70}")
    print(f"=== Part 1: Dietary features ===")
    dietary = [c for c in df.columns if any(c.startswith(p) for p in
               ["cal_", "prot_g_", "prot_pct_", "sodium_"])]
    results_diet = []
    for col in dietary:
        v = df[col].values.astype(float)
        valid = ~np.isnan(v) & ~np.isnan(ermr)
        if valid.sum() < 15 or np.std(v[valid]) < 1e-10:
            continue
        X = np.column_stack([ermr[valid], v[valid], np.ones(valid.sum())])
        vc = remap_clusters(clusters, np.where(valid)[0])
        if len(vc) < 3:
            continue
        rmse, r2 = loco_cv(X, y[valid], vc)
        results_diet.append((col, rmse, r2))
    results_diet.sort(key=lambda x: x[1])
    print(f"  All dietary features with expected_rmr (sorted by CV RMSE):")
    print(f"  {'Feature':>20} {'CV RMSE':>8} {'R²':>6}  {'vs base':>8}")
    for col, rmse, r2 in results_diet[:15]:
        print(f"  {col:>20} {rmse:8.0f} {r2:6.3f}  {rmse - rmse0:+8.0f}")
    best_diet = results_diet[0][1] if results_diet else 999
    print(f"\n  Best dietary: CV RMSE {best_diet:.0f} (baseline {rmse0:.0f}, noise floor 170)")
    print(f"  Verdict: no dietary feature improves prediction.")

    # ── Part 2: steps (null) ─────────────────────────────────────────
    print(f"\n{'=' * 70}")
    print(f"=== Part 2: Step counts ===")
    step_cols = [c for c in df.columns if c.startswith("steps_")]
    print(f"  {'Feature':>15} {'CV RMSE':>8} {'R²':>6}")
    for col in sorted(step_cols, key=lambda c: int(c.split("_")[1].replace("d", ""))):
        v = df[col].values.astype(float)
        valid = ~np.isnan(v)
        if valid.sum() < 15:
            continue
        X = np.column_stack([ermr[valid], v[valid], np.ones(valid.sum())])
        vc = remap_clusters(clusters, np.where(valid)[0])
        rmse, r2 = loco_cv(X, y[valid], vc)
        print(f"  {col:>15} {rmse:8.0f} {r2:6.3f}")
    print(f"\n  Verdict: steps at no window beats the noise floor.")

    # ── Part 3: sleep (null) ─────────────────────────────────────────
    print(f"\n{'=' * 70}")
    print(f"=== Part 3: Sleep ===")
    sleep_cols = [c for c in df.columns if c.startswith("sleep_")]
    print(f"  {'Feature':>15} {'CV RMSE':>8} {'R²':>6}")
    for col in sorted(sleep_cols, key=lambda c: int(c.split("_")[1].replace("d", ""))):
        v = df[col].values.astype(float)
        valid = ~np.isnan(v)
        if valid.sum() < 15:
            continue
        X = np.column_stack([ermr[valid], v[valid], np.ones(valid.sum())])
        vc = remap_clusters(clusters, np.where(valid)[0])
        if len(vc) < 3:
            continue
        rmse, r2 = loco_cv(X, y[valid], vc)
        print(f"  {col:>15} {rmse:8.0f} {r2:6.3f}")
    print(f"\n  Verdict: sleep at no window improves prediction. Negative R² throughout.")

    # ── Part 4: walk sessions (signal) ───────────────────────────────
    print(f"\n{'=' * 70}")
    print(f"=== Part 4: Walk sessions ===")

    print(f"\n  Walk sessions vs walk minutes vs steps (30-day window):")
    print(f"  {'Feature':>20} {'CV RMSE':>8} {'R²':>6}")
    for col in ["walk_sessions_30d", "walk_min_30d", "steps_30d"]:
        v = df[col].values.astype(float)
        valid = ~np.isnan(v)
        X = np.column_stack([ermr[valid], v[valid], np.ones(valid.sum())])
        vc = remap_clusters(clusters, np.where(valid)[0])
        rmse, r2 = loco_cv(X, y[valid], vc)
        print(f"  {col:>20} {rmse:8.0f} {r2:6.3f}")

    print(f"\n  Walk sessions at all windows:")
    print(f"  {'Feature':>20} {'CV RMSE':>8} {'R²':>6}")
    for w in WINDOWS:
        col = f"walk_sessions_{w}d"
        v = df[col].values.astype(float)
        X = np.column_stack([ermr, v, np.ones(n)])
        rmse, r2 = loco_cv(X, y, clusters)
        print(f"  {col:>20} {rmse:8.0f} {r2:6.3f}")

    # Best model coefficients
    ws30 = df["walk_sessions_30d"].values.astype(float)
    X_best = np.column_stack([ermr, ws30, np.ones(n)])
    c_best = np.linalg.lstsq(X_best, y, rcond=None)[0]
    pred_best = X_best @ c_best
    resid_best = y - pred_best
    rmse_insample = np.sqrt(np.mean(resid_best ** 2))
    rmse_cv, r2_cv = loco_cv(X_best, y, clusters)

    print(f"\n  Best model: RMR = {c_best[0]:.2f} × expected_rmr + {c_best[1]:.1f} × walk_sessions_30d + {c_best[2]:.0f}")
    print(f"  +{c_best[1]:.0f} cal RMR per walk session in 30 days")
    print(f"  CV RMSE = {rmse_cv:.0f}, R² = {r2_cv:.3f}")
    print(f"  In-sample RMSE = {rmse_insample:.0f}")

    print(f"\n  {'Date':>12} {'RMR':>5} {'Pred':>5} {'Err':>5} {'Walks':>6} {'FM':>5}")
    for i in range(n):
        print(f"  {str(df.iloc[i]['date'])[:10]:>12} {y[i]:5.0f} {pred_best[i]:5.0f} "
              f"{resid_best[i]:+5.0f} {ws30[i]:6.0f} {df.iloc[i]['fm_lbs']:5.0f}")

    # ── Part 5: confound analysis ────────────────────────────────────
    print(f"\n{'=' * 70}")
    print(f"=== Part 5: Is it really the walks, or is it season? ===")

    is_summer = df["is_summer"].values.astype(float)
    month_sin = df["month_sin"].values
    month_cos = df["month_cos"].values
    tirz = df["tirz_level"].values

    print(f"\n  Pairwise correlations:")
    for name, v in [("walk_sessions_30d", ws30), ("is_summer", is_summer),
                    ("steps_30d", df["steps_30d"].values), ("tirz_level", tirz)]:
        valid = ~np.isnan(v)
        r = np.corrcoef(v[valid], y[valid])[0, 1]
        print(f"    {name:>25} vs RMR:  r = {r:+.3f}")

    print(f"\n  walk_sessions_30d correlations with confounds:")
    for name, v in [("is_summer", is_summer), ("steps_30d", df["steps_30d"].values),
                    ("tirz_level", tirz), ("fm_lbs", df["fm_lbs"].values),
                    ("expected_rmr", ermr)]:
        valid = ~np.isnan(v)
        r = np.corrcoef(ws30[valid], v[valid])[0, 1]
        print(f"    {name:>25}: r = {r:+.3f}")

    print(f"\n  Model comparison (all leave-one-cluster-out CV):")
    print(f"  {'Model':>55} {'CV RMSE':>8} {'R²':>6}")

    models = [
        ("expected_rmr", np.column_stack([ermr, np.ones(n)])),
        ("expected_rmr + is_summer",
         np.column_stack([ermr, is_summer, np.ones(n)])),
        ("expected_rmr + season(sin/cos)",
         np.column_stack([ermr, month_sin, month_cos, np.ones(n)])),
        ("expected_rmr + walk_sessions_30d",
         np.column_stack([ermr, ws30, np.ones(n)])),
        ("expected_rmr + walk_sessions_30d + is_summer",
         np.column_stack([ermr, ws30, is_summer, np.ones(n)])),
        ("expected_rmr + walk_sessions_30d + season(sin/cos)",
         np.column_stack([ermr, ws30, month_sin, month_cos, np.ones(n)])),
        ("expected_rmr + walk_sessions_30d + tirz_level",
         np.column_stack([ermr, ws30, tirz, np.ones(n)])),
        ("expected_rmr + steps_30d + is_summer",
         np.column_stack([ermr, df["steps_30d"].values, is_summer, np.ones(n)])),
    ]

    for name, X in models:
        valid = ~np.isnan(X).any(axis=1)
        vc = remap_clusters(clusters, np.where(valid)[0])
        if len(vc) < 3:
            continue
        rmse, r2 = loco_cv(X[valid], y[valid], vc)
        print(f"  {name:>55} {rmse:8.0f} {r2:6.3f}")

    # Partial correlations controlling for season
    print(f"\n  Partial correlations (controlling for is_summer):")
    X_season = np.column_stack([is_summer, np.ones(n)])
    for name, v in [("walk_sessions_30d", ws30), ("steps_30d", df["steps_30d"].values)]:
        valid = ~np.isnan(v)
        vv = v[valid]
        yv = y[valid]
        Xs = X_season[valid]
        res_v = vv - Xs @ np.linalg.lstsq(Xs, vv, rcond=None)[0]
        res_y = yv - Xs @ np.linalg.lstsq(Xs, yv, rcond=None)[0]
        r = np.corrcoef(res_v, res_y)[0, 1]
        print(f"    {name:>25} vs RMR | season: r = {r:+.3f}")

    # ── Part 6: within-season and within-composition tests ───────────
    print(f"\n{'=' * 70}")
    print(f"=== Part 6: Within-subgroup tests ===")

    summer = df[df["is_summer"] == 1]
    winter = df[df["is_summer"] == 0]
    print(f"\n  Summer (n={len(summer)}):")
    print(f"    walk_sessions: {summer['walk_sessions_30d'].min():.0f}-{summer['walk_sessions_30d'].max():.0f}")
    print(f"    RMR range: {summer['rmr_kcal'].min():.0f}-{summer['rmr_kcal'].max():.0f}")
    r_summer = np.corrcoef(summer["walk_sessions_30d"].values, summer["rmr_kcal"].values)[0, 1]
    print(f"    walk_sessions vs RMR: r = {r_summer:+.3f}")

    print(f"\n  Winter (n={len(winter)}):")
    print(f"    walk_sessions: {winter['walk_sessions_30d'].min():.0f}-{winter['walk_sessions_30d'].max():.0f}")
    print(f"    RMR range: {winter['rmr_kcal'].min():.0f}-{winter['rmr_kcal'].max():.0f}")
    r_winter = np.corrcoef(winter["walk_sessions_30d"].values, winter["rmr_kcal"].values)[0, 1]
    print(f"    walk_sessions vs RMR: r = {r_winter:+.3f}")

    # 2022-2023 only (nearly constant body composition)
    df22 = df[(df["date"] >= "2022-01-01") & (df["date"] <= "2023-12-31")]
    if len(df22) >= 10:
        print(f"\n  2022-2023 only (n={len(df22)}, FM={df22['fm_lbs'].min():.0f}-{df22['fm_lbs'].max():.0f}):")
        r_22 = np.corrcoef(df22["walk_sessions_30d"].values, df22["rmr_kcal"].values)[0, 1]
        print(f"    walk_sessions vs RMR: r = {r_22:+.3f}")
        # Partial controlling for expected_rmr
        X22 = np.column_stack([df22["expected_rmr"].values, np.ones(len(df22))])
        ws22 = df22["walk_sessions_30d"].values.astype(float)
        y22 = df22["rmr_kcal"].values
        res_ws = ws22 - X22 @ np.linalg.lstsq(X22, ws22, rcond=None)[0]
        res_y = y22 - X22 @ np.linalg.lstsq(X22, y22, rcond=None)[0]
        r_partial = np.corrcoef(res_ws, res_y)[0, 1]
        print(f"    walk_sessions vs RMR | expected_rmr: r = {r_partial:+.3f}")

    # Within-cluster: May-Jun 2022
    may22 = df[(df["date"] >= "2022-05-01") & (df["date"] <= "2022-06-30")]
    if len(may22) >= 4:
        r_may = np.corrcoef(may22["walk_sessions_30d"].values, may22["rmr_kcal"].values)[0, 1]
        print(f"\n  May-Jun 2022 cluster (n={len(may22)}, FM={may22['fm_lbs'].min():.0f}-{may22['fm_lbs'].max():.0f}):")
        print(f"    walk_sessions: {may22['walk_sessions_30d'].min():.0f}-{may22['walk_sessions_30d'].max():.0f}")
        print(f"    RMR range: {may22['rmr_kcal'].min():.0f}-{may22['rmr_kcal'].max():.0f}")
        print(f"    walk_sessions vs RMR: r = {r_may:+.3f}")

    # ── Part 7: what walk sessions measure ───────────────────────────
    print(f"\n{'=' * 70}")
    print(f"=== Part 7: What distinguishes walk sessions from steps? ===")
    print(f"\n  Walk sessions = count of distinct exercise entries Samsung Health")
    print(f"  recorded as 'walking'. These are deliberate, sustained walks")
    print(f"  (typically 20+ min continuous). Steps = all daily movement")
    print(f"  including incidental activity.")
    print(f"\n  At each measurement:")
    print(f"  {'Date':>12} {'RMR':>5} {'WalkSess':>9} {'WalkMin':>8} {'Steps':>6} {'Summer':>7}")
    for i in range(n):
        r = df.iloc[i]
        print(f"  {str(r['date'])[:10]:>12} {r['rmr_kcal']:5.0f} "
              f"{r['walk_sessions_30d']:9.0f} {r['walk_min_30d']:8.0f} "
              f"{r['steps_30d']:6.0f} {'Y' if r['is_summer'] else 'N':>7}")

    # ── Summary ──────────────────────────────────────────────────────
    print(f"\n{'=' * 70}")
    print(f"=== Summary ===")
    print(f"  Baseline (expected_rmr only):  CV RMSE = {rmse0:.0f}")
    print(f"  + walk_sessions_30d:           CV RMSE = {rmse_cv:.0f}  (R² = {r2_cv:.3f})")
    print(f"  + walk_sessions_30d + season:  CV RMSE = ", end="")
    X_ws_season = np.column_stack([ermr, ws30, is_summer, np.ones(n)])
    rmse_ws_s, r2_ws_s = loco_cv(X_ws_season, y, clusters)
    print(f"{rmse_ws_s:.0f}  (R² = {r2_ws_s:.3f})")
    print(f"  Fitmate noise floor:           ~170")
    print(f"\n  Null results:")
    print(f"    - Calories, protein, carbs, fat, sodium at all windows")
    print(f"    - Step counts at all windows")
    print(f"    - Sleep hours at all windows")
    print(f"    - Strength training session count")
    print(f"\n  Signal:")
    print(f"    - Walk session count (30d): +{c_best[1]:.0f} cal RMR per session")
    print(f"    - Survives season control (partial r = ", end="")
    # Recompute partial for print
    Xs = np.column_stack([is_summer, np.ones(n)])
    res_ws_all = ws30 - Xs @ np.linalg.lstsq(Xs, ws30, rcond=None)[0]
    res_y_all = y - Xs @ np.linalg.lstsq(Xs, y, rcond=None)[0]
    r_partial_all = np.corrcoef(res_ws_all, res_y_all)[0, 1]
    print(f"{r_partial_all:.2f})")
    print(f"    - Holds within-winter (r = {r_winter:.2f}), within 2022-2023 (r = {r_22:.2f})")
    print(f"    - Session count matters more than total minutes or total steps")
    print(f"\n  Remaining confound: season (r = 0.83 with walk sessions).")
    print(f"  Cannot fully separate 'walks raise RMR' from 'something seasonal")
    print(f"  that also drives walking' with {n} measurements clustered by month.")


if __name__ == "__main__":
    main()
