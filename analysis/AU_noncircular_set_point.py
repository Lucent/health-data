#!/usr/bin/env python3
"""AU. Non-circular set-point test using observed glycogen-smoothed weight only.

Goal:
  Steelman the set-point idea while removing the strongest circularity:
  - no Kalman TDEE
  - no RMR
  - no intake-fed weight interpolation

Method:
  1. Use only observed morning weigh-ins, corrected for glycogen/sodium water
     via P1_smoothed_weight.csv.
  2. Build a trailing EMA "defended state" from those observed weights only.
  3. On each weigh-in date, predict future intake over the next 7/14/30 days.
  4. Compare moving-state distance vs absolute weight, recent weight change,
     and a fixed defended weight.

This is not a full physiological identification test. It is a stronger,
less endogenous intake-side test than the surplus-vs-fat-mass analyses.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
ABS_BINGE_THRESHOLD = 2800
MIN_TRAIN_ROWS = 120


def timed_ema(values, dates, half_life_days):
    """Trailing EMA over irregularly spaced observations."""
    out = np.empty(len(values))
    out[:] = np.nan
    if len(values) == 0:
        return out
    out[0] = values[0]
    for i in range(1, len(values)):
        dt = (dates[i] - dates[i - 1]) / np.timedelta64(1, "D")
        alpha = 1 - np.exp(-np.log(2) * dt / half_life_days)
        out[i] = out[i - 1] + alpha * (values[i] - out[i - 1])
    return out


def fit_linear(train_x, train_y):
    x = np.asarray(train_x, dtype=float)
    y = np.asarray(train_y, dtype=float)
    x_mean = x.mean(axis=0)
    x_std = x.std(axis=0)
    x_std[x_std == 0] = 1.0
    xz = (x - x_mean) / x_std
    X = np.column_stack([xz, np.ones(len(xz))])
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    return beta, x_mean, x_std


def predict_linear(x, beta, x_mean, x_std):
    x = np.asarray(x, dtype=float)
    xz = (x - x_mean) / x_std
    X = np.column_stack([xz, np.ones(len(xz))])
    return X @ beta


def roc_auc_score_np(y_true, scores):
    y_true = np.asarray(y_true, dtype=float)
    scores = np.asarray(scores, dtype=float)
    pos = y_true == 1
    neg = y_true == 0
    n_pos = pos.sum()
    n_neg = neg.sum()
    if n_pos == 0 or n_neg == 0:
        return np.nan
    order = np.argsort(scores)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(scores) + 1)
    _, inv, counts = np.unique(scores, return_inverse=True, return_counts=True)
    for group in np.where(counts > 1)[0]:
        idx = np.where(inv == group)[0]
        ranks[idx] = ranks[idx].mean()
    rank_sum = ranks[pos].sum()
    return (rank_sum - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def load_data():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    weight = pd.read_csv(ROOT / "analysis" / "P1_smoothed_weight.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    tirz = tirz[["date", "effective_level"]]
    tirz["effective_level"] = tirz["effective_level"].fillna(0)
    return intake, weight, tirz


def build_dataset():
    intake, weight, tirz = load_data()

    daily = intake[["date", "calories"]].merge(tirz, on="date", how="left")
    daily["effective_level"] = daily["effective_level"].fillna(0)
    daily["on_tirz"] = (daily["effective_level"] > 0).astype(int)
    daily["binge_abs"] = (daily["calories"] > ABS_BINGE_THRESHOLD).astype(int)
    daily["year"] = daily["date"].dt.year

    obs = weight[["date", "smoothed_weight_lbs"]].dropna().sort_values("date").reset_index(drop=True)
    obs = obs.rename(columns={"smoothed_weight_lbs": "body_state_lbs"})
    obs = obs.merge(daily[["date", "effective_level"]], on="date", how="left")
    obs["effective_level"] = obs["effective_level"].fillna(0)
    obs["on_tirz"] = (obs["effective_level"] > 0).astype(int)
    obs["year"] = obs["date"].dt.year

    # Strictly trailing recent change features.
    obs["prev_state"] = obs["body_state_lbs"].shift(1)
    obs["delta_30_obs"] = obs["body_state_lbs"] - obs["body_state_lbs"].shift(30)
    obs["delta_10_obs"] = obs["body_state_lbs"] - obs["body_state_lbs"].shift(10)

    for horizon in [7, 14, 30]:
        future = []
        binge = []
        prior = []
        for d in obs["date"]:
            next_days = daily[(daily["date"] > d) & (daily["date"] <= d + pd.Timedelta(days=horizon))]
            prev_days = daily[(daily["date"] < d) & (daily["date"] >= d - pd.Timedelta(days=horizon))]
            if len(next_days) >= horizon - 1:
                future.append(next_days["calories"].mean())
                binge.append(float(next_days["binge_abs"].any()))
            else:
                future.append(np.nan)
                binge.append(np.nan)
            prior.append(prev_days["calories"].mean() if len(prev_days) >= horizon - 1 else np.nan)
        obs[f"next_{horizon}d_cal"] = future
        obs[f"next_{horizon}d_binge"] = binge
        obs[f"prev_{horizon}d_cal"] = prior
        obs[f"delta_{horizon}d_cal"] = obs[f"next_{horizon}d_cal"] - obs[f"prev_{horizon}d_cal"]

    return obs


def eval_linear_cv(df, predictors, outcome):
    valid = df.dropna(subset=predictors + [outcome]).copy()
    preds = []
    ys = []
    for year in sorted(valid["year"].unique()):
        train = valid[valid["year"] < year]
        test = valid[valid["year"] == year]
        if len(train) < MIN_TRAIN_ROWS or len(test) < 20:
            continue
        beta, x_mean, x_std = fit_linear(train[predictors].values, train[outcome].values)
        pred = predict_linear(test[predictors].values, beta, x_mean, x_std)
        preds.extend(pred.tolist())
        ys.extend(test[outcome].tolist())
    if len(ys) < 80:
        return None
    preds = np.asarray(preds)
    ys = np.asarray(ys)
    r = np.corrcoef(preds, ys)[0, 1]
    rmse = np.sqrt(np.mean((preds - ys) ** 2))
    return {"r_pred": r, "rmse": rmse, "n": len(ys)}


def eval_auc_cv(df, predictors, outcome):
    valid = df.dropna(subset=predictors + [outcome]).copy()
    scores = []
    ys = []
    for year in sorted(valid["year"].unique()):
        train = valid[valid["year"] < year]
        test = valid[valid["year"] == year]
        if len(train) < MIN_TRAIN_ROWS or len(test) < 20:
            continue
        beta, x_mean, x_std = fit_linear(train[predictors].values, train[outcome].values)
        score = predict_linear(test[predictors].values, beta, x_mean, x_std)
        scores.extend(score.tolist())
        ys.extend(test[outcome].tolist())
    if len(ys) < 80 or sum(ys) < 10:
        return None
    auc = roc_auc_score_np(np.asarray(ys), np.asarray(scores))
    return {"auc": auc, "n": len(ys), "events": int(sum(ys))}


def main():
    obs = build_dataset()
    pre = obs[obs["on_tirz"] == 0].copy().reset_index(drop=True)

    print("=" * 70)
    print("NON-CIRCULAR SET-POINT TEST")
    print("=" * 70)
    print(f"\nObserved weigh-ins with glycogen smoothing only: {len(obs)}")
    print(f"Pre-tirzepatide weigh-ins: {len(pre)}")

    # Half-life sweep for future intake
    print("\n--- Moving state half-life sweep: distance -> next 30d calories ---")
    print(f"{'HL':>6} {'r(dist,outcome)':>16} {'partial|wt':>12}")

    sweep = []
    for hl in range(10, 241, 5):
        ema = timed_ema(pre["body_state_lbs"].values, pre["date"].values, hl)
        dist = ema - pre["body_state_lbs"].values
        sub = pre.copy()
        sub["sp_ema"] = ema
        sub["sp_dist"] = dist
        valid = sub.dropna(subset=["sp_dist", "next_30d_cal"])
        if len(valid) < 150:
            continue
        r = np.corrcoef(valid["sp_dist"], valid["next_30d_cal"])[0, 1]
        X = np.column_stack([valid["body_state_lbs"].values, np.ones(len(valid))])
        res_d = valid["sp_dist"].values - X @ np.linalg.lstsq(X, valid["sp_dist"].values, rcond=None)[0]
        res_y = valid["next_30d_cal"].values - X @ np.linalg.lstsq(X, valid["next_30d_cal"].values, rcond=None)[0]
        pr = np.corrcoef(res_d, res_y)[0, 1]
        sweep.append((hl, r, pr))
        if hl % 20 == 0 or hl in [45, 50, 55]:
            print(f"{hl:6d} {r:16.4f} {pr:12.4f}")

    best_hl, best_r, best_pr = max(sweep, key=lambda x: abs(x[1]))
    pre["sp_ema"] = timed_ema(pre["body_state_lbs"].values, pre["date"].values, best_hl)
    pre["sp_dist"] = pre["sp_ema"] - pre["body_state_lbs"]

    # Fixed-weight comparison
    fixed_rows = []
    valid30 = pre.dropna(subset=["next_30d_cal"])
    for fixed in np.arange(160, 270, 1):
        dist = fixed - valid30["body_state_lbs"].values
        r = np.corrcoef(dist, valid30["next_30d_cal"].values)[0, 1]
        fixed_rows.append((fixed, r))
    best_fixed, best_fixed_r = max(fixed_rows, key=lambda x: abs(x[1]))

    print(f"\nBest moving HL: {best_hl}d  r={best_r:+.4f}  partial|weight={best_pr:+.4f}")
    print(f"Best fixed weight: {best_fixed:.0f} lbs  r={best_fixed_r:+.4f}")

    print("\n--- Future intake prediction from observed weigh-ins only ---")
    models = [
        ("SP distance", ["sp_dist"]),
        ("Absolute weight", ["body_state_lbs"]),
        ("Recent 30-ob diff", ["delta_30_obs"]),
        ("Prev 30d calories", ["prev_30d_cal"]),
        ("SP + prev cal", ["sp_dist", "prev_30d_cal"]),
        ("Weight + prev cal", ["body_state_lbs", "prev_30d_cal"]),
        ("SP + weight + prev cal", ["sp_dist", "body_state_lbs", "prev_30d_cal"]),
    ]
    for horizon in [7, 14, 30]:
        print(f"\nNext {horizon}d mean calories:")
        for label, preds in models:
            res = eval_linear_cv(pre, preds, f"next_{horizon}d_cal")
            if res:
                print(f"  {label:24s} r_pred={res['r_pred']:.4f}  RMSE={res['rmse']:.1f}  n={res['n']}")

    print("\n--- Intake shift prediction (future minus prior window) ---")
    shift_models = [
        ("SP distance", ["sp_dist"]),
        ("Absolute weight", ["body_state_lbs"]),
        ("Recent 30-ob diff", ["delta_30_obs"]),
        ("SP + weight", ["sp_dist", "body_state_lbs"]),
    ]
    for horizon in [7, 14, 30]:
        valid = pre.dropna(subset=["sp_dist", "body_state_lbs", f"delta_{horizon}d_cal"]).copy()
        r = np.corrcoef(valid["sp_dist"], valid[f"delta_{horizon}d_cal"])[0, 1]
        print(f"\nDelta {horizon}d calories: raw corr(sp_dist, delta)={r:+.4f}")
        for label, preds in shift_models:
            res = eval_linear_cv(pre, preds, f"delta_{horizon}d_cal")
            if res:
                print(f"  {label:24s} r_pred={res['r_pred']:.4f}  RMSE={res['rmse']:.1f}  n={res['n']}")

    print("\n--- Future binge prediction (absolute threshold 2800 cal) ---")
    for horizon in [7, 14, 30]:
        print(f"\nAny >{ABS_BINGE_THRESHOLD} cal in next {horizon}d:")
        for label, preds in models:
            res = eval_auc_cv(pre, preds, f"next_{horizon}d_binge")
            if res:
                print(f"  {label:24s} AUC={res['auc']:.4f}  events={res['events']}/{res['n']}")

    print("\n--- Directional summaries at best half-life ---")
    valid = pre.dropna(subset=["sp_dist", "next_30d_cal", "next_30d_binge", "delta_30d_cal"]).copy()
    valid["band"] = pd.cut(valid["sp_dist"], bins=[-50, -5, -2, 0, 2, 5, 50])
    for band, g in valid.groupby("band", observed=True):
        if len(g) < 20:
            continue
        print(
            f"  {str(band):18s} n={len(g):4d}  "
            f"next30_cal={g['next_30d_cal'].mean():.0f}  "
            f"delta30_cal={g['delta_30d_cal'].mean():+.0f}  "
            f"next30_binge={g['next_30d_binge'].mean()*100:.1f}%"
        )


if __name__ == "__main__":
    main()
