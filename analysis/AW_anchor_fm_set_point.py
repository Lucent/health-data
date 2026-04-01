#!/usr/bin/env python3
"""AW. Non-circular set-point test using anchor-only daily fat mass.

Consumes AV_anchor_daily_composition.csv, which is built from:
  - corrected scale weight observations
  - direct body-composition anchors
  - smoothness priors

No intake, TDEE, or RMR enter the FM/FFM state interpolation.
"""

from pathlib import Path
import argparse
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
ABS_BINGE_THRESHOLD = 2800
MIN_TRAIN_ROWS = 365


def ema(series, half_life):
    alpha = 1 - np.exp(-np.log(2) / half_life)
    return pd.Series(series).ewm(alpha=alpha, min_periods=30).mean().values


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


def eval_linear_cv(df, predictors, outcome):
    valid = df.dropna(subset=predictors + [outcome]).copy()
    preds = []
    ys = []
    for year in sorted(valid["year"].unique()):
        train = valid[valid["year"] < year]
        test = valid[valid["year"] == year]
        if len(train) < MIN_TRAIN_ROWS or len(test) < 60:
            continue
        beta, x_mean, x_std = fit_linear(train[predictors].values, train[outcome].values)
        pred = predict_linear(test[predictors].values, beta, x_mean, x_std)
        preds.extend(pred.tolist())
        ys.extend(test[outcome].tolist())
    if len(ys) < 365:
        return None
    preds = np.asarray(preds)
    ys = np.asarray(ys)
    return {"r_pred": np.corrcoef(preds, ys)[0, 1], "rmse": np.sqrt(np.mean((preds - ys) ** 2)), "n": len(ys)}


def eval_auc_cv(df, predictors, outcome):
    valid = df.dropna(subset=predictors + [outcome]).copy()
    scores = []
    ys = []
    for year in sorted(valid["year"].unique()):
        train = valid[valid["year"] < year]
        test = valid[valid["year"] == year]
        if len(train) < MIN_TRAIN_ROWS or len(test) < 60:
            continue
        beta, x_mean, x_std = fit_linear(train[predictors].values, train[outcome].values)
        score = predict_linear(test[predictors].values, beta, x_mean, x_std)
        scores.extend(score.tolist())
        ys.extend(test[outcome].tolist())
    if len(ys) < 365 or sum(ys) < 10:
        return None
    return {"auc": roc_auc_score_np(np.asarray(ys), np.asarray(scores)), "n": len(ys), "events": int(sum(ys))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--body-file",
        default=str(ROOT / "analysis" / "AV_anchor_daily_composition.csv"),
        help="CSV with daily FM column and date",
    )
    parser.add_argument(
        "--fm-col",
        default="fm_lbs_anchor",
        help="Column name for daily fat mass in body-file",
    )
    args = parser.parse_args()

    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    body = pd.read_csv(args.body_file, parse_dates=["date"])

    daily = intake[["date", "calories"]].merge(body[["date", args.fm_col]], on="date", how="left")
    daily = daily.merge(tirz[["date", "effective_level"]], on="date", how="left")
    daily["effective_level"] = daily["effective_level"].fillna(0)
    daily["on_tirz"] = (daily["effective_level"] > 0).astype(int)
    daily["year"] = daily["date"].dt.year
    daily["binge_abs"] = (daily["calories"] > ABS_BINGE_THRESHOLD).astype(int)
    daily = daily.rename(columns={args.fm_col: "fm_lbs"})

    pre = daily[daily["on_tirz"] == 0].copy().reset_index(drop=True)

    # Future intake features.
    for horizon in [7, 14, 30]:
        pre[f"next_{horizon}d_cal"] = (
            pre["calories"].shift(-1).rolling(horizon, min_periods=horizon).mean().shift(-(horizon - 1))
        )
        pre[f"prev_{horizon}d_cal"] = pre["calories"].shift(1).rolling(horizon, min_periods=horizon).mean()
        pre[f"delta_{horizon}d_cal"] = pre[f"next_{horizon}d_cal"] - pre[f"prev_{horizon}d_cal"]

        future_binge = []
        vals = pre["binge_abs"].values
        n = len(pre)
        for i in range(n):
            j = i + 1
            k = min(n, i + 1 + horizon)
            future_binge.append(float(vals[j:k].max()) if k - j >= horizon else np.nan)
        pre[f"next_{horizon}d_binge"] = future_binge

    # FM state features.
    best_rows = []
    for hl in range(10, 241, 5):
        sp = ema(pre["fm_lbs"].values, hl)
        dist = sp - pre["fm_lbs"].values
        valid = ~np.isnan(dist) & pre["delta_30d_cal"].notna()
        if valid.sum() < 365:
            continue
        r = np.corrcoef(dist[valid], pre.loc[valid, "delta_30d_cal"])[0, 1]
        best_rows.append((hl, r))
    best_hl, best_corr = max(best_rows, key=lambda x: abs(x[1]))

    pre["sp_fm"] = ema(pre["fm_lbs"].values, best_hl)
    pre["sp_dist"] = pre["sp_fm"] - pre["fm_lbs"]
    pre["fm_delta_30d"] = pre["fm_lbs"] - pre["fm_lbs"].shift(30)

    print("=" * 70)
    print("ANCHOR-FM SET-POINT TEST")
    print("=" * 70)
    print(f"\nPre-tirzepatide daily rows: {len(pre)}")
    print(f"Best half-life for delta-30d intake: {best_hl}d  r={best_corr:+.4f}")

    valid30 = pre.dropna(subset=["sp_dist", "delta_30d_cal", "next_30d_cal"])
    r_next = np.corrcoef(valid30["sp_dist"], valid30["next_30d_cal"])[0, 1]
    r_shift = np.corrcoef(valid30["sp_dist"], valid30["delta_30d_cal"])[0, 1]
    print(f"Raw corr(sp_dist, next30_cal):  {r_next:+.4f}")
    print(f"Raw corr(sp_dist, delta30_cal): {r_shift:+.4f}")

    fixed_rows = []
    for fixed in np.arange(20, 110, 1):
        dist = fixed - valid30["fm_lbs"].values
        fixed_rows.append((fixed, np.corrcoef(dist, valid30["delta_30d_cal"].values)[0, 1]))
    best_fixed, best_fixed_r = max(fixed_rows, key=lambda x: abs(x[1]))
    print(f"Best fixed-FM target on delta30: {best_fixed:.0f} lbs  r={best_fixed_r:+.4f}")

    models = [
        ("SP distance", ["sp_dist"]),
        ("Absolute FM", ["fm_lbs"]),
        ("FM recent change", ["fm_delta_30d"]),
        ("Prev 30d calories", ["prev_30d_cal"]),
        ("SP + prev cal", ["sp_dist", "prev_30d_cal"]),
        ("FM + prev cal", ["fm_lbs", "prev_30d_cal"]),
    ]

    print("\n--- Future intake level prediction ---")
    for horizon in [7, 14, 30]:
        print(f"\nNext {horizon}d calories:")
        for label, predictors in models:
            res = eval_linear_cv(pre, predictors, f"next_{horizon}d_cal")
            if res:
                print(f"  {label:20s} r_pred={res['r_pred']:.4f}  RMSE={res['rmse']:.1f}  n={res['n']}")

    print("\n--- Intake shift prediction ---")
    for horizon in [7, 14, 30]:
        print(f"\nDelta {horizon}d calories:")
        for label, predictors in models[:-1]:
            res = eval_linear_cv(pre, predictors, f"delta_{horizon}d_cal")
            if res:
                print(f"  {label:20s} r_pred={res['r_pred']:.4f}  RMSE={res['rmse']:.1f}  n={res['n']}")

    print("\n--- Future binge prediction ---")
    for horizon in [7, 14, 30]:
        print(f"\nAny >{ABS_BINGE_THRESHOLD} cal in next {horizon}d:")
        for label, predictors in models:
            res = eval_auc_cv(pre, predictors, f"next_{horizon}d_binge")
            if res:
                print(f"  {label:20s} AUC={res['auc']:.4f}  events={res['events']}/{res['n']}")


if __name__ == "__main__":
    main()
