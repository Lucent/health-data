"""Binge prediction from set point distance.

Tests whether cumulative distance below a slowly-drifting set point
predicts binge probability (>2800 cal/day) better than dietary variables.

The tirzepatide natural experiment: at the same distance below set point,
do binges disappear on the drug?
"""

import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

BINGE_THRESHOLD = 2800
MIN_TRAIN_ROWS = 365
L2_PENALTY = 1.0


def load():
    kf = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    daily = intake.merge(
        kf[["date", "fat_mass_lbs_filtered", "fat_mass_lbs", "tdee_filtered", "tdee"]],
        on="date",
        how="left",
    )
    daily = daily.merge(tirz[["date", "effective_level"]], on="date", how="left")
    daily["effective_level"] = daily["effective_level"].fillna(0)
    return daily


def compute_features(daily):
    """Compute causal candidate binge predictors.

    Every predictive feature must be available before today's eating is known.
    """
    d = daily.copy()

    d["binge"] = (d["calories"] > BINGE_THRESHOLD).astype(int)
    d["fat_mass_for_prediction"] = d["fat_mass_lbs_filtered"].fillna(d["fat_mass_lbs"])
    d["tdee_for_prediction"] = d["tdee_filtered"].fillna(d["tdee"])

    # Set point estimates must be trailing-only and based on prior state.
    prev_fm = d["fat_mass_for_prediction"].shift(1)
    for window in [90, 180, 365]:
        d[f"sp_{window}d"] = prev_fm.rolling(window, min_periods=window // 2).mean()
        d[f"dist_{window}d"] = prev_fm - d[f"sp_{window}d"]

    # Cumulative deficit should only use prior days.
    d["deficit"] = d["calories"] - d["tdee_for_prediction"]
    prev_deficit = d["deficit"].shift(1)
    d["cum_deficit_7d"] = prev_deficit.rolling(7, min_periods=4).sum()
    d["cum_deficit_30d"] = prev_deficit.rolling(30, min_periods=15).sum()

    # Lagged intake and lagged macro composition.
    d["prev_cal"] = d["calories"].shift(1)
    d["prev_3d_cal"] = d["calories"].rolling(3).mean().shift(1)

    # Same-day macro percentages are partly tautological with a same-day binge.
    prot_pct = d["protein_g"] * 4 / d["calories"] * 100
    fat_pct = d["fat_g"] * 9 / d["calories"] * 100
    carb_pct = d["carbs_g"] * 4 / d["calories"] * 100
    d["prev_prot_pct"] = prot_pct.shift(1)
    d["prev_fat_pct"] = fat_pct.shift(1)
    d["prev_carb_pct"] = carb_pct.shift(1)

    # Day of week (weekend binges?)
    d["dow"] = d["date"].dt.dayofweek
    d["weekend"] = (d["dow"] >= 5).astype(int)

    # Season (quarter)
    d["quarter"] = d["date"].dt.quarter

    # On drug
    d["on_tirz"] = (d["effective_level"] > 0).astype(int)

    return d


def roc_auc_score_np(y_true, scores):
    """Rank-based ROC AUC without external dependencies."""
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


def fit_logistic_regression(X, y, l2_penalty=L2_PENALTY, max_iter=50):
    """Fit logistic regression with IRLS using only numpy."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    X_design = np.column_stack([np.ones(len(X)), X])
    beta = np.zeros(X_design.shape[1], dtype=float)
    penalty = np.eye(X_design.shape[1]) * l2_penalty
    penalty[0, 0] = 0.0

    for _ in range(max_iter):
        eta = X_design @ beta
        p = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
        w = np.clip(p * (1.0 - p), 1e-6, None)
        z = eta + (y - p) / w
        XtW = X_design.T * w
        hessian = XtW @ X_design + penalty
        rhs = XtW @ z
        beta_new = np.linalg.solve(hessian, rhs)
        if np.max(np.abs(beta_new - beta)) < 1e-6:
            beta = beta_new
            break
        beta = beta_new

    return beta


def predict_logistic_regression(X, beta):
    X = np.asarray(X, dtype=float)
    X_design = np.column_stack([np.ones(len(X)), X])
    eta = X_design @ beta
    return 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))


def standardize_train_test(train_X, test_X):
    mean = train_X.mean(axis=0)
    std = train_X.std(axis=0)
    std[std == 0] = 1.0
    return (train_X - mean) / std, (test_X - mean) / std


def evaluate_predictor(d, predictor_cols, label="", subset=None):
    """Walk-forward yearly AUC for binge prediction."""
    if subset is not None:
        d = d[subset]
    valid = d.dropna(subset=predictor_cols + ["binge"])
    if valid["binge"].sum() < 10 or len(valid) < 100:
        return None

    years = sorted(valid["date"].dt.year.unique())
    probs = []
    ys = []
    last_beta = None

    for year in years:
        train = valid[valid["date"].dt.year < year]
        test = valid[valid["date"].dt.year == year]
        if len(train) < MIN_TRAIN_ROWS or len(test) < 30:
            continue
        if train["binge"].sum() == 0 or test["binge"].sum() == 0:
            continue

        train_X = train[predictor_cols].values
        test_X = test[predictor_cols].values
        train_X_scaled, test_X_scaled = standardize_train_test(train_X, test_X)
        beta = fit_logistic_regression(train_X_scaled, train["binge"].values)
        last_beta = beta
        probs.extend(predict_logistic_regression(test_X_scaled, beta))
        ys.extend(test["binge"].values)

    if len(ys) < 100 or sum(ys) < 10:
        return None

    auc = roc_auc_score_np(ys, probs)
    coeffs = {}
    if last_beta is not None:
        coeffs = dict(zip(predictor_cols, last_beta[1:]))

    return {
        "label": label,
        "auc": auc,
        "n": len(ys),
        "n_binge": int(sum(ys)),
        "rate": np.mean(ys),
        "coeffs": coeffs,
    }


def main():
    print("Loading data...")
    daily = load()
    d = compute_features(daily)

    valid = d.dropna(subset=["fat_mass_for_prediction", "dist_180d", "prev_cal", "prev_prot_pct"])
    print(f"Days with all features: {len(valid)}")
    print(f"Binges: {valid['binge'].sum()} ({valid['binge'].mean()*100:.1f}%)")

    # === 1. Which set point window works best? ===
    print(f"\n=== Set point window comparison (walk-forward AUC) ===")
    for window in [90, 180, 365]:
        r = evaluate_predictor(valid, [f"dist_{window}d"], f"Distance ({window}d SP)")
        if r:
            print(f"  {r['label']:30s}  AUC={r['auc']:.4f}  n={r['n']}")

    # === 2. Compare predictor groups ===
    print(f"\n=== Predictor comparison (walk-forward AUC) ===")

    groups = [
        ("Distance from 180d SP", ["dist_180d"]),
        ("Cum deficit 30d", ["cum_deficit_30d"]),
        ("Yesterday's calories", ["prev_cal"]),
        ("Trailing 3d calories", ["prev_3d_cal"]),
        ("Yesterday's protein %", ["prev_prot_pct"]),
        ("Weekend", ["weekend"]),
        ("Fat mass (absolute)", ["fat_mass_for_prediction"]),
        ("Tirzepatide level", ["effective_level"]),
        ("Distance + deficit", ["dist_180d", "cum_deficit_30d"]),
        ("All dietary", ["prev_cal", "prev_3d_cal", "prev_prot_pct", "prev_fat_pct", "weekend"]),
        ("Distance + dietary", ["dist_180d", "cum_deficit_30d", "prev_cal", "prev_prot_pct"]),
        ("Kitchen sink", ["dist_180d", "cum_deficit_30d", "prev_cal", "prev_3d_cal",
                          "prev_prot_pct", "prev_fat_pct", "weekend", "fat_mass_for_prediction"]),
    ]

    for label, cols in groups:
        r = evaluate_predictor(valid, cols, label)
        if r:
            print(f"  {r['label']:30s}  AUC={r['auc']:.4f}  n={r['n']}  binges={r['n_binge']}")

    # === 3. Binge rate by distance from set point ===
    print(f"\n=== Binge rate by distance from 180d set point ===")
    print(f"(Negative = below set point, positive = above)")
    for lo, hi in [(-30, -10), (-10, -5), (-5, -2), (-2, 0), (0, 2), (2, 5), (5, 10), (10, 30)]:
        m = (valid["dist_180d"] >= lo) & (valid["dist_180d"] < hi)
        if m.sum() > 30:
            rate = valid.loc[m, "binge"].mean() * 100
            cal = valid.loc[m, "calories"].mean()
            print(f"  {lo:+3d} to {hi:+3d} lbs: binge rate={rate:5.1f}%  "
                  f"mean cal={cal:.0f}  n={m.sum()}")

    # === 4. The tirzepatide natural experiment ===
    print(f"\n=== Tirzepatide: same distance, different binge rate? ===")
    pre = valid[valid["on_tirz"] == 0]
    post = valid[valid["on_tirz"] == 1]

    print(f"\nPre-tirzepatide:")
    for lo, hi in [(-10, -2), (-2, 2), (2, 10)]:
        m = (pre["dist_180d"] >= lo) & (pre["dist_180d"] < hi)
        if m.sum() > 20:
            rate = pre.loc[m, "binge"].mean() * 100
            print(f"  Distance {lo:+d} to {hi:+d}: binge rate={rate:.1f}%  n={m.sum()}")

    print(f"\nOn tirzepatide:")
    for lo, hi in [(-10, -2), (-2, 2), (2, 10)]:
        m = (post["dist_180d"] >= lo) & (post["dist_180d"] < hi)
        if m.sum() > 10:
            rate = post.loc[m, "binge"].mean() * 100
            print(f"  Distance {lo:+d} to {hi:+d}: binge rate={rate:.1f}%  n={m.sum()}")

    # === 5. Logistic regression with tirzepatide interaction ===
    print(f"\n=== Distance × tirzepatide interaction ===")
    r_pre = evaluate_predictor(valid, ["dist_180d"], "Distance (pre-tirz)",
                               subset=valid["on_tirz"] == 0)
    r_post = evaluate_predictor(valid, ["dist_180d"], "Distance (on-tirz)",
                                subset=valid["on_tirz"] == 1)
    if r_pre:
        print(f"  Pre:  AUC={r_pre['auc']:.4f}  binges={r_pre['n_binge']}/{r_pre['n']}")
    if r_post:
        print(f"  Post: AUC={r_post['auc']:.4f}  binges={r_post['n_binge']}/{r_post['n']}")

    # Combined model with drug
    r_combined = evaluate_predictor(valid, ["dist_180d", "effective_level"],
                                    "Distance + tirz level")
    if r_combined:
        print(f"  Combined: AUC={r_combined['auc']:.4f}")
        for k, v in r_combined["coeffs"].items():
            print(f"    {k}: β={v:+.4f}")


if __name__ == "__main__":
    main()
