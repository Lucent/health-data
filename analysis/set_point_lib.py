"""Shared helpers for set-point comparison artifacts."""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
ABS_BINGE_THRESHOLD = 2800


def ema(series, half_life):
    alpha = 1 - np.exp(-np.log(2) / half_life)
    return pd.Series(series).ewm(alpha=alpha, min_periods=30).mean().values


def timed_ema(values, dates, half_life_days):
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


def build_daily_base():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    daily = intake[["date", "calories"]].merge(
        tirz[["date", "effective_level"]], on="date", how="left"
    )
    daily["effective_level"] = daily["effective_level"].fillna(0)
    daily["on_tirz"] = (daily["effective_level"] > 0).astype(int)
    daily["binge_abs"] = (daily["calories"] > ABS_BINGE_THRESHOLD).astype(int)
    return daily


def add_window_features(daily, horizon):
    d = daily.copy()
    d[f"next_{horizon}d_cal"] = (
        d["calories"].shift(-1).rolling(horizon, min_periods=horizon).mean().shift(-(horizon - 1))
    )
    d[f"prev_{horizon}d_cal"] = d["calories"].shift(1).rolling(horizon, min_periods=horizon).mean()
    d[f"delta_{horizon}d_cal"] = d[f"next_{horizon}d_cal"] - d[f"prev_{horizon}d_cal"]
    future_binge = []
    vals = d["binge_abs"].astype(float).values
    n = len(d)
    for i in range(n):
        j = i + 1
        k = min(n, i + 1 + horizon)
        future_binge.append(float(vals[j:k].max()) if k - j >= horizon else np.nan)
    d[f"next_{horizon}d_binge"] = future_binge
    return d


def add_rolling_outcomes(daily):
    d = daily.copy()
    d["mean_90d_cal"] = d["calories"].rolling(90, min_periods=90).mean()
    return d


def load_state_daily(path, col):
    body = pd.read_csv(path, parse_dates=["date"])
    base = build_daily_base()
    return body[["date", col]].rename(columns={col: "state"}).merge(base, on="date", how="left")


def evaluate_fixed_outcome_half_life(df, state_col, outcome, irregular=False, half_life_grid=range(10, 241, 5)):
    pre = df[df["on_tirz"] == 0].copy().reset_index(drop=True)
    best_hl = None
    best_r = None
    r40 = r45 = r50 = np.nan
    for hl in half_life_grid:
        sp = timed_ema(pre[state_col].values, pre["date"].values, hl) if irregular else ema(pre[state_col].values, hl)
        dist = sp - pre[state_col].values
        valid = ~np.isnan(dist) & pre[outcome].notna()
        if valid.sum() < 365:
            continue
        r = np.corrcoef(dist[valid], pre.loc[valid, outcome])[0, 1]
        if hl == 40:
            r40 = r
        if hl == 45:
            r45 = r
        if hl == 50:
            r50 = r
        if best_r is None or abs(r) > abs(best_r):
            best_hl = hl
            best_r = r
    return {
        "best_hl": best_hl,
        "best_r": best_r,
        "r40": r40,
        "r45": r45,
        "r50": r50,
    }


def evaluate_primary_signal(df, state_col, model_name, irregular=False, half_life_grid=range(10, 241, 5)):
    pre = df[df["on_tirz"] == 0].copy().reset_index(drop=True)
    pre = add_window_features(pre, 30)
    best = evaluate_fixed_outcome_half_life(pre, state_col, "delta_30d_cal", irregular=irregular, half_life_grid=half_life_grid)
    sp = timed_ema(pre[state_col].values, pre["date"].values, best["best_hl"]) if irregular else ema(pre[state_col].values, best["best_hl"])
    pre["sp_dist"] = sp - pre[state_col].values
    valid = pre.dropna(subset=["sp_dist", "delta_30d_cal", "next_30d_cal", "next_30d_binge"]).copy()
    return {
        "model": model_name,
        "best_hl_days": best["best_hl"],
        "corr_spdist_delta30": np.corrcoef(valid["sp_dist"], valid["delta_30d_cal"])[0, 1],
        "corr_spdist_next30": np.corrcoef(valid["sp_dist"], valid["next_30d_cal"])[0, 1],
        "auc_next30_binge": roc_auc_score_np(valid["next_30d_binge"].values, valid["sp_dist"].values),
        "n_pre": len(pre),
    }
