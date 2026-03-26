"""Search for a hidden control-capacity model that predicts next-day binges.

The model treats "willpower" as a latent control resource: recent restriction,
monotony, travel, sleep loss, and similar pressures deplete it; easier days
partially restore it. We search for the best short-memory latent model and
compare it with the obvious baseline of yesterday's calories.

Outputs:
    analysis/latent_control_model_summary.csv
    analysis/latent_control_feature_search.csv
    analysis/latent_control_daily.csv
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.append(str(ROOT))

from C_binge_analysis import (  # noqa: E402
    fit_logistic_regression,
    predict_logistic_regression,
    roc_auc_score_np,
    standardize_train_test,
)

BINGE_THRESHOLD = 2800
START_DATE = pd.Timestamp("2017-01-01")
MAX_STRAIN_FEATURES = 4
HALF_LIFE_CANDIDATES = [2, 3, 5, 7, 10, 14, 21]
MIN_TRAIN_ROWS = 365

STRAIN_FEATURES = [
    "restriction_pressure",
    "low_protein_pressure",
    "sleep_debt",
    "high_steps",
    "travel_flag",
    "low_calorie_day",
    "restriction_streak",
    "monotony",
]
RECOVERY_FEATURES = ["recovery_ease", "binge_recovery"]


def parse_travel_windows() -> pd.DataFrame:
    rows = []
    for line in (ROOT / "travel" / "trips.md").read_text().splitlines():
        line = line.strip()
        if not (line.startswith("| 201") or line.startswith("| 202")):
            continue
        parts = [part.strip() for part in line.strip("|").split("|")]
        rows.append(
            {
                "start": pd.Timestamp(parts[0]),
                "end": pd.Timestamp(parts[1]),
                "days": int(parts[2]),
                "location": parts[3],
            }
        )
    return pd.DataFrame(rows)


def load_data() -> pd.DataFrame:
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    steps = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])
    sleep = pd.read_csv(ROOT / "steps-sleep" / "sleep.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    travel = parse_travel_windows()

    daily = intake.merge(kalman[["date", "tdee_filtered", "tdee"]], on="date", how="left")
    daily = daily.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    daily = daily.merge(steps[["date", "steps"]], on="date", how="left")
    daily = daily.merge(sleep[["date", "sleep_hours"]], on="date", how="left")
    daily = daily.merge(tirz[["date", "effective_level"]], on="date", how="left")
    daily["effective_level"] = daily["effective_level"].fillna(0)
    daily["on_tirz"] = (daily["effective_level"] > 0).astype(int)
    daily["tdee_for_prediction"] = daily["tdee_filtered"].fillna(daily["tdee"])
    daily["binge"] = (daily["calories"] > BINGE_THRESHOLD).astype(int)
    daily["protein_pct"] = 4 * daily["protein_g"] / daily["calories"]
    daily["deficit"] = daily["tdee_for_prediction"] - daily["calories"]
    daily["travel"] = 0
    for _, row in travel.iterrows():
        mask = (daily["date"] >= row["start"]) & (daily["date"] <= row["end"])
        daily.loc[mask, "travel"] = 1
    return daily.sort_values("date").reset_index(drop=True)


def compute_sources(daily: pd.DataFrame) -> pd.DataFrame:
    d = daily.copy()
    d["prev_cal"] = d["calories"].shift(1)
    d["restriction_pressure"] = np.maximum(0, d["deficit"].shift(1) / 500)
    d["low_protein_pressure"] = np.maximum(0, (0.16 - d["protein_pct"].shift(1)) / 0.04)
    d["sleep_debt"] = np.maximum(0, (7.5 - d["sleep_hours"].shift(1)) / 1.5)
    d["high_steps"] = np.maximum(0, (d["steps"].shift(1) - 6000) / 4000)
    d["travel_flag"] = d["travel"].shift(1).fillna(0)
    d["low_calorie_day"] = np.maximum(0, (2200 - d["calories"].shift(1)) / 400)
    d["binge_prev"] = (d["calories"].shift(1) > BINGE_THRESHOLD).astype(float)

    restrict_prev = (d["calories"].shift(1) < 1800).fillna(False)
    streak = []
    current = 0
    for value in restrict_prev:
        current = current + 1 if value else 0
        streak.append(current)
    d["restriction_streak"] = np.array(streak) / 3

    trailing_std = d["calories"].rolling(7, min_periods=4).std().shift(1)
    d["monotony"] = np.maximum(0, (500 - trailing_std) / 250)
    d["recovery_ease"] = np.maximum(0, (d["prev_cal"] - 2200) / 600)
    d["binge_recovery"] = d["binge_prev"]
    return d


def add_ewm_features(df: pd.DataFrame, half_life: int) -> pd.DataFrame:
    out = df.copy()
    alpha = 1 - np.exp(-np.log(2) / half_life)
    for feature in STRAIN_FEATURES + RECOVERY_FEATURES:
        out[f"ew_{feature}"] = out[feature].fillna(0).ewm(alpha=alpha, adjust=False).mean()
    return out


def evaluate_predictor(df: pd.DataFrame, cols: list[str]) -> tuple[float, int, int]:
    valid = df.dropna(subset=cols + ["binge"]).copy()
    probs = []
    ys = []
    years = sorted(valid["date"].dt.year.unique())

    for year in years:
        train = valid[valid["date"].dt.year < year]
        test = valid[valid["date"].dt.year == year]
        if len(train) < MIN_TRAIN_ROWS or len(test) < 30:
            continue
        if train["binge"].sum() == 0 or test["binge"].sum() == 0:
            continue

        train_x, test_x = standardize_train_test(train[cols].values, test[cols].values)
        beta = fit_logistic_regression(train_x, train["binge"].values)
        probs.extend(predict_logistic_regression(test_x, beta))
        ys.extend(test["binge"].values)

    return roc_auc_score_np(np.array(ys), np.array(probs)), len(ys), int(sum(ys))


def forward_select(df: pd.DataFrame) -> list[str]:
    chosen = []
    remaining = [f"ew_{feature}" for feature in STRAIN_FEATURES]
    current_auc = 0.0

    for _ in range(MAX_STRAIN_FEATURES):
        best_auc = current_auc
        best_feature = None
        for feature in remaining:
            auc, _, _ = evaluate_predictor(df, chosen + [feature])
            if auc > best_auc + 1e-6:
                best_auc = auc
                best_feature = feature
        if best_feature is None:
            break
        chosen.append(best_feature)
        remaining.remove(best_feature)
        current_auc = best_auc

    return chosen


def fit_search(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    rows = []
    best = None

    for half_life in HALF_LIFE_CANDIDATES:
        tmp = add_ewm_features(df, half_life)
        chosen = forward_select(tmp)

        for recovery_set in [[], ["ew_recovery_ease"], ["ew_binge_recovery"], ["ew_recovery_ease", "ew_binge_recovery"]]:
            cols = chosen + recovery_set
            auc, n_rows, n_binges = evaluate_predictor(tmp, cols)
            row = {
                "half_life_days": half_life,
                "latent_auc": round(auc, 4),
                "n_rows": n_rows,
                "n_binges": n_binges,
                "strain_features": ",".join(chosen),
                "recovery_features": ",".join(recovery_set),
                "all_features": ",".join(cols),
            }
            rows.append(row)
            if best is None or auc > best["latent_auc"]:
                best = {
                    "half_life_days": half_life,
                    "latent_auc": auc,
                    "cols": cols,
                    "strain_features": chosen,
                    "recovery_features": recovery_set,
                }

    return pd.DataFrame(rows).sort_values(["latent_auc", "half_life_days"], ascending=[False, True]), best


def fit_final_score(df: pd.DataFrame, cols: list[str]) -> tuple[np.ndarray, np.ndarray]:
    valid = df.dropna(subset=cols + ["binge"]).copy()
    x = valid[cols].values
    x_scaled, _ = standardize_train_test(x, x)
    beta = fit_logistic_regression(x_scaled, valid["binge"].values)
    score = x_scaled @ beta[1:]
    return score, beta[1:]


def build_outputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    daily = compute_sources(load_data())
    daily = daily[(daily["on_tirz"] == 0) & (daily["date"] >= START_DATE)].copy()
    daily = daily.dropna(subset=["prev_cal", "sleep_hours"]).reset_index(drop=True)

    search_df, best = fit_search(daily)
    best_daily = add_ewm_features(daily, int(best["half_life_days"]))

    baseline_auc, n_rows, n_binges = evaluate_predictor(best_daily, ["prev_cal"])
    latent_auc, _, _ = evaluate_predictor(best_daily, best["cols"])
    combo_auc, _, _ = evaluate_predictor(best_daily, ["prev_cal"] + best["cols"])

    score, beta = fit_final_score(best_daily, best["cols"])
    valid = best_daily.dropna(subset=best["cols"] + ["binge"]).copy()
    valid["latent_control_debt"] = score

    feature_rows = []
    for feature, coef in zip(best["cols"], beta):
        feature_rows.append({"metric": feature, "value": round(coef, 4)})

    summary_rows = [
        {"metric": "analysis_start", "value": START_DATE.strftime("%Y-%m-%d")},
        {"metric": "days", "value": len(best_daily)},
        {"metric": "binges", "value": int(best_daily["binge"].sum())},
        {"metric": "baseline_prev_cal_auc", "value": round(baseline_auc, 4)},
        {"metric": "best_half_life_days", "value": int(best["half_life_days"])},
        {"metric": "latent_auc", "value": round(latent_auc, 4)},
        {"metric": "combo_prev_cal_plus_latent_auc", "value": round(combo_auc, 4)},
        {"metric": "selected_strain_features", "value": ",".join(best["strain_features"])},
        {"metric": "selected_recovery_features", "value": ",".join(best["recovery_features"])},
        {"metric": "mean_latent_before_binge", "value": round(valid.loc[valid["binge"] == 1, "latent_control_debt"].mean(), 4)},
        {"metric": "mean_latent_before_non_binge", "value": round(valid.loc[valid["binge"] == 0, "latent_control_debt"].mean(), 4)},
    ]
    summary_rows.extend(feature_rows)

    daily_cols = ["date", "calories", "binge", "prev_cal", "sleep_hours"] + best["cols"] + ["latent_control_debt"]
    return pd.DataFrame(summary_rows), search_df, valid[daily_cols]


def print_report(summary_df: pd.DataFrame) -> None:
    values = dict(zip(summary_df["metric"], summary_df["value"]))
    print("\n=== Latent Control Capacity Search ===")
    print(
        f"Window: {values['analysis_start']} onward, pre-tirzepatide only | "
        f"days={values['days']}  binges={values['binges']}"
    )
    print(
        f"Baseline prev_cal AUC={values['baseline_prev_cal_auc']} | "
        f"latent AUC={values['latent_auc']} | "
        f"prev_cal + latent AUC={values['combo_prev_cal_plus_latent_auc']}"
    )
    print(
        f"Best half-life={values['best_half_life_days']} days | "
        f"strain={values['selected_strain_features']} | "
        f"recovery={values['selected_recovery_features']}"
    )
    print(
        f"Mean latent debt before binge={values['mean_latent_before_binge']} vs "
        f"non-binge={values['mean_latent_before_non_binge']}"
    )


def save_outputs(summary_df: pd.DataFrame, search_df: pd.DataFrame, daily_df: pd.DataFrame) -> None:
    summary_df.to_csv(ROOT / "analysis" / "Q_latent_control_model_summary.csv", index=False)
    search_df.to_csv(ROOT / "analysis" / "Q_latent_control_feature_search.csv", index=False)
    daily_df.to_csv(ROOT / "analysis" / "Q_latent_control_daily.csv", index=False)


def main() -> None:
    summary_df, search_df, daily_df = build_outputs()
    save_outputs(summary_df, search_df, daily_df)
    print_report(summary_df)


if __name__ == "__main__":
    main()
