"""Audit long-run weight trajectories using direct observations only.

This script deliberately avoids P1/P2/P3/P4/P4b and any inferred daily body
composition or expenditure. It combines the canonical daily intake ledger with
observed scale weights, checked-in travel dates, canonical strength dates, and
a deterministic food-name repetition feature.

The episode catalogue is descriptive and retrospective. The quarterly warning
test is time ordered: all scaling and regression coefficients are fit through
2020, then applied unchanged to 2021 through the pre-tirzepatide boundary.
Because the repository has already been studied extensively, that later period
is an internal locked test, not a pristine prospective holdout.

Run:

    uv run --with numpy --with pandas python analysis/EL_long_run_trajectory_audit.py

Outputs:
    analysis/EL_long_run_episodes.csv
    analysis/EL_quarterly_warning_anchors.csv
    analysis/EL_locked_warning_results.csv
    analysis/EL_weight_gaps.csv
    analysis/EL_weekly_weight_swings.csv
"""

from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent

INTAKE_START = pd.Timestamp("2011-04-20")
INTAKE_END = pd.Timestamp("2026-02-28")
PRE_TIRZ_END = pd.Timestamp("2024-09-16")
LOCKED_TEST_START = pd.Timestamp("2021-01-01")

HIGH_INTAKE_KCAL = 3000
REPEAT_WINDOW_DAYS = 14
EPISODE_CONTEXT_DAYS = 56
QUARTER_DAYS = 91
QUARTER_ANCHOR_TOLERANCE_DAYS = 21
QUARTER_FEATURE_DAYS = 28
GAIN_THRESHOLD_LBS = 2.0
LONG_WEIGHT_GAP_DAYS = 30
VERY_LONG_WEIGHT_GAP_DAYS = 90
WEEKLY_DENSE_START = pd.Timestamp("2019-10-21")

# These are observed anchors chosen to summarize the conspicuous long-run
# movements in the raw scale record. They are not change points estimated by a
# model and are not independent hypothesis tests.
EPISODES = [
    (
        "willpower-era long loss",
        "loss",
        "2011-04-21",
        "2013-09-23",
        "Long observed loss ending at the last pre-gap weigh-in.",
    ),
    (
        "sparsely observed rebound",
        "regain",
        "2013-09-23",
        "2015-10-17",
        "Crosses a 517-day weigh-in gap and the first annotated travel binge.",
    ),
    (
        "2015-16 renewed loss",
        "loss",
        "2015-10-17",
        "2016-04-08",
        "Dense enough to show a sustained multi-month decline.",
    ),
    (
        "2016-17 regain",
        "regain",
        "2016-04-08",
        "2017-09-11",
        "Observed regain before the high-protein/keto experiments.",
    ),
    (
        "2017-19 sparse regain",
        "regain",
        "2017-10-06",
        "2019-10-23",
        "Only one intervening weight; magnitude is observed, path is not.",
    ),
    (
        "2019-21 broad maintenance",
        "maintenance",
        "2019-10-23",
        "2021-12-21",
        "Dense modern record with fluctuations but little net change.",
    ),
    (
        "early-2022 loss",
        "loss",
        "2022-01-04",
        "2022-06-06",
        "Multi-month loss; endpoint is a local low and single-day magnitude is noisy.",
    ),
    (
        "mid-2022 to 2023 regain",
        "regain",
        "2022-06-06",
        "2023-04-01",
        "Reversal of the early-2022 decline, sustained into spring 2023.",
    ),
    (
        "2023 to pre-treatment drift",
        "regain",
        "2023-04-01",
        "2024-09-16",
        "Slow observed rise before tirzepatide.",
    ),
    (
        "tirzepatide-era loss",
        "loss",
        "2024-09-16",
        "2026-02-28",
        "Observed treatment-era loss through the canonical intake boundary.",
    ),
    (
        "late-treatment maintenance",
        "maintenance",
        "2025-12-02",
        "2026-02-28",
        "Overlapping final interval documenting the low-weight plateau.",
    ),
]


def normalize_food_name(value: str) -> str:
    """Conservatively collapse portion-specific strings to a food key."""
    text = str(value).lower().replace("’", "'").strip()
    text = re.sub(r"\([^)]*\)", " ", text)
    text = re.sub(r",\s*\d+(?:\.\d+)?\s+.*$", "", text)
    text = re.sub(
        r"\b\d+(?:\.\d+)?\s*(?:g|kg|ml|oz|cup|cups|tbsp|tsp)\b",
        "",
        text,
    )
    return re.sub(r"\s+", " ", text).strip(" ,-")


def build_repeat_share() -> pd.DataFrame:
    """Calculate calorie share from foods seen in the prior 14 calendar days."""
    foods = pd.read_csv(
        ROOT / "intake" / "intake_foods.csv",
        parse_dates=["date"],
        usecols=["date", "food", "calories"],
    )
    foods["food_key"] = foods["food"].map(normalize_food_name)
    foods["calories"] = pd.to_numeric(
        foods["calories"], errors="coerce"
    ).fillna(0.0)
    grouped = (
        foods.groupby(["date", "food_key"], as_index=False)["calories"]
        .sum()
        .sort_values(["date", "food_key"])
    )

    last_seen: dict[str, pd.Timestamp] = {}
    rows = []
    for date, day in grouped.groupby("date", sort=True):
        calories = day["calories"].clip(lower=0)
        total = float(calories.sum())
        repeat_calories = 0.0
        for food_key, energy in zip(day["food_key"], calories):
            previous = last_seen.get(food_key)
            if previous is not None and (date - previous).days <= REPEAT_WINDOW_DAYS:
                repeat_calories += float(energy)
        rows.append(
            {
                "date": date,
                "repeat_share_14d": (
                    repeat_calories / total if total > 0 else np.nan
                ),
            }
        )
        for food_key in day["food_key"]:
            last_seen[food_key] = date
    return pd.DataFrame(rows)


def parse_travel_dates() -> set[pd.Timestamp]:
    """Parse checked-in trips from the repository table."""
    travel_dates: set[pd.Timestamp] = set()
    line_pattern = re.compile(
        r"^\|\s*((?:19|20)\d{2}-\d{2}-\d{2})\s*\|\s*"
        r"((?:19|20)\d{2}-\d{2}-\d{2})\s*\|"
    )
    for line in (ROOT / "travel" / "trips.md").read_text().splitlines():
        match = line_pattern.match(line)
        if match is None:
            continue
        start, end = (pd.Timestamp(value) for value in match.groups())
        travel_dates.update(pd.date_range(start, end, freq="D"))
    return travel_dates


def build_daily() -> pd.DataFrame:
    intake = pd.read_csv(
        ROOT / "intake" / "intake_daily.csv", parse_dates=["date"]
    ).sort_values("date")
    repeat = build_repeat_share()
    strength = pd.read_csv(
        ROOT / "workout" / "strength.csv", parse_dates=["date"]
    )
    strength_daily = (
        strength.groupby("date")
        .agg(
            strength_sessions=("date", "size"),
            strength_minutes=("duration_min", "sum"),
        )
        .reset_index()
    )
    steps = pd.read_csv(
        ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"]
    )[["date", "steps"]]
    travel_dates = parse_travel_dates()

    daily = (
        intake.merge(repeat, on="date", how="left")
        .merge(strength_daily, on="date", how="left")
        .merge(steps, on="date", how="left")
        .sort_values("date")
    )
    daily["strength_sessions"] = daily["strength_sessions"].fillna(0).astype(int)
    daily["strength_minutes"] = daily["strength_minutes"].fillna(0.0)
    daily["travel"] = daily["date"].isin(travel_dates).astype(int)
    daily["high_intake"] = (daily["calories"] >= HIGH_INTAKE_KCAL).astype(int)
    return daily


def observed_anchor(weights: pd.DataFrame, date: pd.Timestamp) -> pd.Series:
    row = weights.loc[weights["date"] == date]
    if len(row) != 1:
        raise ValueError(f"Expected one observed weight on {date.date()}, got {len(row)}")
    return row.iloc[0]


def local_weight_summary(
    weights: pd.DataFrame, date: pd.Timestamp, radius_days: int = 7
) -> tuple[float, int]:
    frame = weights[
        weights["date"].between(
            date - pd.Timedelta(days=radius_days),
            date + pd.Timedelta(days=radius_days),
        )
    ]
    return float(frame["weight_lbs"].median()), int(len(frame))


def summarize_window(frame: pd.DataFrame) -> dict[str, float]:
    if frame.empty:
        return {
            "days": 0,
            "mean_calories": np.nan,
            "high_intake_rate": np.nan,
            "travel_days": 0,
            "travel_rate": np.nan,
            "mean_repeat_share": np.nan,
            "strength_sessions_per_week": np.nan,
            "steps_coverage": np.nan,
            "mean_steps_when_observed": np.nan,
        }
    return {
        "days": int(len(frame)),
        "mean_calories": float(frame["calories"].mean()),
        "high_intake_rate": float(frame["high_intake"].mean()),
        "travel_days": int(frame["travel"].sum()),
        "travel_rate": float(frame["travel"].mean()),
        "mean_repeat_share": float(frame["repeat_share_14d"].mean()),
        "strength_sessions_per_week": float(
            frame["strength_sessions"].sum() * 7 / len(frame)
        ),
        "steps_coverage": float(frame["steps"].notna().mean()),
        "mean_steps_when_observed": float(frame["steps"].mean()),
    }


def prefix_values(prefix: str, values: dict[str, float]) -> dict[str, float]:
    return {f"{prefix}_{key}": value for key, value in values.items()}


def build_episode_catalog(
    daily: pd.DataFrame, weights: pd.DataFrame
) -> pd.DataFrame:
    rows = []
    for name, direction, start_text, end_text, note in EPISODES:
        start = pd.Timestamp(start_text)
        end = pd.Timestamp(end_text)
        start_row = observed_anchor(weights, start)
        end_row = observed_anchor(weights, end)
        start_median, start_local_n = local_weight_summary(weights, start)
        end_median, end_local_n = local_weight_summary(weights, end)

        pre = daily[
            daily["date"].between(
                start - pd.Timedelta(days=EPISODE_CONTEXT_DAYS),
                start - pd.Timedelta(days=1),
            )
        ]
        during = daily[daily["date"].between(start, end)]
        post = daily[
            daily["date"].between(
                end + pd.Timedelta(days=1),
                end + pd.Timedelta(days=EPISODE_CONTEXT_DAYS),
            )
        ]
        observed = weights[weights["date"].between(start, end)].copy()
        observed_gaps = observed["date"].diff().dt.days
        duration_days = int((end - start).days)
        exact_change = float(end_row["weight_lbs"] - start_row["weight_lbs"])

        row = {
            "episode": name,
            "direction_label": direction,
            "start": start,
            "end": end,
            "duration_days": duration_days,
            "start_observed_weight_lbs": float(start_row["weight_lbs"]),
            "end_observed_weight_lbs": float(end_row["weight_lbs"]),
            "observed_change_lbs": exact_change,
            "annualized_observed_change_lbs": (
                exact_change * 365.25 / duration_days
            ),
            "start_local_15d_median_lbs": start_median,
            "start_local_weight_n": start_local_n,
            "end_local_15d_median_lbs": end_median,
            "end_local_weight_n": end_local_n,
            "local_median_change_lbs": end_median - start_median,
            "weight_observations_during": int(len(observed)),
            "max_between_weight_gap_days": (
                float(observed_gaps.max()) if len(observed) > 1 else np.nan
            ),
            "selection_note": note,
        }
        row.update(prefix_values("pre56", summarize_window(pre)))
        row.update(prefix_values("during", summarize_window(during)))
        row.update(prefix_values("post56", summarize_window(post)))
        rows.append(row)
    return pd.DataFrame(rows)


def nearest_weight(
    weights: pd.DataFrame, target: pd.Timestamp, tolerance_days: int
) -> pd.Series | None:
    distance = (weights["date"] - target).abs().dt.days
    if distance.empty or int(distance.min()) > tolerance_days:
        return None
    return weights.loc[distance.idxmin()]


def previous_weight_gap_days(
    weights: pd.DataFrame, anchor_date: pd.Timestamp
) -> float:
    prior = weights.loc[weights["date"] < anchor_date, "date"]
    if prior.empty:
        return np.nan
    return float((anchor_date - prior.max()).days)


def build_quarterly_anchors(
    daily: pd.DataFrame, weights: pd.DataFrame
) -> pd.DataFrame:
    rows = []
    quarter_targets = pd.date_range(
        "2011-04-01", "2024-07-01", freq="QS"
    )
    for quarter in quarter_targets:
        end_target = quarter + pd.Timedelta(days=QUARTER_DAYS)
        start_row = nearest_weight(
            weights, quarter, QUARTER_ANCHOR_TOLERANCE_DAYS
        )
        end_row = nearest_weight(
            weights, end_target, QUARTER_ANCHOR_TOLERANCE_DAYS
        )
        if start_row is None or end_row is None:
            continue
        start_date = pd.Timestamp(start_row["date"])
        end_date = pd.Timestamp(end_row["date"])
        if end_date > PRE_TIRZ_END or end_date <= start_date:
            continue

        prior = daily[
            daily["date"].between(
                start_date - pd.Timedelta(days=QUARTER_FEATURE_DAYS),
                start_date - pd.Timedelta(days=1),
            )
        ]
        if len(prior) != QUARTER_FEATURE_DAYS:
            continue
        prior90_weight_count = int(
            weights["date"]
            .between(start_date - pd.Timedelta(days=90), start_date - pd.Timedelta(days=1))
            .sum()
        )
        previous_gap = previous_weight_gap_days(weights, start_date)
        rows.append(
            {
                "quarter_target": quarter,
                "split": (
                    "internal_locked_test"
                    if quarter >= LOCKED_TEST_START
                    else "development"
                ),
                "start_weight_date": start_date,
                "start_anchor_offset_days": int((start_date - quarter).days),
                "start_weight_lbs": float(start_row["weight_lbs"]),
                "end_weight_date": end_date,
                "end_anchor_offset_days": int((end_date - end_target).days),
                "end_weight_lbs": float(end_row["weight_lbs"]),
                "outcome_days": int((end_date - start_date).days),
                "future_weight_change_lbs": float(
                    end_row["weight_lbs"] - start_row["weight_lbs"]
                ),
                "future_gain_ge_2lb": int(
                    end_row["weight_lbs"] - start_row["weight_lbs"]
                    >= GAIN_THRESHOLD_LBS
                ),
                "prior28_mean_calories": float(prior["calories"].mean()),
                "prior28_high_intake_rate": float(prior["high_intake"].mean()),
                "prior28_travel_days": int(prior["travel"].sum()),
                "prior28_mean_repeat_share": float(
                    prior["repeat_share_14d"].mean()
                ),
                "prior28_strength_sessions_per_week": float(
                    prior["strength_sessions"].sum() * 7 / len(prior)
                ),
                "days_since_previous_weight": previous_gap,
                "prior90_weight_observations": prior90_weight_count,
            }
        )
    return pd.DataFrame(rows)


def correlation(left: np.ndarray, right: np.ndarray) -> float:
    if len(left) < 3 or np.std(left) == 0 or np.std(right) == 0:
        return np.nan
    return float(np.corrcoef(left, right)[0, 1])


def rank_auc(score: np.ndarray, outcome: np.ndarray) -> float:
    positives = outcome == 1
    negatives = outcome == 0
    if positives.sum() == 0 or negatives.sum() == 0:
        return np.nan
    ranks = pd.Series(score).rank(method="average").to_numpy()
    rank_sum = float(ranks[positives].sum())
    n_positive = int(positives.sum())
    n_negative = int(negatives.sum())
    return (
        rank_sum - n_positive * (n_positive + 1) / 2
    ) / (n_positive * n_negative)


def prediction_metrics(
    outcome: np.ndarray, prediction: np.ndarray, binary: np.ndarray
) -> dict[str, float]:
    error = prediction - outcome
    return {
        "mae_lbs": float(np.mean(np.abs(error))),
        "rmse_lbs": float(np.sqrt(np.mean(error**2))),
        "prediction_outcome_r": correlation(prediction, outcome),
        "gain_auc": rank_auc(prediction, binary),
    }


def build_locked_warning_results(anchors: pd.DataFrame) -> pd.DataFrame:
    development = anchors[anchors["split"] == "development"].copy()
    locked = anchors[anchors["split"] == "internal_locked_test"].copy()
    outcome_column = "future_weight_change_lbs"
    binary_column = "future_gain_ge_2lb"
    baseline = float(development[outcome_column].mean())

    rows = []
    baseline_prediction = np.repeat(baseline, len(locked))
    baseline_metrics = prediction_metrics(
        locked[outcome_column].to_numpy(dtype=float),
        baseline_prediction,
        locked[binary_column].to_numpy(dtype=int),
    )
    rows.append(
        {
            "warning_feature": "development mean baseline",
            "expected_risk_direction": "none",
            "development_n": len(development),
            "locked_test_n": len(locked),
            "development_mean": np.nan,
            "development_sd": np.nan,
            "development_slope_lbs_per_sd": 0.0,
            "locked_feature_outcome_r": np.nan,
            "locked_high_risk_minus_low_risk_lbs": np.nan,
            "locked_high_risk_n": np.nan,
            "locked_low_risk_n": np.nan,
            "development_matches_expected_direction": np.nan,
            "locked_matches_expected_direction": np.nan,
            "locked_replicates_development_direction": np.nan,
            "mae_improvement_vs_baseline_lbs": 0.0,
            **baseline_metrics,
        }
    )

    feature_specs = [
        ("prior28_mean_calories", 1, "higher"),
        ("prior28_high_intake_rate", 1, "higher"),
        ("prior28_travel_days", 1, "higher"),
        ("prior28_mean_repeat_share", -1, "lower"),
        ("prior28_strength_sessions_per_week", -1, "lower"),
        ("days_since_previous_weight", 1, "higher"),
        ("prior90_weight_observations", -1, "lower"),
    ]

    standardized_development: dict[str, np.ndarray] = {}
    standardized_locked: dict[str, np.ndarray] = {}
    for feature, expected_sign, direction in feature_specs:
        mean = float(development[feature].mean())
        sd = max(float(development[feature].std(ddof=0)), 1e-9)
        z_development = (
            development[feature].to_numpy(dtype=float) - mean
        ) / sd
        z_locked = (locked[feature].to_numpy(dtype=float) - mean) / sd
        standardized_development[feature] = z_development
        standardized_locked[feature] = z_locked

        x = np.column_stack([np.ones(len(development)), z_development])
        beta = np.linalg.lstsq(
            x,
            development[outcome_column].to_numpy(dtype=float),
            rcond=None,
        )[0]
        prediction = beta[0] + beta[1] * z_locked
        locked_correlation = correlation(
            locked[feature].to_numpy(dtype=float),
            locked[outcome_column].to_numpy(dtype=float),
        )
        risk_score = expected_sign * z_locked
        high = risk_score >= 0
        high_minus_low = (
            float(
                locked.loc[high, outcome_column].mean()
                - locked.loc[~high, outcome_column].mean()
            )
            if high.any() and (~high).any()
            else np.nan
        )
        metrics = prediction_metrics(
            locked[outcome_column].to_numpy(dtype=float),
            prediction,
            locked[binary_column].to_numpy(dtype=int),
        )
        rows.append(
            {
                "warning_feature": feature,
                "expected_risk_direction": direction,
                "development_n": len(development),
                "locked_test_n": len(locked),
                "development_mean": mean,
                "development_sd": sd,
                "development_slope_lbs_per_sd": float(beta[1]),
                "locked_feature_outcome_r": locked_correlation,
                "locked_high_risk_minus_low_risk_lbs": high_minus_low,
                "locked_high_risk_n": int(high.sum()),
                "locked_low_risk_n": int((~high).sum()),
                "development_matches_expected_direction": bool(
                    beta[1] * expected_sign > 0
                ),
                "locked_matches_expected_direction": bool(
                    locked_correlation * expected_sign > 0
                ),
                "locked_replicates_development_direction": bool(
                    beta[1] * locked_correlation > 0
                ),
                "mae_improvement_vs_baseline_lbs": (
                    baseline_metrics["mae_lbs"] - metrics["mae_lbs"]
                ),
                **metrics,
            }
        )

    # Equal-weight, direction-prespecified score. It uses no locked outcomes and
    # intentionally avoids fitting feature weights in the small development set.
    composite_features = [
        ("prior28_mean_calories", 1),
        ("prior28_high_intake_rate", 1),
        ("prior28_travel_days", 1),
        ("prior28_mean_repeat_share", -1),
        ("prior28_strength_sessions_per_week", -1),
        ("days_since_previous_weight", 1),
    ]
    composite_development = np.mean(
        [
            sign * standardized_development[feature]
            for feature, sign in composite_features
        ],
        axis=0,
    )
    composite_locked = np.mean(
        [
            sign * standardized_locked[feature]
            for feature, sign in composite_features
        ],
        axis=0,
    )
    x = np.column_stack([np.ones(len(development)), composite_development])
    beta = np.linalg.lstsq(
        x,
        development[outcome_column].to_numpy(dtype=float),
        rcond=None,
    )[0]
    prediction = beta[0] + beta[1] * composite_locked
    locked_correlation = correlation(
        composite_locked,
        locked[outcome_column].to_numpy(dtype=float),
    )
    high = composite_locked >= 0
    metrics = prediction_metrics(
        locked[outcome_column].to_numpy(dtype=float),
        prediction,
        locked[binary_column].to_numpy(dtype=int),
    )
    rows.append(
        {
            "warning_feature": "equal_weight_context_tracking_composite",
            "expected_risk_direction": "higher",
            "development_n": len(development),
            "locked_test_n": len(locked),
            "development_mean": float(composite_development.mean()),
            "development_sd": float(composite_development.std(ddof=0)),
            "development_slope_lbs_per_sd": float(beta[1]),
            "locked_feature_outcome_r": locked_correlation,
            "locked_high_risk_minus_low_risk_lbs": (
                float(
                    locked.loc[high, outcome_column].mean()
                    - locked.loc[~high, outcome_column].mean()
                )
                if high.any() and (~high).any()
                else np.nan
            ),
            "locked_high_risk_n": int(high.sum()),
            "locked_low_risk_n": int((~high).sum()),
            "development_matches_expected_direction": bool(beta[1] > 0),
            "locked_matches_expected_direction": bool(locked_correlation > 0),
            "locked_replicates_development_direction": bool(
                beta[1] * locked_correlation > 0
            ),
            "mae_improvement_vs_baseline_lbs": (
                baseline_metrics["mae_lbs"] - metrics["mae_lbs"]
            ),
            **metrics,
        }
    )
    return pd.DataFrame(rows)


def build_weight_gaps(
    daily: pd.DataFrame, weights: pd.DataFrame
) -> pd.DataFrame:
    in_scope = weights[
        weights["date"].between(INTAKE_START, INTAKE_END)
    ].sort_values("date")
    rows = []
    for previous, current in zip(
        in_scope.iloc[:-1].itertuples(), in_scope.iloc[1:].itertuples()
    ):
        gap_days = int((current.date - previous.date).days)
        if gap_days < LONG_WEIGHT_GAP_DAYS:
            continue
        during = daily[
            daily["date"].between(
                previous.date + pd.Timedelta(days=1),
                current.date - pd.Timedelta(days=1),
            )
        ]
        pre = daily[
            daily["date"].between(
                previous.date - pd.Timedelta(days=EPISODE_CONTEXT_DAYS),
                previous.date - pd.Timedelta(days=1),
            )
        ]
        rows.append(
            {
                "previous_weight_date": previous.date,
                "next_weight_date": current.date,
                "gap_days": gap_days,
                "gap_class": (
                    "90+ days"
                    if gap_days >= VERY_LONG_WEIGHT_GAP_DAYS
                    else "30-89 days"
                ),
                "previous_weight_lbs": float(previous.weight_lbs),
                "next_weight_lbs": float(current.weight_lbs),
                "observed_change_lbs": float(
                    current.weight_lbs - previous.weight_lbs
                ),
                "annualized_change_lbs": float(
                    (current.weight_lbs - previous.weight_lbs)
                    * 365.25
                    / gap_days
                ),
                "gap_mean_calories": float(during["calories"].mean()),
                "pre56_mean_calories": float(pre["calories"].mean()),
                "gap_minus_pre56_calories": float(
                    during["calories"].mean() - pre["calories"].mean()
                ),
                "gap_high_intake_rate": float(during["high_intake"].mean()),
                "gap_travel_days": int(during["travel"].sum()),
            }
        )
    return pd.DataFrame(rows)


def build_weekly_weight_swings(weights: pd.DataFrame) -> pd.DataFrame:
    modern = weights[
        weights["date"].between(WEEKLY_DENSE_START, INTAKE_END)
    ].copy()
    weekly = (
        modern.set_index("date")["weight_lbs"]
        .resample("W-SUN")
        .agg(["median", "count"])
        .reset_index()
        .rename(
            columns={
                "date": "week_end",
                "median": "weekly_median_weight_lbs",
                "count": "weight_observations",
            }
        )
    )
    weekly.loc[
        weekly["weight_observations"] < 2, "weekly_median_weight_lbs"
    ] = np.nan
    weekly["one_week_change_lbs"] = weekly[
        "weekly_median_weight_lbs"
    ].diff()
    weekly["four_week_change_lbs"] = (
        weekly["weekly_median_weight_lbs"]
        - weekly["weekly_median_weight_lbs"].shift(4)
    )
    weekly["next_week_change_lbs"] = weekly[
        "one_week_change_lbs"
    ].shift(-1)
    weekly["large_one_week_swing"] = (
        weekly["one_week_change_lbs"].abs() >= 5.0
    )
    weekly["next_week_opposes_swing"] = (
        weekly["one_week_change_lbs"] * weekly["next_week_change_lbs"] < 0
    )
    return weekly


def print_report(
    episodes: pd.DataFrame,
    anchors: pd.DataFrame,
    results: pd.DataFrame,
    gaps: pd.DataFrame,
    weekly: pd.DataFrame,
) -> None:
    print("LONG-RUN DIRECT TRAJECTORY AUDIT")
    print("=" * 104)
    print(
        "Inputs: canonical intake, observed scale weights, checked-in travel, "
        "canonical strength dates, steps, and deterministic food repetition."
    )
    print("No P1/P2/P3/P4/P4b states or inferred TDEE are read.")

    print("\nEpisode catalogue")
    print(
        episodes[
            [
                "episode",
                "start",
                "end",
                "observed_change_lbs",
                "local_median_change_lbs",
                "during_mean_calories",
                "during_high_intake_rate",
                "during_travel_days",
                "during_mean_repeat_share",
                "during_strength_sessions_per_week",
                "weight_observations_during",
                "max_between_weight_gap_days",
            ]
        ].to_string(index=False)
    )

    development = anchors[anchors["split"] == "development"]
    locked = anchors[anchors["split"] == "internal_locked_test"]
    print(
        f"\nQuarterly internal validation: {len(development)} development anchors "
        f"through 2020; {len(locked)} locked-test anchors in 2021-pre-treatment."
    )
    baseline_mae = float(
        results.loc[
            results["warning_feature"] == "development mean baseline", "mae_lbs"
        ].iloc[0]
    )
    print(f"Locked-test development-mean baseline MAE: {baseline_mae:.2f} lb")
    print(
        results[
            [
                "warning_feature",
                "development_slope_lbs_per_sd",
                "locked_feature_outcome_r",
                "locked_high_risk_minus_low_risk_lbs",
                "mae_lbs",
                "mae_improvement_vs_baseline_lbs",
                "gain_auc",
                "development_matches_expected_direction",
                "locked_matches_expected_direction",
            ]
        ].to_string(index=False)
    )

    for label, frame in [
        ("30+ day gaps", gaps),
        ("90+ day gaps", gaps[gaps["gap_days"] >= VERY_LONG_WEIGHT_GAP_DAYS]),
    ]:
        if frame.empty:
            continue
        print(
            f"\n{label}: n={len(frame)}, "
            f"higher endpoint={(frame['observed_change_lbs'] > 0).sum()}/{len(frame)}, "
            f"median change={frame['observed_change_lbs'].median():+.1f} lb."
        )

    valid_weekly = weekly["one_week_change_lbs"].dropna()
    large = weekly[weekly["large_one_week_swing"]]
    print(
        "\nDense-era weekly medians: "
        f"n={len(valid_weekly)}, median absolute one-week change="
        f"{valid_weekly.abs().median():.2f} lb, 95th percentile="
        f"{valid_weekly.abs().quantile(0.95):.2f} lb."
    )
    print(
        f"One-week swings >=5 lb: {len(large)}; "
        f"opposed the following week: "
        f"{int(large['next_week_opposes_swing'].sum())}/{len(large)}."
    )


def main() -> None:
    daily = build_daily()
    weights = pd.read_csv(
        ROOT / "weight" / "weight.csv", parse_dates=["date"]
    ).sort_values("date")
    weights = weights[
        weights["date"].between(INTAKE_START, INTAKE_END)
    ].copy()

    episodes = build_episode_catalog(daily, weights)
    anchors = build_quarterly_anchors(daily, weights)
    results = build_locked_warning_results(anchors)
    gaps = build_weight_gaps(daily, weights)
    weekly = build_weekly_weight_swings(weights)

    episodes.to_csv(ROOT / "analysis" / "EL_long_run_episodes.csv", index=False)
    anchors.to_csv(
        ROOT / "analysis" / "EL_quarterly_warning_anchors.csv", index=False
    )
    results.to_csv(
        ROOT / "analysis" / "EL_locked_warning_results.csv", index=False
    )
    gaps.to_csv(ROOT / "analysis" / "EL_weight_gaps.csv", index=False)
    weekly.to_csv(
        ROOT / "analysis" / "EL_weekly_weight_swings.csv", index=False
    )
    print_report(episodes, anchors, results, gaps, weekly)


if __name__ == "__main__":
    main()
