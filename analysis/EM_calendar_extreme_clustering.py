"""Audit calendar effects and clustering of directly observed extreme-intake days.

This script uses only the canonical daily intake ledger and checked-in travel
dates. It does not read latent body-composition, TDEE, set-point, or modeled
drug-exposure files.

Run:

    uv run --with numpy --with pandas python analysis/EM_calendar_extreme_clustering.py

Outputs:
    analysis/EM_calendar_dayofweek.csv
    analysis/EM_epoch_extreme_summary.csv
    analysis/EM_extreme_transitions.csv
    analysis/EM_extreme_runs.csv
    analysis/EM_forward_prediction.csv
    analysis/EM_daily_features.csv
"""

from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
INTAKE_PATH = ROOT / "intake" / "intake_daily.csv"
TRAVEL_PATH = ROOT / "travel" / "trips.md"
EPOCH_PATH = ROOT / "intake" / "diet_epochs.csv"

EXTREME_KCAL = 3000
TIRZEPATIDE_START = pd.Timestamp("2024-09-17")

DISCOVERY_END = pd.Timestamp("2017-12-31")
VALIDATION_END = pd.Timestamp("2021-12-31")
CONFIRMATION_END = TIRZEPATIDE_START - pd.Timedelta(days=1)

MIN_HISTORY_DAYS = 14
HAC_LAGS = 7
LOGISTIC_MAX_ITERATIONS = 200
LOGISTIC_TOLERANCE = 1e-10
NUMERICAL_RIDGE = 1e-8

WEEKDAY_NAMES = {
    0: "Monday",
    1: "Tuesday",
    2: "Wednesday",
    3: "Thursday",
    4: "Friday",
    5: "Saturday",
    6: "Sunday",
}


def parse_travel_dates() -> set[pd.Timestamp]:
    """Parse the repository's checked-in trip table."""
    dates: set[pd.Timestamp] = set()
    line_pattern = re.compile(
        r"^\|\s*((?:19|20)\d{2}-\d{2}-\d{2})\s*\|\s*"
        r"((?:19|20)\d{2}-\d{2}-\d{2})\s*\|"
    )
    for line in TRAVEL_PATH.read_text().splitlines():
        match = line_pattern.match(line)
        if not match:
            continue
        start, end = (pd.Timestamp(value) for value in match.groups())
        dates.update(pd.date_range(start, end, freq="D"))
    return dates


def days_since_prior_event(events: pd.Series) -> pd.Series:
    """Return days since the most recent earlier event, capped later."""
    result = np.full(len(events), np.nan)
    most_recent: int | None = None
    values = events.to_numpy(dtype=int)
    for index in range(len(values)):
        if most_recent is not None:
            result[index] = index - most_recent
        if values[index] == 1:
            most_recent = index
    return pd.Series(result, index=events.index)


def annotate_diet_epochs(daily: pd.DataFrame) -> pd.DataFrame:
    """Attach overlapping declared epoch labels without forcing a partition."""
    epochs = pd.read_csv(EPOCH_PATH, parse_dates=["start", "end"])
    labels: list[list[str]] = [[] for _ in range(len(daily))]
    dates = daily["date"]
    for epoch in epochs.itertuples():
        positions = np.flatnonzero(dates.between(epoch.start, epoch.end).to_numpy())
        for position in positions:
            labels[position].append(str(epoch.label))
    daily["diet_epoch_labels"] = [";".join(values) for values in labels]
    daily["diet_epoch_count"] = [len(values) for values in labels]
    return daily


def build_daily() -> pd.DataFrame:
    daily = pd.read_csv(INTAKE_PATH, parse_dates=["date"]).sort_values("date")
    daily = daily.reset_index(drop=True)
    expected_dates = pd.date_range(daily["date"].min(), daily["date"].max(), freq="D")
    if not daily["date"].equals(pd.Series(expected_dates, name="date")):
        raise ValueError("Canonical intake ledger is not one row per consecutive day")

    travel_dates = parse_travel_dates()
    daily["travel_today"] = daily["date"].isin(travel_dates).astype(int)
    daily["weekday"] = daily["date"].dt.dayofweek
    daily["weekday_name"] = daily["weekday"].map(WEEKDAY_NAMES)
    daily["weekend"] = daily["weekday"].isin([5, 6]).astype(int)
    daily["is_extreme"] = (daily["calories"] >= EXTREME_KCAL).astype(int)
    daily["isolated_extreme"] = (
        (daily["is_extreme"] == 1)
        & (daily["is_extreme"].shift(1).fillna(0) == 0)
        & (daily["is_extreme"].shift(-1).fillna(0) == 0)
    ).astype(int)
    daily["clustered_extreme"] = (
        (daily["is_extreme"] == 1) & (daily["isolated_extreme"] == 0)
    ).astype(int)
    daily["medication_era"] = np.where(
        daily["date"] >= TIRZEPATIDE_START,
        "tirzepatide",
        "pre_tirzepatide",
    )

    daily["split"] = "outside"
    daily.loc[daily["date"] <= DISCOVERY_END, "split"] = "discovery"
    daily.loc[
        daily["date"].between(DISCOVERY_END + pd.Timedelta(days=1), VALIDATION_END),
        "split",
    ] = "validation"
    daily.loc[
        daily["date"].between(VALIDATION_END + pd.Timedelta(days=1), CONFIRMATION_END),
        "split",
    ] = "confirmation"
    daily.loc[daily["date"] >= TIRZEPATIDE_START, "split"] = "treatment_shift"

    daily["time_days"] = (daily["date"] - daily["date"].min()).dt.days
    daily["day_of_year"] = daily["date"].dt.dayofyear
    daily["season_sin"] = np.sin(2 * np.pi * daily["day_of_year"] / 365.25)
    daily["season_cos"] = np.cos(2 * np.pi * daily["day_of_year"] / 365.25)

    prior_calories = daily["calories"].shift(1)
    prior_extreme = daily["is_extreme"].shift(1)
    daily["prior_1d_calories"] = prior_calories
    daily["prior_1d_extreme"] = prior_extreme
    daily["prior_7d_mean_calories"] = prior_calories.rolling(
        7, min_periods=7
    ).mean()
    daily["prior_7d_max_calories"] = prior_calories.rolling(
        7, min_periods=7
    ).max()
    daily["prior_7d_extreme_count"] = prior_extreme.rolling(
        7, min_periods=7
    ).sum()
    daily["prior_28d_mean_calories"] = prior_calories.rolling(
        28, min_periods=MIN_HISTORY_DAYS
    ).mean()
    daily["prior_1d_vs_28d"] = (
        daily["prior_1d_calories"] - daily["prior_28d_mean_calories"]
    )
    daily["prior_7d_vs_28d"] = (
        daily["prior_7d_mean_calories"] - daily["prior_28d_mean_calories"]
    )
    daily["days_since_prior_extreme"] = days_since_prior_event(daily["is_extreme"])
    daily["days_since_prior_extreme_capped"] = daily[
        "days_since_prior_extreme"
    ].clip(upper=30)

    daily["next_day_calories"] = daily["calories"].shift(-1)
    daily["next_day_extreme"] = daily["is_extreme"].shift(-1)
    daily["travel_tomorrow"] = daily["travel_today"].shift(-1)
    daily["target_weekday"] = daily["weekday"].shift(-1)
    daily["next_day_minus_prior_28d"] = (
        daily["next_day_calories"] - daily["prior_28d_mean_calories"]
    )
    days_2_to_7 = pd.concat(
        [daily["calories"].shift(-offset) for offset in range(2, 8)],
        axis=1,
    )
    daily["days_2_to_7_mean_calories"] = days_2_to_7.mean(
        axis=1, skipna=False
    )
    daily["days_2_to_7_minus_prior_28d"] = (
        daily["days_2_to_7_mean_calories"]
        - daily["prior_28d_mean_calories"]
    )

    week_start = daily["date"] - pd.to_timedelta(daily["weekday"], unit="D")
    daily["week_start"] = week_start
    daily["weekly_mean_calories"] = daily.groupby("week_start")["calories"].transform(
        "mean"
    )
    daily["weekly_extreme_rate"] = daily.groupby("week_start")[
        "is_extreme"
    ].transform("mean")
    daily["calories_minus_week_mean"] = (
        daily["calories"] - daily["weekly_mean_calories"]
    )
    daily["extreme_minus_week_rate"] = (
        daily["is_extreme"] - daily["weekly_extreme_rate"]
    )
    return annotate_diet_epochs(daily)


def calendar_samples(daily: pd.DataFrame) -> list[tuple[str, pd.DataFrame]]:
    pre = daily[daily["medication_era"] == "pre_tirzepatide"]
    treatment = daily[daily["medication_era"] == "tirzepatide"]
    samples = [
        ("all_days", daily),
        ("all_home", daily[daily["travel_today"] == 0]),
        ("pre_tirzepatide", pre),
        ("pre_tirzepatide_home", pre[pre["travel_today"] == 0]),
        ("tirzepatide", treatment),
        ("tirzepatide_home", treatment[treatment["travel_today"] == 0]),
        ("travel_days", daily[daily["travel_today"] == 1]),
    ]
    for split in [
        "discovery",
        "validation",
        "confirmation",
        "treatment_shift",
    ]:
        frame = daily[daily["split"] == split]
        samples.append((split, frame))
        samples.append((f"{split}_home", frame[frame["travel_today"] == 0]))
    return samples


def summarize_calendar(daily: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for sample, frame in calendar_samples(daily):
        bins: list[tuple[str, pd.DataFrame]] = [
            (WEEKDAY_NAMES[weekday], frame[frame["weekday"] == weekday])
            for weekday in range(7)
        ]
        bins.extend(
            [
                ("weekday_Mon_Fri", frame[frame["weekend"] == 0]),
                ("weekend_Sat_Sun", frame[frame["weekend"] == 1]),
            ]
        )
        for calendar_bin, subset in bins:
            rows.append(
                {
                    "sample": sample,
                    "calendar_bin": calendar_bin,
                    "n_days": len(subset),
                    "mean_calories": subset["calories"].mean(),
                    "median_calories": subset["calories"].median(),
                    "extreme_days": int(subset["is_extreme"].sum()),
                    "extreme_rate": subset["is_extreme"].mean(),
                    "mean_calories_minus_week_mean": subset[
                        "calories_minus_week_mean"
                    ].mean(),
                    "mean_extreme_minus_week_rate": subset[
                        "extreme_minus_week_rate"
                    ].mean(),
                }
            )
    return pd.DataFrame(rows)


def summarize_epochs(daily: pd.DataFrame) -> pd.DataFrame:
    """Summarize each overlapping annotation independently."""
    epochs = pd.read_csv(EPOCH_PATH, parse_dates=["start", "end"])
    rows: list[dict[str, object]] = []
    for epoch in epochs.itertuples():
        frame = daily[daily["date"].between(epoch.start, epoch.end)]
        extreme_days = int(frame["is_extreme"].sum())
        rows.append(
            {
                "epoch_label": epoch.label,
                "start": epoch.start,
                "end": epoch.end,
                "declared_days": (epoch.end - epoch.start).days + 1,
                "observed_days": len(frame),
                "mean_calories": frame["calories"].mean(),
                "extreme_days": extreme_days,
                "extreme_rate": frame["is_extreme"].mean(),
                "clustered_extreme_days": int(frame["clustered_extreme"].sum()),
                "clustered_share_of_extreme_days": (
                    frame["clustered_extreme"].sum() / extreme_days
                    if extreme_days
                    else np.nan
                ),
                "travel_days": int(frame["travel_today"].sum()),
                "detail": epoch.detail,
            }
        )
    return pd.DataFrame(rows)


def newey_west_covariance(
    design: np.ndarray,
    residuals: np.ndarray,
    lag_count: int,
) -> np.ndarray:
    """Newey-West covariance with Bartlett weights and fixed lag count."""
    n_rows, n_columns = design.shape
    bread = np.linalg.pinv(design.T @ design)
    meat = np.zeros((n_columns, n_columns))
    score = design * residuals[:, None]
    meat += score.T @ score
    for lag in range(1, min(lag_count, n_rows - 1) + 1):
        weight = 1 - lag / (lag_count + 1)
        cross = score[lag:].T @ score[:-lag]
        meat += weight * (cross + cross.T)
    correction = n_rows / max(n_rows - n_columns, 1)
    return correction * bread @ meat @ bread


def adjusted_event_effect(
    frame: pd.DataFrame,
    outcome: str,
) -> tuple[float, float, float, float, int]:
    required = [
        outcome,
        "is_extreme",
        "prior_28d_mean_calories",
        "target_weekday",
        "season_sin",
        "season_cos",
        "time_days",
        "travel_today",
        "travel_tomorrow",
    ]
    model = frame.dropna(subset=required).copy()
    if len(model) < 50 or model["is_extreme"].nunique() < 2:
        return np.nan, np.nan, np.nan, np.nan, len(model)

    columns = [
        np.ones(len(model)),
        model["is_extreme"].to_numpy(dtype=float),
        model["prior_28d_mean_calories"].to_numpy(dtype=float) / 1000,
        model["season_sin"].to_numpy(dtype=float),
        model["season_cos"].to_numpy(dtype=float),
        model["time_days"].to_numpy(dtype=float) / 1000,
        model["travel_today"].to_numpy(dtype=float),
        model["travel_tomorrow"].to_numpy(dtype=float),
    ]
    for weekday in range(1, 7):
        columns.append(
            (model["target_weekday"].to_numpy(dtype=float) == weekday).astype(float)
        )
    design = np.column_stack(columns)
    values = model[outcome].to_numpy(dtype=float)
    beta = np.linalg.lstsq(design, values, rcond=None)[0]
    residuals = values - design @ beta
    covariance = newey_west_covariance(design, residuals, HAC_LAGS)
    standard_error = float(np.sqrt(max(covariance[1, 1], 0)))
    coefficient = float(beta[1])
    return (
        coefficient,
        standard_error,
        coefficient - 1.96 * standard_error,
        coefficient + 1.96 * standard_error,
        len(model),
    )


def transition_samples(daily: pd.DataFrame) -> list[tuple[str, pd.DataFrame]]:
    samples: list[tuple[str, pd.DataFrame]] = []
    for split in [
        "discovery",
        "validation",
        "confirmation",
        "treatment_shift",
    ]:
        frame = daily[daily["split"] == split]
        samples.append((split, frame))
        samples.append(
            (f"{split}_home", frame[(frame["travel_today"] == 0) & (frame["travel_tomorrow"] == 0)])
        )
    pre = daily[daily["medication_era"] == "pre_tirzepatide"]
    samples.extend(
        [
            ("pre_tirzepatide", pre),
            (
                "pre_tirzepatide_home",
                pre[(pre["travel_today"] == 0) & (pre["travel_tomorrow"] == 0)],
            ),
            ("all_days", daily),
        ]
    )
    return samples


def summarize_transitions(daily: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for sample, frame in transition_samples(daily):
        eligible = frame.dropna(
            subset=[
                "next_day_calories",
                "next_day_extreme",
                "prior_28d_mean_calories",
            ]
        )
        after_extreme = eligible[eligible["is_extreme"] == 1]
        after_other = eligible[eligible["is_extreme"] == 0]
        risk_extreme = after_extreme["next_day_extreme"].mean()
        risk_other = after_other["next_day_extreme"].mean()
        risk_ratio = (
            risk_extreme / risk_other
            if pd.notna(risk_other) and risk_other > 0
            else np.nan
        )

        calorie_effect = adjusted_event_effect(eligible, "next_day_calories")
        probability_effect = adjusted_event_effect(eligible, "next_day_extreme")

        in_sample_events = frame["is_extreme"].sum()
        isolated = frame[frame["isolated_extreme"] == 1]
        rows.append(
            {
                "sample": sample,
                "n_eligible_anchors": len(eligible),
                "current_extreme_days": len(after_extreme),
                "next_extreme_rate_after_extreme": risk_extreme,
                "next_extreme_rate_after_other": risk_other,
                "next_extreme_risk_difference": risk_extreme - risk_other,
                "next_extreme_risk_ratio": risk_ratio,
                "next_mean_calories_after_extreme": after_extreme[
                    "next_day_calories"
                ].mean(),
                "next_mean_calories_after_other": after_other[
                    "next_day_calories"
                ].mean(),
                "next_minus_prior28_after_extreme": after_extreme[
                    "next_day_minus_prior_28d"
                ].mean(),
                "next_minus_prior28_after_other": after_other[
                    "next_day_minus_prior_28d"
                ].mean(),
                "days2_7_minus_prior28_after_extreme": after_extreme[
                    "days_2_to_7_minus_prior_28d"
                ].mean(),
                "days2_7_minus_prior28_after_other": after_other[
                    "days_2_to_7_minus_prior_28d"
                ].mean(),
                "extreme_days_in_sample": int(in_sample_events),
                "isolated_extreme_days": len(isolated),
                "isolated_extreme_share": (
                    len(isolated) / in_sample_events if in_sample_events else np.nan
                ),
                "adjusted_next_calorie_effect": calorie_effect[0],
                "adjusted_next_calorie_se_hac7": calorie_effect[1],
                "adjusted_next_calorie_ci_low": calorie_effect[2],
                "adjusted_next_calorie_ci_high": calorie_effect[3],
                "adjusted_next_extreme_probability_effect": probability_effect[0],
                "adjusted_next_extreme_probability_se_hac7": probability_effect[1],
                "adjusted_next_extreme_probability_ci_low": probability_effect[2],
                "adjusted_next_extreme_probability_ci_high": probability_effect[3],
            }
        )
    return pd.DataFrame(rows)


def build_extreme_runs(daily: pd.DataFrame) -> pd.DataFrame:
    runs: list[dict[str, object]] = []
    event_dates = daily.loc[daily["is_extreme"] == 1, "date"]
    if event_dates.empty:
        return pd.DataFrame()

    run_start = event_dates.iloc[0]
    prior_date = event_dates.iloc[0]
    for date in event_dates.iloc[1:]:
        if (date - prior_date).days > 1:
            runs.append(summarize_run(daily, run_start, prior_date))
            run_start = date
        prior_date = date
    runs.append(summarize_run(daily, run_start, prior_date))
    return pd.DataFrame(runs)


def summarize_run(
    daily: pd.DataFrame,
    start: pd.Timestamp,
    end: pd.Timestamp,
) -> dict[str, object]:
    frame = daily[daily["date"].between(start, end)]
    overlapping_labels = sorted(
        {
            label
            for value in frame["diet_epoch_labels"]
            for label in value.split(";")
            if label
        }
    )
    return {
        "run_start": start,
        "run_end": end,
        "run_length_days": len(frame),
        "mean_calories": frame["calories"].mean(),
        "maximum_calories": frame["calories"].max(),
        "split_at_start": frame.iloc[0]["split"],
        "medication_era_at_start": frame.iloc[0]["medication_era"],
        "travel_days_in_run": int(frame["travel_today"].sum()),
        "any_travel": int(frame["travel_today"].any()),
        "overlapping_diet_epoch_labels": ";".join(overlapping_labels),
    }


def logistic_fit(design: np.ndarray, outcome: np.ndarray) -> np.ndarray:
    """Fit deterministic near-unpenalized logistic regression by IRLS."""
    beta = np.zeros(design.shape[1])
    penalty = np.eye(design.shape[1]) * NUMERICAL_RIDGE
    penalty[0, 0] = 0
    for _ in range(LOGISTIC_MAX_ITERATIONS):
        linear = np.clip(design @ beta, -30, 30)
        probability = 1 / (1 + np.exp(-linear))
        weights = np.clip(probability * (1 - probability), 1e-8, None)
        working = linear + (outcome - probability) / weights
        lhs = design.T @ (design * weights[:, None]) + penalty
        rhs = design.T @ (weights * working)
        updated = np.linalg.solve(lhs, rhs)
        if np.max(np.abs(updated - beta)) < LOGISTIC_TOLERANCE:
            beta = updated
            break
        beta = updated
    return beta


def prepare_design(
    frame: pd.DataFrame,
    features: list[str],
    training_statistics: dict[str, tuple[float, float]] | None = None,
) -> tuple[np.ndarray, dict[str, tuple[float, float]]]:
    if training_statistics is None:
        training_statistics = {}
        for feature in features:
            mean = float(frame[feature].mean())
            standard_deviation = float(frame[feature].std(ddof=0))
            training_statistics[feature] = (
                mean,
                standard_deviation if standard_deviation > 1e-9 else 1.0,
            )
    columns = [np.ones(len(frame))]
    for feature in features:
        mean, standard_deviation = training_statistics[feature]
        columns.append(
            (frame[feature].to_numpy(dtype=float) - mean) / standard_deviation
        )
    for weekday in range(1, 7):
        columns.append((frame["weekday"].to_numpy() == weekday).astype(float))
    return np.column_stack(columns), training_statistics


def roc_auc(outcome: np.ndarray, prediction: np.ndarray) -> float:
    positives = outcome == 1
    n_positive = int(positives.sum())
    n_negative = len(outcome) - n_positive
    if n_positive == 0 or n_negative == 0:
        return np.nan
    ranks = pd.Series(prediction).rank(method="average").to_numpy()
    return float(
        (ranks[positives].sum() - n_positive * (n_positive + 1) / 2)
        / (n_positive * n_negative)
    )


def average_precision(outcome: np.ndarray, prediction: np.ndarray) -> float:
    positives = int(outcome.sum())
    if positives == 0:
        return np.nan
    order = np.argsort(-prediction, kind="mergesort")
    sorted_outcome = outcome[order]
    precision = np.cumsum(sorted_outcome) / np.arange(1, len(outcome) + 1)
    return float((precision * sorted_outcome).sum() / positives)


def prediction_metrics(
    frame: pd.DataFrame,
    prediction: np.ndarray,
    high_risk_threshold: float,
) -> dict[str, float | int]:
    outcome = frame["is_extreme"].to_numpy(dtype=float)
    clipped = np.clip(prediction, 1e-8, 1 - 1e-8)
    high_risk = prediction >= high_risk_threshold
    return {
        "n_days": len(frame),
        "extreme_days": int(outcome.sum()),
        "observed_extreme_rate": float(outcome.mean()),
        "mean_predicted_probability": float(prediction.mean()),
        "brier_score": float(np.mean((prediction - outcome) ** 2)),
        "log_loss": float(
            -np.mean(outcome * np.log(clipped) + (1 - outcome) * np.log(1 - clipped))
        ),
        "roc_auc": roc_auc(outcome, prediction),
        "average_precision": average_precision(outcome, prediction),
        "discovery_probability_p90": high_risk_threshold,
        "high_risk_days": int(high_risk.sum()),
        "high_risk_extreme_rate": (
            float(outcome[high_risk].mean()) if high_risk.any() else np.nan
        ),
        "other_days_extreme_rate": (
            float(outcome[~high_risk].mean()) if (~high_risk).any() else np.nan
        ),
    }


def forward_prediction_audit(daily: pd.DataFrame) -> pd.DataFrame:
    context_features = [
        "prior_28d_mean_calories",
        "season_sin",
        "season_cos",
        "time_days",
        "travel_today",
    ]
    extreme_history_features = context_features + [
        "prior_1d_extreme",
        "prior_7d_extreme_count",
        "days_since_prior_extreme_capped",
    ]
    calorie_history_features = context_features + [
        "prior_1d_vs_28d",
        "prior_7d_vs_28d",
    ]
    history_features = list(
        dict.fromkeys(extreme_history_features + calorie_history_features)
    )
    required = list(dict.fromkeys(history_features + ["weekday", "is_extreme"]))
    eligible = daily[
        daily["split"].isin(["discovery", "validation", "confirmation"])
    ].dropna(subset=required)

    rows: list[dict[str, object]] = []
    for sample in ["all_days", "home_only"]:
        sample_frame = eligible.copy()
        if sample == "home_only":
            sample_frame = sample_frame[sample_frame["travel_today"] == 0]
        training = sample_frame[sample_frame["split"] == "discovery"]
        for model_name, features in [
            ("context_only", context_features),
            ("context_plus_extreme_history", extreme_history_features),
            ("context_plus_calorie_history", calorie_history_features),
            ("context_plus_recent_history", history_features),
        ]:
            train_design, statistics = prepare_design(training, features)
            beta = logistic_fit(
                train_design,
                training["is_extreme"].to_numpy(dtype=float),
            )
            train_probability = 1 / (
                1 + np.exp(-np.clip(train_design @ beta, -30, 30))
            )
            high_risk_threshold = float(np.quantile(train_probability, 0.90))

            for split in ["discovery", "validation", "confirmation"]:
                frame = sample_frame[sample_frame["split"] == split]
                design, _ = prepare_design(frame, features, statistics)
                probability = 1 / (
                    1 + np.exp(-np.clip(design @ beta, -30, 30))
                )
                row = {
                    "sample": sample,
                    "model": model_name,
                    "split": split,
                    **prediction_metrics(frame, probability, high_risk_threshold),
                }
                rows.append(row)
    return pd.DataFrame(rows)


def print_report(
    daily: pd.DataFrame,
    calendar: pd.DataFrame,
    epochs: pd.DataFrame,
    transitions: pd.DataFrame,
    runs: pd.DataFrame,
    predictions: pd.DataFrame,
) -> None:
    print("CALENDAR AND EXTREME-INTAKE CLUSTERING AUDIT")
    print("=" * 110)
    print(
        f"Ledger: {len(daily):,} consecutive days; "
        f"{int(daily['is_extreme'].sum())} days at or above {EXTREME_KCAL:,} kcal."
    )

    print("\nWeekend contrasts")
    weekend = calendar[
        calendar["calendar_bin"].isin(["weekday_Mon_Fri", "weekend_Sat_Sun"])
    ]
    for sample in [
        "pre_tirzepatide",
        "pre_tirzepatide_home",
        "tirzepatide",
        "tirzepatide_home",
    ]:
        view = weekend[weekend["sample"] == sample].set_index("calendar_bin")
        weekday = view.loc["weekday_Mon_Fri"]
        weekend_row = view.loc["weekend_Sat_Sun"]
        print(
            f"{sample:26s} "
            f"weekend-weekday={weekend_row['mean_calories'] - weekday['mean_calories']:+.0f} kcal; "
            f"extreme={weekend_row['extreme_rate']:.1%} vs {weekday['extreme_rate']:.1%}"
        )

    print("\nExact next-day transition")
    for sample in [
        "discovery_home",
        "validation_home",
        "confirmation_home",
        "treatment_shift_home",
    ]:
        row = transitions[transitions["sample"] == sample].iloc[0]
        print(
            f"{sample:26s} "
            f"P(extreme tomorrow|extreme)={row['next_extreme_rate_after_extreme']:.1%}; "
            f"P(extreme tomorrow|other)={row['next_extreme_rate_after_other']:.1%}; "
            f"adjusted kcal={row['adjusted_next_calorie_effect']:+.0f} "
            f"[{row['adjusted_next_calorie_ci_low']:+.0f}, {row['adjusted_next_calorie_ci_high']:+.0f}]"
        )

    multi_day_runs = runs[runs["run_length_days"] >= 2]
    print(
        f"\nExtreme runs: {len(runs)} episodes; "
        f"{len(multi_day_runs)} multi-day episodes; "
        f"longest={int(runs['run_length_days'].max())} days."
    )
    high_rate_epochs = epochs[epochs["observed_days"] >= 10].nlargest(
        5, "extreme_rate"
    )
    print("\nHighest declared-epoch extreme rates (10+ observed days)")
    print(
        high_rate_epochs[
            ["epoch_label", "observed_days", "extreme_days", "extreme_rate"]
        ].to_string(index=False)
    )

    print("\nLocked one-day-ahead prediction")
    metric_columns = [
        "sample",
        "model",
        "split",
        "n_days",
        "extreme_days",
        "brier_score",
        "roc_auc",
        "average_precision",
        "high_risk_extreme_rate",
        "other_days_extreme_rate",
    ]
    print(predictions[metric_columns].to_string(index=False))


def main() -> None:
    daily = build_daily()
    calendar = summarize_calendar(daily)
    epochs = summarize_epochs(daily)
    transitions = summarize_transitions(daily)
    runs = build_extreme_runs(daily)
    predictions = forward_prediction_audit(daily)

    calendar.to_csv(ROOT / "analysis" / "EM_calendar_dayofweek.csv", index=False)
    epochs.to_csv(ROOT / "analysis" / "EM_epoch_extreme_summary.csv", index=False)
    transitions.to_csv(ROOT / "analysis" / "EM_extreme_transitions.csv", index=False)
    runs.to_csv(ROOT / "analysis" / "EM_extreme_runs.csv", index=False)
    predictions.to_csv(ROOT / "analysis" / "EM_forward_prediction.csv", index=False)
    daily.to_csv(ROOT / "analysis" / "EM_daily_features.csv", index=False)
    print_report(daily, calendar, epochs, transitions, runs, predictions)


if __name__ == "__main__":
    main()
