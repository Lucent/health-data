#!/usr/bin/env python3
"""Re-audit sleep and activity with exact calendar joins and direct outcomes.

This analysis repairs three problems in the legacy sleep/activity scripts:

1. ``shift(-1)`` on a sleep-anchored table means "next recorded sleep row",
   not necessarily the next calendar day.
2. ``shift(-1).rolling(horizon)`` is mostly a trailing window, not a forward
   outcome window.
3. Early Samsung sleep coverage includes short daytime episodes that look like
   naps.  The primary analysis therefore uses the nearly complete 2022+ regime
   and a broad 2-14 hour main-sleep plausibility range.

Only direct intake, Samsung steps/sleep/exercise, and observed weights are used.
No P2/P3/P4 latent states enter this analysis.  The script is deterministic and
takes no command-line arguments.

Outputs:
  analysis/EI_sleep_activity_quality.csv
  analysis/EI_sleep_activity_results.csv
  analysis/EI_sleep_activity_holdout.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
START_DATE = pd.Timestamp("2022-01-01")
END_DATE = pd.Timestamp("2026-02-28")
TIRZEPATIDE_START = pd.Timestamp("2024-09-17")
HOLDOUT_START = pd.Timestamp("2025-01-01")
MIN_MAIN_SLEEP_HOURS = 2.0
MAX_MAIN_SLEEP_HOURS = 14.0
BASELINE_DAYS = 28
MIN_BASELINE_DAYS = 14

QUALITY_PATH = ROOT / "analysis" / "EI_sleep_activity_quality.csv"
RESULTS_PATH = ROOT / "analysis" / "EI_sleep_activity_results.csv"
HOLDOUT_PATH = ROOT / "analysis" / "EI_sleep_activity_holdout.csv"


def night_clock(series: pd.Series) -> pd.Series:
    """Map post-midnight bedtimes onto a monotonic night clock."""
    hour = pd.to_timedelta(series + ":00").dt.total_seconds() / 3600
    return pd.Series(np.where(hour < 12, hour + 24, hour), index=series.index)


def exact_future_mean(series: pd.Series, horizon: int) -> pd.Series:
    """Mean of dates t+1 through t+horizon on a complete calendar index."""
    shifted = [series.shift(-lag) for lag in range(1, horizon + 1)]
    return pd.concat(shifted, axis=1).mean(axis=1)


def load_daily() -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    intake = pd.read_csv(
        ROOT / "intake" / "intake_daily.csv", parse_dates=["date"]
    ).set_index("date")
    sleep = pd.read_csv(
        ROOT / "steps-sleep" / "sleep.csv", parse_dates=["date"]
    ).set_index("date")
    steps = pd.read_csv(
        ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"]
    ).set_index("date")
    samsung_steps = pd.read_csv(
        ROOT / "steps-sleep" / "steps_samsung.csv", parse_dates=["date"]
    ).set_index("date")
    weight = (
        pd.read_csv(ROOT / "weight" / "weight.csv", parse_dates=["date"])
        .sort_values(["date", "time"])
        .drop_duplicates("date", keep="first")
        .set_index("date")
    )
    exercise = pd.read_csv(
        ROOT / "steps-sleep" / "exercises_samsung.csv", parse_dates=["date"]
    )
    exercise["duration_min"] = pd.to_numeric(
        exercise["duration_min"], errors="coerce"
    )
    exercise["count"] = pd.to_numeric(exercise["count"], errors="coerce")
    exercise_daily = exercise.groupby("date").agg(
        exercise_sessions=("datauuid", "size"),
        exercise_minutes=("duration_min", "sum"),
        walk_sessions=("exercise_label", lambda values: (values == "walking").sum()),
    )

    index = pd.date_range(START_DATE, END_DATE, name="date")
    daily = (
        pd.DataFrame(index=index)
        .join(intake[["calories"]])
        .join(sleep[["sleep_hours", "sleep_start", "sleep_end"]])
        .join(steps[["steps"]])
        .join(weight[["weight_lbs"]])
        .join(exercise_daily)
    )
    event_columns = ["exercise_sessions", "exercise_minutes", "walk_sessions"]
    daily[event_columns] = daily[event_columns].fillna(0)

    daily["bed_hour"] = night_clock(daily["sleep_start"])
    daily["wake_hour"] = (
        pd.to_timedelta(daily["sleep_end"] + ":00").dt.total_seconds() / 3600
    )
    plausible = daily["sleep_hours"].between(
        MIN_MAIN_SLEEP_HOURS, MAX_MAIN_SLEEP_HOURS
    )
    daily.loc[
        ~plausible, ["sleep_hours", "bed_hour", "wake_hour"]
    ] = np.nan

    for column in ["calories", "sleep_hours", "bed_hour", "wake_hour", "steps"]:
        prior = daily[column].shift(1).rolling(
            BASELINE_DAYS, min_periods=MIN_BASELINE_DAYS
        ).mean()
        daily[f"{column}_prior28"] = prior
        daily[f"{column}_dev"] = daily[column] - prior

    for outcome in ["calories", "steps"]:
        for horizon in [1, 3, 7, 14]:
            future = exact_future_mean(daily[outcome], horizon)
            daily[f"{outcome}_future{horizon}"] = future
            daily[f"{outcome}_future{horizon}_dev"] = (
                future - daily[f"{outcome}_prior28"]
            )

    for horizon in [1, 7, 14, 28]:
        daily[f"weight_delta{horizon}"] = (
            daily["weight_lbs"].shift(-horizon) - daily["weight_lbs"]
        )

    daily["steps_dev_k"] = daily["steps_dev"] / 1000
    daily["steps_7d"] = daily["steps"].rolling(7, min_periods=7).mean()
    daily["steps_preceding_28d"] = daily["steps"].shift(7).rolling(
        28, min_periods=21
    ).mean()
    daily["steps_7d_dev_k"] = (
        daily["steps_7d"] - daily["steps_preceding_28d"]
    ) / 1000

    daily["tirzepatide_era"] = (
        daily.index >= TIRZEPATIDE_START
    ).astype(float)
    daily["time_years"] = (
        (daily.index - daily.index.min()).days / 365.25
    )
    daily["season_sin"] = np.sin(
        2 * np.pi * daily.index.dayofyear / 365.25
    )
    daily["season_cos"] = np.cos(
        2 * np.pi * daily.index.dayofyear / 365.25
    )
    for day in range(1, 7):
        daily[f"dow_{day}"] = (daily.index.dayofweek == day).astype(float)

    sources = {
        "sleep": sleep,
        "steps": steps,
        "samsung_steps": samsung_steps,
        "weight": weight,
        "exercise": exercise,
    }
    return daily, sources


def partial_estimate(
    data: pd.DataFrame,
    outcome: str,
    exposure: str,
    controls: list[str],
) -> tuple[int, float, float]:
    """Return complete-case n, original-unit coefficient, and partial r."""
    frame = data[[outcome, exposure] + controls].replace(
        [np.inf, -np.inf], np.nan
    ).dropna()
    if len(frame) < 30:
        return len(frame), np.nan, np.nan

    control_matrix = np.column_stack(
        [np.ones(len(frame)), frame[controls].to_numpy(float)]
    )
    exposure_values = frame[exposure].to_numpy(float)
    outcome_values = frame[outcome].to_numpy(float)
    design = np.column_stack([control_matrix, exposure_values])
    coefficient = np.linalg.lstsq(
        design, outcome_values, rcond=None
    )[0][-1]

    exposure_residual = exposure_values - control_matrix @ np.linalg.lstsq(
        control_matrix, exposure_values, rcond=None
    )[0]
    outcome_residual = outcome_values - control_matrix @ np.linalg.lstsq(
        control_matrix, outcome_values, rcond=None
    )[0]
    partial_r = np.corrcoef(exposure_residual, outcome_residual)[0, 1]
    return len(frame), float(coefficient), float(partial_r)


def quality_table(
    daily: pd.DataFrame, sources: dict[str, pd.DataFrame]
) -> pd.DataFrame:
    sleep = sources["sleep"]
    post_sleep = sleep.loc[
        (sleep.index >= START_DATE) & (sleep.index <= END_DATE)
    ]
    post_steps = sources["steps"].loc[
        (sources["steps"].index >= START_DATE)
        & (sources["steps"].index <= END_DATE)
    ]
    samsung_dates = set(sources["samsung_steps"].index)
    exercise_steps = (
        sources["exercise"].groupby("date")["count"].sum().rename("exercise_steps")
    )
    step_consistency = post_steps[["steps"]].join(exercise_steps, how="inner")
    step_consistency = step_consistency[step_consistency["exercise_steps"] > 0]
    rows = [
        {
            "metric": "sleep_all_rows",
            "value": len(sleep),
            "denominator": "",
            "pct": "",
            "note": "one longest Samsung sleep episode per wake date",
        },
        {
            "metric": "sleep_all_nap_like_under_2h",
            "value": int((sleep["sleep_hours"] < MIN_MAIN_SLEEP_HOURS).sum()),
            "denominator": len(sleep),
            "pct": 100 * (sleep["sleep_hours"] < MIN_MAIN_SLEEP_HOURS).mean(),
            "note": "all occur before 2022; legacy AE and BS did not remove them",
        },
        {
            "metric": "sleep_post2022_calendar_coverage",
            "value": len(post_sleep),
            "denominator": len(daily),
            "pct": 100 * len(post_sleep) / len(daily),
            "note": "near-complete, stable Samsung regime",
        },
        {
            "metric": "sleep_post2022_plausible_2_to_14h",
            "value": int(daily["sleep_hours"].notna().sum()),
            "denominator": len(daily),
            "pct": 100 * daily["sleep_hours"].notna().mean(),
            "note": "primary analysis set",
        },
        {
            "metric": "sleep_post2022_under_5h",
            "value": int((post_sleep["sleep_hours"] < 5).sum()),
            "denominator": len(post_sleep),
            "pct": 100 * (post_sleep["sleep_hours"] < 5).mean(),
            "note": "only six observations; extreme-sleep inference is underpowered",
        },
        {
            "metric": "steps_post2022_calendar_coverage",
            "value": len(post_steps),
            "denominator": len(daily),
            "pct": 100 * len(post_steps) / len(daily),
            "note": "complete canonical daily series",
        },
        {
            "metric": "steps_post2022_direct_samsung_source",
            "value": int(sum(date in samsung_dates for date in post_steps.index)),
            "denominator": len(post_steps),
            "pct": 100
            * sum(date in samsung_dates for date in post_steps.index)
            / len(post_steps),
            "note": "no MFP/calendar step estimates in primary regime",
        },
        {
            "metric": "steps_post2022_under_100",
            "value": int((post_steps["steps"] < 100).sum()),
            "denominator": len(post_steps),
            "pct": 100 * (post_steps["steps"] < 100).mean(),
            "note": "likely non-wear or partial-wear days; retain an explicit wear-quality flag",
        },
        {
            "metric": "exercise_days_steps_exceed_daily_total",
            "value": int(
                (
                    step_consistency["exercise_steps"]
                    > step_consistency["steps"] + 10
                ).sum()
            ),
            "denominator": len(step_consistency),
            "pct": 100
            * (
                step_consistency["exercise_steps"]
                > step_consistency["steps"] + 10
            ).mean(),
            "note": "device reconciliation failures; do not assume every daily total is complete",
        },
        {
            "metric": "weight_post2022_observed_dates",
            "value": int(daily["weight_lbs"].notna().sum()),
            "denominator": len(daily),
            "pct": 100 * daily["weight_lbs"].notna().mean(),
            "note": "exact observed-date pairs only; no interpolation",
        },
        {
            "metric": "exercise_post2022_sessions",
            "value": int(
                (
                    (sources["exercise"]["date"] >= START_DATE)
                    & (sources["exercise"]["date"] <= END_DATE)
                ).sum()
            ),
            "denominator": "",
            "pct": "",
            "note": "Samsung session records",
        },
    ]
    return pd.DataFrame(rows)


def result_table(daily: pd.DataFrame) -> pd.DataFrame:
    calendar_controls = [
        "tirzepatide_era",
        "time_years",
        "season_sin",
        "season_cos",
    ] + [f"dow_{day}" for day in range(1, 7)]
    periods = {
        "post2022": daily,
        "pre_tirzepatide": daily[daily.index < TIRZEPATIDE_START],
        "tirzepatide_era": daily[daily.index >= TIRZEPATIDE_START],
        "recent_2025_plus": daily[daily.index >= HOLDOUT_START],
    }
    rows: list[dict[str, object]] = []

    def add_result(
        period_name: str,
        data: pd.DataFrame,
        family: str,
        outcome: str,
        exposure: str,
        controls: list[str],
        exposure_unit: str,
        outcome_unit: str,
    ) -> None:
        n_rows, coefficient, partial_r = partial_estimate(
            data, outcome, exposure, controls
        )
        rows.append(
            {
                "period": period_name,
                "family": family,
                "exposure": exposure,
                "exposure_unit": exposure_unit,
                "outcome": outcome,
                "outcome_unit": outcome_unit,
                "n": n_rows,
                "coefficient": coefficient,
                "partial_r": partial_r,
                "controls": ",".join(controls),
            }
        )

    for period_name, period in periods.items():
        for horizon in [1, 3, 7, 14]:
            add_result(
                period_name,
                period,
                "sleep_to_future_intake",
                f"calories_future{horizon}_dev",
                "sleep_hours_dev",
                calendar_controls + ["bed_hour_dev"],
                "hour above prior-28d mean",
                "kcal/day above prior-28d mean",
            )
        for horizon in [1, 3, 7]:
            add_result(
                period_name,
                period,
                "sleep_to_future_steps",
                f"steps_future{horizon}_dev",
                "sleep_hours_dev",
                calendar_controls + ["bed_hour_dev"],
                "hour above prior-28d mean",
                "steps/day above prior-28d mean",
            )
        for horizon in [1, 7]:
            add_result(
                period_name,
                period,
                "bedtime_to_future_intake",
                f"calories_future{horizon}_dev",
                "bed_hour_dev",
                calendar_controls + ["sleep_hours_dev"],
                "hour later than prior-28d mean",
                "kcal/day above prior-28d mean",
            )

        activity_controls = calendar_controls + [
            "sleep_hours_dev",
            "bed_hour_dev",
            "calories_dev",
        ]
        for exposure in ["steps_dev_k", "steps_7d_dev_k"]:
            for horizon in [1, 3, 7, 14]:
                add_result(
                    period_name,
                    period,
                    "steps_to_future_intake",
                    f"calories_future{horizon}_dev",
                    exposure,
                    activity_controls,
                    "1,000 steps above local baseline",
                    "kcal/day above prior-28d mean",
                )

        for exposure, unit in [
            ("exercise_sessions", "session"),
            ("exercise_minutes", "minute"),
            ("walk_sessions", "walk session"),
        ]:
            add_result(
                period_name,
                period,
                "exercise_to_next_day_intake",
                "calories_future1_dev",
                exposure,
                activity_controls + ["steps_dev_k"],
                unit,
                "kcal/day above prior-28d mean",
            )

    weight_controls = calendar_controls + ["weight_lbs"]
    for horizon in [1, 7, 14, 28]:
        add_result(
            "post2022",
            daily,
            "sleep_to_observed_weight",
            f"weight_delta{horizon}",
            "sleep_hours_dev",
            weight_controls + ["bed_hour_dev"],
            "hour above prior-28d mean",
            f"observed lb change over {horizon}d",
        )
        add_result(
            "post2022",
            daily,
            "steps_to_observed_weight",
            f"weight_delta{horizon}",
            "steps_7d_dev_k",
            weight_controls + ["sleep_hours_dev", "bed_hour_dev"],
            "1,000 steps/day above preceding baseline",
            f"observed lb change over {horizon}d",
        )

    return pd.DataFrame(rows)


def holdout_table(daily: pd.DataFrame) -> pd.DataFrame:
    """Locked temporal check: fit through 2024, score on 2025+ dates."""
    base_features = [
        "calories_dev",
        "tirzepatide_era",
        "time_years",
        "season_sin",
        "season_cos",
    ] + [f"dow_{day}" for day in range(1, 7)]
    feature_sets = {
        "baseline": base_features,
        "baseline_plus_sleep": base_features
        + ["sleep_hours_dev", "bed_hour_dev"],
        "baseline_plus_steps": base_features + ["steps_dev_k"],
        "baseline_plus_steps7": base_features + ["steps_7d_dev_k"],
        "full": base_features
        + [
            "sleep_hours_dev",
            "bed_hour_dev",
            "steps_dev_k",
            "steps_7d_dev_k",
        ],
    }
    all_features = list(
        dict.fromkeys(
            feature for features in feature_sets.values() for feature in features
        )
    )
    rows = []
    for horizon in [1, 7]:
        outcome = f"calories_future{horizon}_dev"
        complete = daily[[outcome] + all_features].dropna()
        train = complete[complete.index < HOLDOUT_START]
        test = complete[complete.index >= HOLDOUT_START]
        for model_name, features in feature_sets.items():
            train_matrix = np.column_stack(
                [np.ones(len(train)), train[features].to_numpy(float)]
            )
            coefficients = np.linalg.lstsq(
                train_matrix, train[outcome].to_numpy(float), rcond=None
            )[0]
            test_matrix = np.column_stack(
                [np.ones(len(test)), test[features].to_numpy(float)]
            )
            errors = (
                test_matrix @ coefficients - test[outcome].to_numpy(float)
            )
            rows.append(
                {
                    "outcome": outcome,
                    "model": model_name,
                    "train_start": train.index.min().date().isoformat(),
                    "train_end": train.index.max().date().isoformat(),
                    "test_start": test.index.min().date().isoformat(),
                    "test_end": test.index.max().date().isoformat(),
                    "train_n": len(train),
                    "test_n": len(test),
                    "mae_kcal": float(np.mean(np.abs(errors))),
                    "rmse_kcal": float(np.sqrt(np.mean(errors**2))),
                    "bias_kcal": float(np.mean(errors)),
                    "features": ",".join(features),
                }
            )
    return pd.DataFrame(rows)


def print_summary(
    daily: pd.DataFrame,
    quality: pd.DataFrame,
    results: pd.DataFrame,
    holdout: pd.DataFrame,
) -> None:
    print("=== EI. Sleep and activity re-audit ===")
    print(
        f"Primary calendar: {START_DATE.date()} to {END_DATE.date()} "
        f"({len(daily)} days)"
    )
    print(
        f"Valid sleep: {daily['sleep_hours'].notna().sum()} days; "
        f"mean {daily['sleep_hours'].mean():.2f} h; "
        f"median {daily['sleep_hours'].median():.2f} h"
    )
    print(
        f"Mean bedtime: {daily['bed_hour'].mean() % 24:.2f} clock hour; "
        f"wake: {daily['wake_hour'].mean():.2f}"
    )

    selected = results[
        (results["period"] == "post2022")
        & (
            (
                results["family"] == "sleep_to_future_intake"
            )
            | (results["family"] == "steps_to_future_intake")
        )
        & results["outcome"].isin(
            ["calories_future1_dev", "calories_future7_dev"]
        )
    ]
    print("\nSelected direct associations:")
    print(
        selected[
            [
                "family",
                "exposure",
                "outcome",
                "n",
                "coefficient",
                "partial_r",
            ]
        ].to_string(index=False, float_format=lambda value: f"{value:.3f}")
    )

    print("\nLocked 2025+ temporal holdout:")
    print(
        holdout[
            ["outcome", "model", "test_n", "mae_kcal", "rmse_kcal"]
        ].to_string(index=False, float_format=lambda value: f"{value:.1f}")
    )
    print(f"\nWrote {QUALITY_PATH}")
    print(f"Wrote {RESULTS_PATH}")
    print(f"Wrote {HOLDOUT_PATH}")


def main() -> None:
    daily, sources = load_daily()
    quality = quality_table(daily, sources)
    results = result_table(daily)
    holdout = holdout_table(daily)
    quality.to_csv(QUALITY_PATH, index=False)
    results.to_csv(RESULTS_PATH, index=False)
    holdout.to_csv(HOLDOUT_PATH, index=False)
    print_summary(daily, quality, results, holdout)


if __name__ == "__main__":
    main()
