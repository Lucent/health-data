"""Test food-pattern antecedents of future intake without latent-state models.

The analysis is intentionally restricted to the modern pre-tirzepatide MHTML
era and uses calendar-ordered discovery, validation, and confirmation periods.
It never reads P2, P3, P4, or modeled drug levels.

Run:

    uv run --with numpy --with pandas python analysis/EB_direct_future_intake.py

Outputs:
    analysis/EB_direct_daily_features.csv
    analysis/EB_direct_future_intake_results.csv
    analysis/EB_direct_event_summaries.csv
    analysis/EB_direct_event_year_summaries.csv
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import math
import re

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent

ANALYSIS_START = pd.Timestamp("2018-03-01")
DISCOVERY_END = pd.Timestamp("2021-12-31")
VALIDATION_END = pd.Timestamp("2022-12-31")
CONFIRMATION_END = pd.Timestamp("2024-09-16")

REPEAT_WINDOW_DAYS = 14
NOVEL_WINDOW_DAYS = 30
TRAILING_BASELINE_DAYS = 28
HIGH_INTAKE_KCAL = 3000

DESSERT_PATTERN = re.compile(
    r"\b(?:ice cream|gelato|cake|cheesecake|cookie|cookies|brownie|brownies|"
    r"candy|chocolate|chocolates|pie|donut|doughnut|scone|cupcake|pudding|"
    r"fudge|truffle|tiramisu)\b",
    re.IGNORECASE,
)
POTATO_PATTERN = re.compile(r"\bpotato(?:es)?\b", re.IGNORECASE)
POTATO_EXCLUSIONS = re.compile(
    r"\b(?:chip|chips|french fries|hash brown|vodka)\b",
    re.IGNORECASE,
)
DOMINOS_PATTERN = re.compile(r"domino", re.IGNORECASE)


def normalize_food_name(value: str) -> str:
    """Reduce portion-specific diary strings to a conservative food key."""
    text = str(value).lower().replace("’", "'").strip()
    text = re.sub(r"\([^)]*\)", " ", text)
    text = re.sub(r",\s*\d+(?:\.\d+)?\s+.*$", "", text)
    text = re.sub(r"\b\d+(?:\.\d+)?\s*(?:g|kg|ml|oz|cup|cups|tbsp|tsp)\b", "", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip(" ,-")


def parse_travel_dates() -> set[pd.Timestamp]:
    """Parse the checked-in trip table; the unavailable raw Timeline is not used."""
    travel_dates: set[pd.Timestamp] = set()
    line_pattern = re.compile(
        r"^\|\s*((?:19|20)\d{2}-\d{2}-\d{2})\s*\|\s*"
        r"((?:19|20)\d{2}-\d{2}-\d{2})\s*\|"
    )
    for line in (ROOT / "travel" / "trips.md").read_text().splitlines():
        match = line_pattern.match(line)
        if not match:
            continue
        start, end = (pd.Timestamp(value) for value in match.groups())
        travel_dates.update(pd.date_range(start, end, freq="D"))
    return travel_dates


def build_food_features() -> pd.DataFrame:
    foods = pd.read_csv(
        ROOT / "intake" / "intake_foods.csv",
        parse_dates=["date"],
        usecols=["date", "food", "calories"],
    )
    foods["food_key"] = foods["food"].map(normalize_food_name)
    foods["calories"] = pd.to_numeric(foods["calories"], errors="coerce").fillna(0)
    foods["is_dessert"] = foods["food"].str.contains(DESSERT_PATTERN, na=False)
    foods["is_potato"] = (
        foods["food"].str.contains(POTATO_PATTERN, na=False)
        & ~foods["food"].str.contains(POTATO_EXCLUSIONS, na=False)
    )
    foods["is_dominos"] = foods["food"].str.contains(DOMINOS_PATTERN, na=False)

    by_key = (
        foods.groupby(["date", "food_key"], as_index=False)
        .agg(
            calories=("calories", "sum"),
            is_dessert=("is_dessert", "max"),
            is_potato=("is_potato", "max"),
            is_dominos=("is_dominos", "max"),
        )
        .sort_values(["date", "food_key"])
    )

    rows = []
    last_seen: dict[str, int] = {}
    for day_index, (date, group) in enumerate(by_key.groupby("date", sort=True)):
        energy = group["calories"].clip(lower=0).to_numpy(dtype=float)
        total = energy.sum()
        shares = energy / total if total > 0 else np.zeros(len(group))
        positive = shares > 0
        entropy = -float(np.sum(shares[positive] * np.log(shares[positive])))
        n_clusters = int(len(group))
        normalized_entropy = (
            entropy / math.log(n_clusters) if n_clusters > 1 else 0.0
        )

        repeat_energy = 0.0
        novel_energy = 0.0
        for key, calories in zip(group["food_key"], energy):
            age = day_index - last_seen[key] if key in last_seen else None
            if age is not None and age <= REPEAT_WINDOW_DAYS:
                repeat_energy += calories
            if age is None or age > NOVEL_WINDOW_DAYS:
                novel_energy += calories

        dessert_energy = float(
            group.loc[group["is_dessert"], "calories"].clip(lower=0).sum()
        )
        potato_energy = float(
            group.loc[group["is_potato"], "calories"].clip(lower=0).sum()
        )

        rows.append(
            {
                "date": date,
                "food_clusters": n_clusters,
                "menu_entropy_normalized": normalized_entropy,
                "largest_food_share": float(shares.max()) if len(shares) else np.nan,
                "repeat_share_14d": repeat_energy / total if total > 0 else np.nan,
                "novel_share_30d": novel_energy / total if total > 0 else np.nan,
                "dessert_share": dessert_energy / total if total > 0 else np.nan,
                "potato_share": potato_energy / total if total > 0 else np.nan,
                "dominos_day": int(group["is_dominos"].any()),
            }
        )

        for key in group["food_key"]:
            last_seen[key] = day_index

    return pd.DataFrame(rows)


def build_daily_table() -> pd.DataFrame:
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    features = build_food_features()
    travel_dates = parse_travel_dates()

    daily = intake.merge(features, on="date", how="left").sort_values("date")
    daily["travel_today"] = daily["date"].isin(travel_dates).astype(int)
    daily["travel_tomorrow"] = (
        daily["date"] + pd.Timedelta(days=1)
    ).isin(travel_dates).astype(int)
    daily["prior_28d_mean"] = (
        daily["calories"]
        .shift(1)
        .rolling(TRAILING_BASELINE_DAYS, min_periods=14)
        .mean()
    )
    daily["target_weekday"] = (daily["date"] + pd.Timedelta(days=1)).dt.dayofweek
    daily["day_of_year"] = daily["date"].dt.dayofyear
    daily["season_sin"] = np.sin(2 * np.pi * daily["day_of_year"] / 365.25)
    daily["season_cos"] = np.cos(2 * np.pi * daily["day_of_year"] / 365.25)
    daily["time_days"] = (daily["date"] - ANALYSIS_START).dt.days

    daily["next_day_calories"] = daily["calories"].shift(-1)
    daily["next_day_high_intake"] = (
        daily["next_day_calories"] >= HIGH_INTAKE_KCAL
    ).astype(float)
    daily.loc[daily["next_day_calories"].isna(), "next_day_high_intake"] = np.nan
    future = pd.concat(
        [daily["calories"].shift(-offset) for offset in (1, 2, 3)],
        axis=1,
    )
    daily["next_3d_mean_calories"] = future.mean(axis=1)
    daily["next_day_residual_vs_prior"] = (
        daily["next_day_calories"] - daily["prior_28d_mean"]
    )

    daily["split"] = "outside"
    daily.loc[
        daily["date"].between(ANALYSIS_START, DISCOVERY_END),
        "split",
    ] = "discovery"
    daily.loc[
        daily["date"].between(DISCOVERY_END + pd.Timedelta(days=1), VALIDATION_END),
        "split",
    ] = "validation"
    daily.loc[
        daily["date"].between(VALIDATION_END + pd.Timedelta(days=1), CONFIRMATION_END),
        "split",
    ] = "confirmation"
    return daily.reset_index(drop=True)


def design_matrix(
    frame: pd.DataFrame,
    train_means: dict[str, float] | None = None,
    train_stds: dict[str, float] | None = None,
) -> tuple[np.ndarray, dict[str, float], dict[str, float]]:
    numeric = ["calories", "prior_28d_mean", "season_sin", "season_cos", "time_days"]
    if train_means is None:
        train_means = {column: float(frame[column].mean()) for column in numeric}
    if train_stds is None:
        train_stds = {
            column: max(float(frame[column].std(ddof=0)), 1e-9)
            for column in numeric
        }

    columns = [np.ones(len(frame))]
    for column in numeric:
        columns.append(
            (frame[column].to_numpy(dtype=float) - train_means[column])
            / train_stds[column]
        )
    for weekday in range(1, 7):
        columns.append((frame["target_weekday"].to_numpy() == weekday).astype(float))
    return np.column_stack(columns), train_means, train_stds


def correlation(left: np.ndarray, right: np.ndarray) -> float:
    if len(left) < 3 or np.std(left) == 0 or np.std(right) == 0:
        return np.nan
    return float(np.corrcoef(left, right)[0, 1])


def partial_outcome_correlations(
    daily: pd.DataFrame,
    feature: str,
    outcome: str,
    sample: str,
) -> list[dict[str, object]]:
    required = [
        feature,
        outcome,
        "calories",
        "prior_28d_mean",
        "season_sin",
        "season_cos",
        "time_days",
        "target_weekday",
    ]
    eligible = daily[
        daily["split"].isin(["discovery", "validation", "confirmation"])
    ].dropna(subset=required)
    if sample == "home_only":
        eligible = eligible[
            (eligible["travel_today"] == 0)
            & (eligible["travel_tomorrow"] == 0)
        ]
    train = eligible[eligible["split"] == "discovery"].copy()
    if len(train) < 100 or train[feature].std(ddof=0) == 0:
        return []

    x_train, means, stds = design_matrix(train)
    feature_mean = float(train[feature].mean())
    feature_std = max(float(train[feature].std(ddof=0)), 1e-9)
    feature_train = (train[feature].to_numpy(dtype=float) - feature_mean) / feature_std
    outcome_train = train[outcome].to_numpy(dtype=float)

    beta_outcome = np.linalg.lstsq(x_train, outcome_train, rcond=None)[0]
    beta_feature = np.linalg.lstsq(x_train, feature_train, rcond=None)[0]

    results = []
    for split in ["discovery", "validation", "confirmation"]:
        frame = eligible[eligible["split"] == split].copy()
        if frame[feature].nunique() < 2:
            results.append(
                {
                    "feature": feature,
                    "outcome": outcome,
                    "sample": sample,
                    "split": split,
                    "n": len(frame),
                    "partial_r": np.nan,
                    "residual_slope_per_discovery_sd": np.nan,
                    "feature_mean": float(frame[feature].mean()),
                    "feature_sd": float(frame[feature].std(ddof=0)),
                }
            )
            continue
        x_split, _, _ = design_matrix(frame, means, stds)
        feature_scaled = (
            frame[feature].to_numpy(dtype=float) - feature_mean
        ) / feature_std
        outcome_values = frame[outcome].to_numpy(dtype=float)
        outcome_residual = outcome_values - x_split @ beta_outcome
        feature_residual = feature_scaled - x_split @ beta_feature
        denominator = float(np.dot(feature_residual, feature_residual))
        slope = (
            float(np.dot(feature_residual, outcome_residual) / denominator)
            if denominator > 0
            else np.nan
        )
        results.append(
            {
                "feature": feature,
                "outcome": outcome,
                "sample": sample,
                "split": split,
                "n": len(frame),
                "partial_r": correlation(feature_residual, outcome_residual),
                "residual_slope_per_discovery_sd": slope,
                "feature_mean": float(frame[feature].mean()),
                "feature_sd": float(frame[feature].std(ddof=0)),
            }
        )
    return results


def event_summaries(daily: pd.DataFrame) -> pd.DataFrame:
    eligible = daily[
        daily["split"].isin(["discovery", "validation", "confirmation"])
    ].copy()
    event_definitions = {
        "travel day": eligible["travel_today"] == 1,
        "Dominos day": eligible["dominos_day"] == 1,
        "potato majority day": eligible["potato_share"] >= 0.50,
        "dessert >=25% calories": eligible["dessert_share"] >= 0.25,
        "high repetition >=90%": eligible["repeat_share_14d"] >= 0.90,
        "high novelty >=50%": eligible["novel_share_30d"] >= 0.50,
    }
    rows = []
    for event, event_mask in event_definitions.items():
        for split in ["discovery", "validation", "confirmation"]:
            split_mask = eligible["split"] == split
            event_frame = eligible[split_mask & event_mask]
            other_frame = eligible[split_mask & ~event_mask]
            rows.append(
                {
                    "event": event,
                    "split": split,
                    "event_days": len(event_frame),
                    "other_days": len(other_frame),
                    "event_day_mean_calories": event_frame["calories"].mean(),
                    "other_day_mean_calories": other_frame["calories"].mean(),
                    "event_next_day_mean_calories": event_frame[
                        "next_day_calories"
                    ].mean(),
                    "other_next_day_mean_calories": other_frame[
                        "next_day_calories"
                    ].mean(),
                    "event_next_day_high_intake_rate": event_frame[
                        "next_day_high_intake"
                    ].mean(),
                    "other_next_day_high_intake_rate": other_frame[
                        "next_day_high_intake"
                    ].mean(),
                }
            )
    return pd.DataFrame(rows)


def event_year_summaries(daily: pd.DataFrame) -> pd.DataFrame:
    eligible = daily[daily["date"].between(ANALYSIS_START, CONFIRMATION_END)].copy()
    event_definitions = {
        "travel day": eligible["travel_today"] == 1,
        "Dominos day": eligible["dominos_day"] == 1,
        "potato majority day": eligible["potato_share"] >= 0.50,
        "dessert >=25% calories": eligible["dessert_share"] >= 0.25,
        "high repetition >=90%": eligible["repeat_share_14d"] >= 0.90,
        "high novelty >=50%": eligible["novel_share_30d"] >= 0.50,
    }
    rows = []
    for event, event_mask in event_definitions.items():
        for year in sorted(eligible["date"].dt.year.unique()):
            year_mask = eligible["date"].dt.year == year
            event_frame = eligible[year_mask & event_mask]
            other_frame = eligible[year_mask & ~event_mask]
            rows.append(
                {
                    "event": event,
                    "year": int(year),
                    "event_days": len(event_frame),
                    "other_days": len(other_frame),
                    "same_day_delta": (
                        event_frame["calories"].mean()
                        - other_frame["calories"].mean()
                    ),
                    "next_day_delta": (
                        event_frame["next_day_calories"].mean()
                        - other_frame["next_day_calories"].mean()
                    ),
                    "event_next_day_high_intake_rate": event_frame[
                        "next_day_high_intake"
                    ].mean(),
                    "other_next_day_high_intake_rate": other_frame[
                        "next_day_high_intake"
                    ].mean(),
                }
            )
    return pd.DataFrame(rows)


def print_report(
    results: pd.DataFrame,
    events: pd.DataFrame,
    event_years: pd.DataFrame,
) -> None:
    print("DIRECT FUTURE-INTAKE ANALYSIS")
    print("=" * 108)
    print(
        f"Discovery: {ANALYSIS_START.date()} to {DISCOVERY_END.date()} | "
        f"Validation: 2022-01-01 to {VALIDATION_END.date()} | "
        f"Confirmation: 2023-01-01 to {CONFIRMATION_END.date()}"
    )
    print("No P2/P3/P4 states and no modeled drug levels are used.")

    for outcome in results["outcome"].unique():
        print(f"\nOutcome: {outcome}")
        subset = results[
            (results["outcome"] == outcome) & (results["sample"] == "all_days")
        ]
        pivot = subset.pivot(
            index="feature",
            columns="split",
            values="partial_r",
        )
        order = pivot.get("confirmation", pd.Series(index=pivot.index)).abs().sort_values(
            ascending=False
        ).index
        print(
            pivot.reindex(order)[
                ["discovery", "validation", "confirmation"]
            ].round(3).to_string()
        )

    print("\nHome-only sensitivity for repetition and novelty")
    home = results[
        (results["sample"] == "home_only")
        & results["feature"].isin(["repeat_share_14d", "novel_share_30d"])
    ]
    for outcome in home["outcome"].unique():
        print(f"\nOutcome: {outcome}")
        print(
            home[home["outcome"] == outcome]
            .pivot(index="feature", columns="split", values="partial_r")[
                ["discovery", "validation", "confirmation"]
            ]
            .round(3)
            .to_string()
        )

    print("\nDirect event summaries")
    display = events.copy()
    display["same_day_delta"] = (
        display["event_day_mean_calories"] - display["other_day_mean_calories"]
    )
    display["next_day_delta"] = (
        display["event_next_day_mean_calories"]
        - display["other_next_day_mean_calories"]
    )
    print(
        display[
            [
                "event",
                "split",
                "event_days",
                "same_day_delta",
                "next_day_delta",
                "event_next_day_high_intake_rate",
            ]
        ].round(3).to_string(index=False)
    )

    print("\nYear-level stress test for the two generic menu signals")
    year_display = event_years[
        event_years["event"].isin(
            ["high repetition >=90%", "high novelty >=50%"]
        )
    ]
    print(
        year_display[
            ["event", "year", "event_days", "same_day_delta", "next_day_delta"]
        ].round(1).to_string(index=False)
    )


def main() -> None:
    daily = build_daily_table()
    feature_columns = [
        "food_clusters",
        "menu_entropy_normalized",
        "largest_food_share",
        "repeat_share_14d",
        "novel_share_30d",
        "dessert_share",
        "potato_share",
        "dominos_day",
        "travel_today",
    ]
    outcome_columns = [
        "next_day_calories",
        "next_3d_mean_calories",
        "next_day_high_intake",
    ]

    result_rows = []
    for sample in ["all_days", "home_only"]:
        for feature in feature_columns:
            for outcome in outcome_columns:
                result_rows.extend(
                    partial_outcome_correlations(
                        daily,
                        feature,
                        outcome,
                        sample,
                    )
                )
    results = pd.DataFrame(result_rows)
    events = event_summaries(daily)
    event_years = event_year_summaries(daily)

    daily.to_csv(ROOT / "analysis" / "EB_direct_daily_features.csv", index=False)
    results.to_csv(
        ROOT / "analysis" / "EB_direct_future_intake_results.csv",
        index=False,
    )
    events.to_csv(
        ROOT / "analysis" / "EB_direct_event_summaries.csv",
        index=False,
    )
    event_years.to_csv(
        ROOT / "analysis" / "EB_direct_event_year_summaries.csv",
        index=False,
    )
    print_report(results, events, event_years)


if __name__ == "__main__":
    main()
