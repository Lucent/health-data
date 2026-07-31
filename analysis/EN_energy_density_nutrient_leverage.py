#!/usr/bin/env python3
"""Audit food-mass coverage, nutrient density, liquids, and meal structure.

This is an observation-first analysis. It reads only the canonical food/day
ledgers, the checked-in MFP enrichment table, and the checked-in travel table.
It does not read P1/P2/P3/P4, modeled TDEE, modeled body composition, or drug
PK variables.

The analysis deliberately does not infer gram weights for generic servings.
Whole-diet energy density is therefore unavailable on most days. Exact-mass
energy density is reported only as the energy density of the directly weighed
solid-food subset, accompanied by its calorie coverage.

Run:

    uv run --with numpy --with pandas python analysis/EN_energy_density_nutrient_leverage.py

Outputs:
    analysis/EN_data_coverage_by_year_source.csv
    analysis/EN_daily_features.csv
    analysis/EN_feature_split_summary.csv
    analysis/EN_same_day_markers.csv
    analysis/EN_locked_prediction_results.csv
    analysis/EN_current_meal_summary.csv
    analysis/EN_liquid_food_summary.csv
"""

from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent

MODERN_START = pd.Timestamp("2018-03-01")
DISCOVERY_END = pd.Timestamp("2021-12-31")
VALIDATION_END = pd.Timestamp("2022-12-31")
CONFIRMATION_END = pd.Timestamp("2024-09-16")
TREATMENT_START = pd.Timestamp("2024-09-17")
TREATMENT_TRAIN_END = pd.Timestamp("2025-08-18")
DATA_END = pd.Timestamp("2026-02-28")
CURRENT_WINDOW_DAYS = 90

MIN_DAY_KCAL = 500.0
MIN_EXACT_MASS_COVERAGE = 0.25
TRAILING_BASELINE_DAYS = 28
HIGH_INTAKE_KCAL = 3000.0
HAC_LAGS = 7

DIRECT_GRAM_UNITS = re.compile(
    r"^(?:g|gram|grams|gram\(s\))\b",
    re.IGNORECASE,
)
DIRECT_OUNCE_UNITS = re.compile(
    r"^(?:oz|oz\.|oz\(s\)|ounce|ounces)\b",
    re.IGNORECASE,
)
DIRECT_KILOGRAM_UNITS = re.compile(
    r"^(?:kg|kilogram|kilograms)\b",
    re.IGNORECASE,
)
LIQUID_UNITS = re.compile(
    r"^(?:fl\s*\.?\s*oz|fluid\s+ounce|ml|milliliter|milliliters|"
    r"l|liter|liters|litre|litres)\b",
    re.IGNORECASE,
)
LIQUID_NAMES = re.compile(
    r"\b(?:beer|lager|wine|cocktail|cider|juice|lemonade|soda|cola|"
    r"kombucha|smoothie|milkshake|protein shake|energy drink|sports drink|"
    r"latte|cappuccino|espresso|coffee|tea|hot chocolate|milk|water)\b",
    re.IGNORECASE,
)
LIQUID_NAME_EXCLUSIONS = re.compile(
    r"\b(?:cake|ice\s+cream|gelato|bar|candy|cheese|dressing|dip|"
    r"syrup|sauce|soup|tomato|muffin|pretzel|cookie|brownie|"
    r"milk\s+duds?|milk\s+chocolate|milk\s+bread|milk\s+bun|"
    r"peanut\s+butter\s+cup|squares?|"
    r"water\s+cracker|water\s+chest|coconut\s+milk)\b",
    re.IGNORECASE,
)
VALID_MEALS = {"Breakfast", "Lunch", "Dinner", "Snacks", "Supper"}

FEATURES = [
    "protein_g_per_1000_kcal",
    "fiber_g_per_1000_kcal",
    "liquid_calorie_share",
    "snack_calorie_share",
    "food_rows",
    "food_clusters",
    "meals_logged",
    "exact_mass_solid_coverage",
    "exact_mass_solid_energy_density",
]
OUTCOMES = ["next_day_calories", "next_3d_mean_calories"]


def normalize_food_name(value: str) -> str:
    """Reduce portion variants to a conservative food identity."""
    text = str(value).lower().replace("’", "'").strip()
    text = re.sub(r"\([^)]*\)", " ", text)
    text = re.sub(r",\s*\d+(?:\.\d+)?\s+.*$", "", text)
    text = re.sub(
        r"\b\d+(?:\.\d+)?\s*(?:g|kg|ml|oz|cup|cups|tbsp|tsp)\b",
        "",
        text,
    )
    text = re.sub(r"\s+", " ", text)
    return text.strip(" ,-")


def parse_travel_dates() -> set[pd.Timestamp]:
    """Parse the auditable checked-in trip table."""
    travel_dates: set[pd.Timestamp] = set()
    pattern = re.compile(
        r"^\|\s*((?:19|20)\d{2}-\d{2}-\d{2})\s*\|\s*"
        r"((?:19|20)\d{2}-\d{2}-\d{2})\s*\|"
    )
    for line in (ROOT / "travel" / "trips.md").read_text().splitlines():
        match = pattern.match(line)
        if not match:
            continue
        start, end = (pd.Timestamp(value) for value in match.groups())
        travel_dates.update(pd.date_range(start, end, freq="D"))
    return travel_dates


def consumed_grams(row: pd.Series) -> tuple[float, str]:
    """Recover mass only when the logged serving unit itself is a mass."""
    quantity = pd.to_numeric(row.get("diary_qty"), errors="coerce")
    unit = str(row.get("diary_unit_string") or "").strip()
    if pd.isna(quantity) or quantity <= 0 or not unit:
        return np.nan, "unavailable"
    if DIRECT_GRAM_UNITS.match(unit):
        return float(quantity), "direct_gram_unit"
    if DIRECT_OUNCE_UNITS.match(unit):
        return float(quantity) * 28.349523125, "direct_ounce_unit"
    if DIRECT_KILOGRAM_UNITS.match(unit):
        return float(quantity) * 1000.0, "direct_kilogram_unit"
    return np.nan, "unavailable"


def is_liquid(row: pd.Series) -> bool:
    """Conservatively identify caloric beverages by unit or diary name."""
    unit = str(row.get("diary_unit_string") or "").strip()
    name = str(row.get("food") or "")
    by_unit = bool(LIQUID_UNITS.match(unit))
    excluded = bool(LIQUID_NAME_EXCLUSIONS.search(name))
    by_name = bool(LIQUID_NAMES.search(name))
    return (by_unit or by_name) and not excluded


def classify_food_rows() -> pd.DataFrame:
    """Attach parsed units and quality flags to every canonical food row."""
    foods = pd.read_csv(
        ROOT / "intake" / "intake_foods.csv",
        parse_dates=["date"],
        low_memory=False,
    )
    food_dict = pd.read_csv(
        ROOT / "intake" / "mfp_food_dict.csv",
        low_memory=False,
    )
    data = foods.merge(
        food_dict[
            [
                "diary_food_name",
                "status",
                "diary_qty",
                "diary_unit_string",
                "match_tier",
            ]
        ],
        left_on="food",
        right_on="diary_food_name",
        how="left",
        validate="many_to_one",
    )
    for column in [
        "calories",
        "protein_g",
        "fiber_g",
        "carbs_g",
        "fat_g",
        "diary_qty",
        "match_tier",
    ]:
        data[column] = pd.to_numeric(data[column], errors="coerce")

    mass = data.apply(consumed_grams, axis=1, result_type="expand")
    data["consumed_grams"] = pd.to_numeric(mass[0], errors="coerce")
    data["mass_source"] = mass[1]
    data["is_liquid"] = data.apply(is_liquid, axis=1)
    data["is_solid"] = ~data["is_liquid"]
    data["exact_mass_solid"] = (
        data["is_solid"]
        & data["consumed_grams"].gt(0)
        & data["calories"].ge(0)
    )
    data["food_key"] = data["food"].map(normalize_food_name)
    data["meal_valid"] = data["meal"].isin(VALID_MEALS)
    data["enrichment_resolved"] = data["status"].str.startswith(
        "resolved",
        na=False,
    )

    # OXPS renders an unknown nutrient as zero. Positive OXPS values remain
    # informative, but zero cannot be distinguished from a true zero.
    data["protein_recorded"] = data["protein_g"].notna()
    data["fiber_recorded"] = data["fiber_g"].notna()
    data["protein_zero_ambiguous"] = (
        data["source"].eq("oxps") & data["protein_g"].eq(0)
    )
    data["fiber_zero_ambiguous"] = (
        data["source"].eq("oxps") & data["fiber_g"].eq(0)
    )
    data["protein_unambiguous"] = (
        data["protein_recorded"] & ~data["protein_zero_ambiguous"]
    )
    data["fiber_unambiguous"] = (
        data["fiber_recorded"] & ~data["fiber_zero_ambiguous"]
    )
    data["calendar_year"] = data["date"].dt.year
    return data


def safe_share(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if denominator > 0 else np.nan


def summarize_coverage(food_rows: pd.DataFrame) -> pd.DataFrame:
    """Summarize source-dependent field availability and exact-mass coverage."""
    rows: list[dict[str, object]] = []
    groupings: list[tuple[str, object, pd.DataFrame]] = [
        ("all_years", "all", food_rows)
    ]
    groupings.extend(
        ("source", source, frame)
        for source, frame in food_rows.groupby("source", sort=True)
    )
    groupings.extend(
        ("year_source", f"{year}_{source}", frame)
        for (year, source), frame in food_rows.groupby(
            ["calendar_year", "source"],
            sort=True,
        )
    )
    for level, label, frame in groupings:
        calories = float(frame["calories"].clip(lower=0).sum())
        solid_calories = float(
            frame.loc[frame["is_solid"], "calories"].clip(lower=0).sum()
        )
        rows.append(
            {
                "level": level,
                "label": label,
                "year": (
                    int(frame["calendar_year"].iloc[0])
                    if level == "year_source"
                    else np.nan
                ),
                "source": (
                    str(frame["source"].iloc[0])
                    if level in {"source", "year_source"}
                    else "all"
                ),
                "days": frame["date"].nunique(),
                "food_rows": len(frame),
                "calories": calories,
                "calorie_value_row_coverage": frame["calories"].notna().mean(),
                "protein_recorded_row_coverage": frame[
                    "protein_recorded"
                ].mean(),
                "protein_unambiguous_row_lower_bound": frame[
                    "protein_unambiguous"
                ].mean(),
                "protein_zero_ambiguous_rows": int(
                    frame["protein_zero_ambiguous"].sum()
                ),
                "fiber_recorded_row_coverage": frame["fiber_recorded"].mean(),
                "fiber_unambiguous_row_lower_bound": frame[
                    "fiber_unambiguous"
                ].mean(),
                "fiber_zero_ambiguous_rows": int(
                    frame["fiber_zero_ambiguous"].sum()
                ),
                "valid_meal_label_row_coverage": frame["meal_valid"].mean(),
                "enrichment_resolved_row_coverage": frame[
                    "enrichment_resolved"
                ].mean(),
                "exact_mass_row_coverage": frame[
                    "consumed_grams"
                ].notna().mean(),
                "exact_mass_solid_calorie_coverage": safe_share(
                    frame.loc[
                        frame["exact_mass_solid"],
                        "calories",
                    ].clip(lower=0).sum(),
                    solid_calories,
                ),
                "liquid_calorie_share_heuristic": safe_share(
                    frame.loc[
                        frame["is_liquid"],
                        "calories",
                    ].clip(lower=0).sum(),
                    calories,
                ),
            }
        )
    return pd.DataFrame(rows)


def aggregate_daily(food_rows: pd.DataFrame) -> pd.DataFrame:
    """Build observed daily features with explicit data-quality columns."""
    records: list[dict[str, object]] = []
    for date, frame in food_rows.groupby("date", sort=True):
        total_kcal = float(frame["calories"].clip(lower=0).sum())
        solid_kcal = float(
            frame.loc[frame["is_solid"], "calories"].clip(lower=0).sum()
        )
        exact = frame[frame["exact_mass_solid"]]
        exact_kcal = float(exact["calories"].clip(lower=0).sum())
        exact_grams = float(exact["consumed_grams"].sum())
        protein = float(frame["protein_g"].fillna(0).sum())
        fiber = float(frame["fiber_g"].fillna(0).sum())
        source_calories = (
            frame.groupby("source")["calories"].sum().clip(lower=0).to_dict()
        )
        meal_calories = (
            frame.groupby("meal")["calories"].sum().clip(lower=0).to_dict()
        )
        records.append(
            {
                "date": date,
                "food_row_calories": total_kcal,
                "food_rows": len(frame),
                "food_clusters": frame["food_key"].nunique(),
                "meals_logged": frame.loc[
                    frame["calories"].gt(0),
                    "meal",
                ].nunique(),
                "meal_labels_all_valid": bool(frame["meal_valid"].all()),
                "protein_g_from_food_rows": protein,
                "fiber_g_from_food_rows": fiber,
                "protein_g_per_1000_kcal": (
                    1000.0 * protein / total_kcal
                    if total_kcal > 0
                    else np.nan
                ),
                "fiber_g_per_1000_kcal": (
                    1000.0 * fiber / total_kcal
                    if total_kcal > 0
                    else np.nan
                ),
                "protein_unambiguous_calorie_coverage": safe_share(
                    frame.loc[
                        frame["protein_unambiguous"],
                        "calories",
                    ].clip(lower=0).sum(),
                    total_kcal,
                ),
                "fiber_unambiguous_calorie_coverage": safe_share(
                    frame.loc[
                        frame["fiber_unambiguous"],
                        "calories",
                    ].clip(lower=0).sum(),
                    total_kcal,
                ),
                "liquid_calories": float(
                    frame.loc[
                        frame["is_liquid"],
                        "calories",
                    ].clip(lower=0).sum()
                ),
                "liquid_calorie_share": safe_share(
                    frame.loc[
                        frame["is_liquid"],
                        "calories",
                    ].clip(lower=0).sum(),
                    total_kcal,
                ),
                "snack_calories": float(meal_calories.get("Snacks", 0.0)),
                "snack_calorie_share": safe_share(
                    meal_calories.get("Snacks", 0.0),
                    total_kcal,
                ),
                "exact_mass_solid_calories": exact_kcal,
                "exact_mass_solid_grams": exact_grams,
                "exact_mass_solid_coverage": safe_share(
                    exact_kcal,
                    solid_kcal,
                ),
                "exact_mass_solid_energy_density": (
                    exact_kcal / exact_grams
                    if exact_grams > 0
                    else np.nan
                ),
                "source_html_calories": float(source_calories.get("html", 0.0)),
                "source_mhtml_calories": float(
                    source_calories.get("mhtml", 0.0)
                ),
                "source_oxps_calories": float(source_calories.get("oxps", 0.0)),
            }
        )
    return pd.DataFrame(records)


def build_daily_table(food_rows: pd.DataFrame) -> pd.DataFrame:
    intake = pd.read_csv(
        ROOT / "intake" / "intake_daily.csv",
        parse_dates=["date"],
    )
    features = aggregate_daily(food_rows)
    daily = (
        intake.merge(features, on="date", how="left", validate="one_to_one")
        .sort_values("date")
        .reset_index(drop=True)
    )
    daily["calorie_reconciliation_error_kcal"] = (
        daily["food_row_calories"] - daily["calories"]
    )
    daily["protein_reconciliation_error_g"] = (
        daily["protein_g_from_food_rows"] - daily["protein_g"]
    )
    daily["fiber_reconciliation_error_g"] = (
        daily["fiber_g_from_food_rows"] - daily["fiber_g"]
    )
    travel_dates = parse_travel_dates()
    daily["travel_today"] = daily["date"].isin(travel_dates).astype(int)
    for offset in range(1, 4):
        daily[f"travel_t_plus_{offset}"] = (
            daily["date"] + pd.Timedelta(days=offset)
        ).isin(travel_dates).astype(int)
    daily["future_3d_travel_share"] = daily[
        [f"travel_t_plus_{offset}" for offset in range(1, 4)]
    ].mean(axis=1)
    daily["target_weekday"] = (
        daily["date"] + pd.Timedelta(days=1)
    ).dt.dayofweek
    daily["current_weekday"] = daily["date"].dt.dayofweek
    daily["day_of_year"] = daily["date"].dt.dayofyear
    daily["season_sin"] = np.sin(
        2.0 * np.pi * daily["day_of_year"] / 365.25
    )
    daily["season_cos"] = np.cos(
        2.0 * np.pi * daily["day_of_year"] / 365.25
    )
    daily["time_days"] = (daily["date"] - MODERN_START).dt.days
    daily["prior_28d_mean"] = (
        daily["calories"]
        .shift(1)
        .rolling(TRAILING_BASELINE_DAYS, min_periods=14)
        .mean()
    )
    daily["next_day_calories"] = daily["calories"].shift(-1)
    future = pd.concat(
        [daily["calories"].shift(-offset) for offset in range(1, 4)],
        axis=1,
    )
    daily["next_3d_mean_calories"] = future.mean(axis=1, skipna=False)
    daily["next_day_high_intake"] = (
        daily["next_day_calories"] >= HIGH_INTAKE_KCAL
    ).astype(float)
    daily.loc[
        daily["next_day_calories"].isna(),
        "next_day_high_intake",
    ] = np.nan

    daily["split"] = "outside"
    daily.loc[
        daily["date"].between(MODERN_START, DISCOVERY_END),
        "split",
    ] = "discovery"
    daily.loc[
        daily["date"].between(
            DISCOVERY_END + pd.Timedelta(days=1),
            VALIDATION_END,
        ),
        "split",
    ] = "validation"
    daily.loc[
        daily["date"].between(
            VALIDATION_END + pd.Timedelta(days=1),
            CONFIRMATION_END,
        ),
        "split",
    ] = "confirmation"
    daily.loc[
        daily["date"].between(TREATMENT_START, TREATMENT_TRAIN_END),
        "split",
    ] = "treatment_train"
    daily.loc[
        daily["date"].between(
            TREATMENT_TRAIN_END + pd.Timedelta(days=1),
            DATA_END,
        ),
        "split",
    ] = "treatment_holdout"
    return daily


def summarize_feature_splits(daily: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    split_order = [
        "discovery",
        "validation",
        "confirmation",
        "treatment_train",
        "treatment_holdout",
    ]
    for split in split_order:
        for context, context_mask in [
            ("all_context", pd.Series(True, index=daily.index)),
            ("home_only", daily["travel_today"].eq(0)),
            ("travel_only", daily["travel_today"].eq(1)),
        ]:
            frame = daily[
                daily["split"].eq(split)
                & context_mask
                & daily["calories"].ge(MIN_DAY_KCAL)
            ]
            if frame.empty:
                continue
            rows.append(
                {
                    "split": split,
                    "context": context,
                    "start": frame["date"].min().date(),
                    "end": frame["date"].max().date(),
                    "days": len(frame),
                    "mean_calories": frame["calories"].mean(),
                    "mean_protein_g": frame["protein_g"].mean(),
                    "mean_fiber_g": frame["fiber_g"].mean(),
                    "mean_protein_g_per_1000_kcal": frame[
                        "protein_g_per_1000_kcal"
                    ].mean(),
                    "mean_fiber_g_per_1000_kcal": frame[
                        "fiber_g_per_1000_kcal"
                    ].mean(),
                    "mean_liquid_calorie_share": frame[
                        "liquid_calorie_share"
                    ].mean(),
                    "mean_snack_calorie_share": frame[
                        "snack_calorie_share"
                    ].mean(),
                    "mean_food_rows": frame["food_rows"].mean(),
                    "mean_food_clusters": frame["food_clusters"].mean(),
                    "mean_meals_logged": frame["meals_logged"].mean(),
                    "calorie_weighted_exact_mass_coverage": safe_share(
                        frame["exact_mass_solid_calories"].sum(),
                        (
                            frame["calories"] - frame["liquid_calories"]
                        ).clip(lower=0).sum(),
                    ),
                    "days_at_least_25pct_exact_mass_coverage": int(
                        frame["exact_mass_solid_coverage"]
                        .ge(MIN_EXACT_MASS_COVERAGE)
                        .sum()
                    ),
                    "mean_exact_subset_energy_density_on_eligible_days": frame.loc[
                        frame["exact_mass_solid_coverage"].ge(
                            MIN_EXACT_MASS_COVERAGE
                        ),
                        "exact_mass_solid_energy_density",
                    ].mean(),
                }
            )
    return pd.DataFrame(rows)


def add_column(
    columns: list[np.ndarray],
    names: list[str],
    name: str,
    values: pd.Series | np.ndarray,
) -> None:
    array = np.asarray(values, dtype=float)
    columns.append(array)
    names.append(name)


def control_columns(outcome: str) -> list[str]:
    columns = [
        "prior_28d_mean",
        "time_days",
        "season_sin",
        "season_cos",
        "travel_today",
    ]
    if outcome in OUTCOMES:
        columns.append("calories")
    if outcome == "next_day_calories":
        columns.append("travel_t_plus_1")
    if outcome == "next_3d_mean_calories":
        columns.append("future_3d_travel_share")
    return columns


def design_matrix(
    frame: pd.DataFrame,
    outcome: str,
    feature: str | None,
    means: dict[str, float] | None = None,
    stds: dict[str, float] | None = None,
) -> tuple[np.ndarray, list[str], dict[str, float], dict[str, float]]:
    """Build fixed controls; scaling can be locked to a training interval."""
    numeric = control_columns(outcome)
    if means is None:
        means = {
            column: float(frame[column].mean())
            for column in numeric + ([feature] if feature else [])
        }
    if stds is None:
        stds = {
            column: max(float(frame[column].std(ddof=0)), 1e-9)
            for column in numeric + ([feature] if feature else [])
        }

    columns: list[np.ndarray] = []
    names: list[str] = []
    add_column(columns, names, "intercept", np.ones(len(frame)))
    for column in numeric:
        add_column(
            columns,
            names,
            column,
            (frame[column] - means[column]) / stds[column],
        )
    weekday_column = (
        "current_weekday" if outcome == "calories" else "target_weekday"
    )
    for weekday in range(1, 7):
        add_column(
            columns,
            names,
            f"weekday_{weekday}",
            frame[weekday_column].eq(weekday).astype(float),
        )
    if feature:
        add_column(
            columns,
            names,
            "feature_per_training_sd",
            (frame[feature] - means[feature]) / stds[feature],
        )
    return np.column_stack(columns), names, means, stds


def newey_west_covariance(
    design: np.ndarray,
    residuals: np.ndarray,
    lags: int,
) -> np.ndarray:
    """Deterministic Newey-West covariance for date-ordered rows."""
    n_rows, n_columns = design.shape
    xtx_inverse = np.linalg.pinv(design.T @ design)
    meat = np.zeros((n_columns, n_columns), dtype=float)
    for index in range(n_rows):
        vector = design[index] * residuals[index]
        meat += np.outer(vector, vector)
    max_lag = min(lags, n_rows - 1)
    for lag in range(1, max_lag + 1):
        weight = 1.0 - lag / (max_lag + 1.0)
        cross = np.zeros_like(meat)
        for index in range(lag, n_rows):
            current = design[index] * residuals[index]
            previous = design[index - lag] * residuals[index - lag]
            cross += np.outer(current, previous)
        meat += weight * (cross + cross.T)
    scale = n_rows / max(n_rows - n_columns, 1)
    return scale * xtx_inverse @ meat @ xtx_inverse


def feature_eligibility(daily: pd.DataFrame, feature: str) -> pd.Series:
    eligible = daily["calories"].ge(MIN_DAY_KCAL) & daily[feature].notna()
    if feature == "exact_mass_solid_energy_density":
        eligible &= daily["exact_mass_solid_coverage"].ge(
            MIN_EXACT_MASS_COVERAGE
        )
    return eligible


def same_day_markers(daily: pd.DataFrame) -> pd.DataFrame:
    """Estimate descriptive same-day associations separately in ordered eras."""
    rows: list[dict[str, object]] = []
    discovery = daily[
        daily["split"].eq("discovery") & daily["calories"].ge(MIN_DAY_KCAL)
    ]
    discovery_sds = {
        feature: float(discovery.loc[
            feature_eligibility(discovery, feature),
            feature,
        ].std(ddof=0))
        for feature in FEATURES
    }
    for sample in ["all_context", "home_only"]:
        for split in [
            "discovery",
            "validation",
            "confirmation",
            "treatment_train",
            "treatment_holdout",
        ]:
            for feature in FEATURES:
                required = (
                    ["calories", feature, "current_weekday"]
                    + control_columns("calories")
                )
                frame = daily[
                    daily["split"].eq(split)
                    & feature_eligibility(daily, feature)
                ].dropna(subset=required)
                if sample == "home_only":
                    frame = frame[frame["travel_today"].eq(0)]
                frame = frame.sort_values("date")
                scale = discovery_sds[feature]
                if (
                    len(frame) < 50
                    or frame[feature].nunique() < 3
                    or not np.isfinite(scale)
                    or scale <= 0
                ):
                    continue
                means = {
                    column: float(frame[column].mean())
                    for column in control_columns("calories")
                }
                means[feature] = float(frame[feature].mean())
                stds = {
                    column: max(float(frame[column].std(ddof=0)), 1e-9)
                    for column in control_columns("calories")
                }
                stds[feature] = scale
                design, names, _, _ = design_matrix(
                    frame,
                    "calories",
                    feature,
                    means,
                    stds,
                )
                outcome = frame["calories"].to_numpy(dtype=float)
                beta = np.linalg.lstsq(design, outcome, rcond=None)[0]
                residuals = outcome - design @ beta
                covariance = newey_west_covariance(
                    design,
                    residuals,
                    HAC_LAGS,
                )
                index = names.index("feature_per_training_sd")
                effect = float(beta[index])
                standard_error = float(
                    np.sqrt(max(covariance[index, index], 0.0))
                )
                rows.append(
                    {
                        "sample": sample,
                        "split": split,
                        "feature": feature,
                        "n": len(frame),
                        "effect_kcal_per_discovery_sd": effect,
                        "hac_se": standard_error,
                        "ci95_low": effect - 1.96 * standard_error,
                        "ci95_high": effect + 1.96 * standard_error,
                        "discovery_sd_in_native_units": scale,
                        "feature_mean": frame[feature].mean(),
                        "feature_sd": frame[feature].std(ddof=0),
                        "note": (
                            "same-day descriptive marker; exposure and calorie "
                            "outcome share rows and may be mathematically coupled"
                        ),
                    }
                )
    return pd.DataFrame(rows)


def target_within_split(frame: pd.DataFrame, outcome: str) -> pd.Series:
    offset = 1 if outcome == "next_day_calories" else 3
    target_date = frame["date"] + pd.Timedelta(days=offset)
    if frame.empty:
        return pd.Series(False, index=frame.index)
    split = str(frame["split"].iloc[0])
    boundaries = {
        "discovery": (MODERN_START, DISCOVERY_END),
        "validation": (
            DISCOVERY_END + pd.Timedelta(days=1),
            VALIDATION_END,
        ),
        "confirmation": (
            VALIDATION_END + pd.Timedelta(days=1),
            CONFIRMATION_END,
        ),
        "treatment_train": (TREATMENT_START, TREATMENT_TRAIN_END),
        "treatment_holdout": (
            TREATMENT_TRAIN_END + pd.Timedelta(days=1),
            DATA_END,
        ),
    }
    start, end = boundaries[split]
    return target_date.between(start, end)


def prediction_results(daily: pd.DataFrame) -> pd.DataFrame:
    """Score discovery-locked and treatment-era-locked one-feature models."""
    rows: list[dict[str, object]] = []
    regimes = [
        (
            "modern_pre_treatment",
            "discovery",
            ["validation", "confirmation"],
        ),
        (
            "treatment_transport",
            "treatment_train",
            ["treatment_holdout"],
        ),
    ]
    for regime, training_split, holdouts in regimes:
        for sample in ["all_context", "home_only"]:
            for feature in FEATURES:
                for outcome in OUTCOMES:
                    required = (
                        [feature, outcome, "target_weekday"]
                        + control_columns(outcome)
                    )
                    train = daily[
                        daily["split"].eq(training_split)
                        & feature_eligibility(daily, feature)
                    ].dropna(subset=required)
                    train = train[target_within_split(train, outcome)]
                    if sample == "home_only":
                        travel_future = (
                            train["travel_t_plus_1"].eq(0)
                            if outcome == "next_day_calories"
                            else train["future_3d_travel_share"].eq(0)
                        )
                        train = train[
                            train["travel_today"].eq(0) & travel_future
                        ]
                    train = train.sort_values("date")
                    if len(train) < 100 or train[feature].nunique() < 3:
                        continue

                    baseline_design, _, baseline_means, baseline_stds = (
                        design_matrix(train, outcome, None)
                    )
                    extended_design, extended_names, extended_means, extended_stds = (
                        design_matrix(train, outcome, feature)
                    )
                    y_train = train[outcome].to_numpy(dtype=float)
                    baseline_beta = np.linalg.lstsq(
                        baseline_design,
                        y_train,
                        rcond=None,
                    )[0]
                    extended_beta = np.linalg.lstsq(
                        extended_design,
                        y_train,
                        rcond=None,
                    )[0]
                    feature_index = extended_names.index(
                        "feature_per_training_sd"
                    )

                    for holdout_split in holdouts:
                        holdout = daily[
                            daily["split"].eq(holdout_split)
                            & feature_eligibility(daily, feature)
                        ].dropna(subset=required)
                        holdout = holdout[
                            target_within_split(holdout, outcome)
                        ]
                        if sample == "home_only":
                            travel_future = (
                                holdout["travel_t_plus_1"].eq(0)
                                if outcome == "next_day_calories"
                                else holdout["future_3d_travel_share"].eq(0)
                            )
                            holdout = holdout[
                                holdout["travel_today"].eq(0) & travel_future
                            ]
                        holdout = holdout.sort_values("date")
                        if len(holdout) < 50:
                            continue
                        baseline_x, _, _, _ = design_matrix(
                            holdout,
                            outcome,
                            None,
                            baseline_means,
                            baseline_stds,
                        )
                        extended_x, _, _, _ = design_matrix(
                            holdout,
                            outcome,
                            feature,
                            extended_means,
                            extended_stds,
                        )
                        observed = holdout[outcome].to_numpy(dtype=float)
                        baseline_error = observed - baseline_x @ baseline_beta
                        extended_error = observed - extended_x @ extended_beta
                        baseline_mae = float(np.mean(np.abs(baseline_error)))
                        extended_mae = float(np.mean(np.abs(extended_error)))
                        baseline_rmse = float(
                            np.sqrt(np.mean(baseline_error**2))
                        )
                        extended_rmse = float(
                            np.sqrt(np.mean(extended_error**2))
                        )
                        rows.append(
                            {
                                "regime": regime,
                                "sample": sample,
                                "feature": feature,
                                "outcome": outcome,
                                "training_split": training_split,
                                "holdout_split": holdout_split,
                                "train_n": len(train),
                                "holdout_n": len(holdout),
                                "training_feature_coefficient_per_sd": float(
                                    extended_beta[feature_index]
                                ),
                                "training_feature_sd_native": extended_stds[
                                    feature
                                ],
                                "baseline_mae": baseline_mae,
                                "extended_mae": extended_mae,
                                "delta_mae_extended_minus_baseline": (
                                    extended_mae - baseline_mae
                                ),
                                "baseline_rmse": baseline_rmse,
                                "extended_rmse": extended_rmse,
                                "delta_rmse_extended_minus_baseline": (
                                    extended_rmse - baseline_rmse
                                ),
                                "feature_improves_mae": (
                                    extended_mae < baseline_mae
                                ),
                            }
                        )
    return pd.DataFrame(rows)


def current_meal_summary(food_rows: pd.DataFrame) -> pd.DataFrame:
    """Show current calories/protein/fiber/liquids by diary meal label."""
    start = DATA_END - pd.Timedelta(days=CURRENT_WINDOW_DAYS - 1)
    frame = food_rows[food_rows["date"].between(start, DATA_END)].copy()
    rows: list[dict[str, object]] = []
    for meal, group in frame.groupby("meal", sort=True):
        daily_meal = (
            group.groupby("date", as_index=False)
            .agg(
                calories=("calories", "sum"),
                protein_g=("protein_g", "sum"),
                fiber_g=("fiber_g", "sum"),
                liquid_calories=(
                    "calories",
                    lambda values: float(
                        values[group.loc[values.index, "is_liquid"]].sum()
                    ),
                ),
                food_rows=("food", "size"),
            )
        )
        meal_kcal = float(group["calories"].sum())
        rows.append(
            {
                "meal": meal,
                "window_start": start.date(),
                "window_end": DATA_END.date(),
                "days_logged": group["date"].nunique(),
                "mean_calories_on_logged_days": daily_meal["calories"].mean(),
                "mean_protein_g_on_logged_days": daily_meal["protein_g"].mean(),
                "mean_fiber_g_on_logged_days": daily_meal["fiber_g"].mean(),
                "protein_g_per_1000_kcal": safe_share(
                    1000.0 * group["protein_g"].fillna(0).sum(),
                    meal_kcal,
                ),
                "fiber_g_per_1000_kcal": safe_share(
                    1000.0 * group["fiber_g"].fillna(0).sum(),
                    meal_kcal,
                ),
                "liquid_calorie_share": safe_share(
                    group.loc[group["is_liquid"], "calories"].sum(),
                    meal_kcal,
                ),
                "mean_food_rows_on_logged_days": daily_meal[
                    "food_rows"
                ].mean(),
                "share_of_window_calories": safe_share(
                    meal_kcal,
                    frame["calories"].sum(),
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(
        "share_of_window_calories",
        ascending=False,
    )


def liquid_food_summary(food_rows: pd.DataFrame) -> pd.DataFrame:
    """Expose the heuristic's largest classified foods for manual audit."""
    liquids = food_rows[food_rows["is_liquid"]].copy()
    return (
        liquids.groupby("food", as_index=False)
        .agg(
            days=("date", "nunique"),
            entries=("date", "size"),
            calories=("calories", "sum"),
            first_date=("date", "min"),
            last_date=("date", "max"),
            sources=("source", lambda values: ",".join(sorted(set(values)))),
            diary_units=(
                "diary_unit_string",
                lambda values: ",".join(
                    sorted(
                        {
                            str(value)
                            for value in values.dropna()
                            if str(value).strip()
                        }
                    )
                ),
            ),
        )
        .sort_values("calories", ascending=False)
    )


def print_report(
    coverage: pd.DataFrame,
    daily: pd.DataFrame,
    split_summary: pd.DataFrame,
    same_day: pd.DataFrame,
    predictions: pd.DataFrame,
    meals: pd.DataFrame,
) -> None:
    overall = coverage[
        (coverage["level"] == "all_years") & (coverage["label"] == "all")
    ].iloc[0]
    modern = split_summary[
        (split_summary["split"].isin(["discovery", "validation", "confirmation"]))
        & split_summary["context"].eq("all_context")
    ]
    current = split_summary[
        split_summary["split"].eq("treatment_holdout")
        & split_summary["context"].eq("all_context")
    ].iloc[0]
    print("ENERGY DENSITY / NUTRIENT LEVERAGE AUDIT")
    print("=" * 110)
    print(
        f"Ledger: {len(daily):,} days. Calories and meal labels are complete; "
        f"exact mass covers {overall['exact_mass_solid_calorie_coverage']:.1%} "
        "of solid calories overall."
    )
    print(
        "Modern pre-treatment exact-mass coverage by ordered split: "
        + ", ".join(
            f"{row.split} {row.calorie_weighted_exact_mass_coverage:.1%}"
            for row in modern.itertuples()
        )
    )
    print(
        "Treatment holdout means: "
        f"{current.mean_calories:.0f} kcal, "
        f"{current.mean_protein_g_per_1000_kcal:.1f} g protein/1000 kcal, "
        f"{current.mean_fiber_g_per_1000_kcal:.1f} g fiber/1000 kcal, "
        f"{current.mean_liquid_calorie_share:.1%} liquid calories."
    )

    print("\nCurrent meal-label summary")
    print(
        meals[
            [
                "meal",
                "days_logged",
                "mean_calories_on_logged_days",
                "mean_protein_g_on_logged_days",
                "mean_fiber_g_on_logged_days",
                "protein_g_per_1000_kcal",
                "fiber_g_per_1000_kcal",
                "share_of_window_calories",
            ]
        ]
        .round(2)
        .to_string(index=False)
    )

    print("\nLocked holdout MAE changes (negative is better)")
    shown = predictions[
        predictions["sample"].eq("all_context")
        & predictions["outcome"].eq("next_day_calories")
    ][
        [
            "regime",
            "feature",
            "holdout_split",
            "holdout_n",
            "delta_mae_extended_minus_baseline",
        ]
    ]
    print(shown.round(2).to_string(index=False))

    print("\nSame-day effects are descriptive and mathematically coupled.")
    print(
        same_day[
            same_day["sample"].eq("all_context")
            & same_day["split"].isin(
                ["discovery", "validation", "confirmation"]
            )
        ][
            [
                "split",
                "feature",
                "n",
                "effect_kcal_per_discovery_sd",
                "ci95_low",
                "ci95_high",
            ]
        ]
        .round(2)
        .to_string(index=False)
    )


def main() -> None:
    food_rows = classify_food_rows()
    coverage = summarize_coverage(food_rows)
    daily = build_daily_table(food_rows)
    split_summary = summarize_feature_splits(daily)
    same_day = same_day_markers(daily)
    predictions = prediction_results(daily)
    meals = current_meal_summary(food_rows)
    liquids = liquid_food_summary(food_rows)

    outputs = {
        "EN_data_coverage_by_year_source.csv": coverage,
        "EN_daily_features.csv": daily,
        "EN_feature_split_summary.csv": split_summary,
        "EN_same_day_markers.csv": same_day,
        "EN_locked_prediction_results.csv": predictions,
        "EN_current_meal_summary.csv": meals,
        "EN_liquid_food_summary.csv": liquids,
    }
    for filename, frame in outputs.items():
        frame.to_csv(ROOT / "analysis" / filename, index=False)
    print_report(
        coverage,
        daily,
        split_summary,
        same_day,
        predictions,
        meals,
    )


if __name__ == "__main__":
    main()
