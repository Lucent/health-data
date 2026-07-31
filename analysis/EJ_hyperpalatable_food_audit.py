#!/usr/bin/env python3
"""Audit Fazzino-style hyper-palatable-food exposure against direct intake.

This is an observation-first analysis. It reads only the raw food/day ledgers,
the checked-in MFP nutrient enrichment, and checked-in travel dates. It does
not read P1, P2, P3, P4, drug PK, modeled TDEE, or modeled body composition.

The prespecified primary exposure is the calorie share of assessable solid
foods meeting the published fat + simple-sugar (FS) criterion:

    >20% kcal from fat AND >20% kcal from total sugar.

The published sodium clusters require consumed food weight. They are computed
only for entries logged directly in grams/ounces and are not used as the
primary exposure because food-weight missingness is extensive and non-random.
Liquids are excluded by a deterministic unit/name heuristic because the
published definition was developed for solid foods.

Run:

    uv run --with numpy --with pandas python analysis/EJ_hyperpalatable_food_audit.py

Outputs:
    analysis/EJ_hpf_daily.csv
    analysis/EJ_hpf_results.csv
    analysis/EJ_hpf_split_summary.csv
    analysis/EJ_hpf_food_summary.csv
"""

from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent

ANALYSIS_START = pd.Timestamp("2018-03-01")
DISCOVERY_END = pd.Timestamp("2021-12-31")
VALIDATION_END = pd.Timestamp("2022-12-31")
CONFIRMATION_END = pd.Timestamp("2024-09-16")

HIGH_INTAKE_KCAL = 3000
MIN_CURRENT_KCAL = 500
PRIMARY_COVERAGE_MIN = 0.80
HIGH_COVERAGE_MIN = 0.90
QUALITY_COVERAGE_MIN = 0.50
ATWATER_DRIFT_MAX = 0.20
HAC_LAGS = 7

FS_FAT_KCAL_PCT = 20.0
FS_SUGAR_KCAL_PCT = 20.0
FSOD_FAT_KCAL_PCT = 25.0
FSOD_SODIUM_WEIGHT_PCT = 0.30
CSOD_STARCH_KCAL_PCT = 40.0
CSOD_SODIUM_WEIGHT_PCT = 0.20

RESOLVED_PREFIX = "resolved"

DIRECT_GRAM_UNITS = re.compile(r"^(?:g|gram|grams|gram\(s\))\b", re.IGNORECASE)
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
    r"latte|cappuccino|espresso|coffee|tea|hot chocolate|water)\b",
    re.IGNORECASE,
)
LIQUID_NAME_EXCLUSIONS = re.compile(
    r"\b(?:coffee cake|tea cake|water cracker|water chestnut)\b",
    re.IGNORECASE,
)


def parse_travel_dates() -> set[pd.Timestamp]:
    """Parse the auditable checked-in trip table."""
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


def consumed_grams(row: pd.Series) -> tuple[float, str]:
    """Recover consumed grams only from explicit mass units in the diary."""
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
    """Apply a conservative physical-form heuristic required by the index."""
    unit = str(row.get("diary_unit_string") or "").strip()
    name = str(row.get("food") or "")
    by_unit = bool(LIQUID_UNITS.match(unit))
    by_name = bool(LIQUID_NAMES.search(name)) and not bool(
        LIQUID_NAME_EXCLUSIONS.search(name)
    )
    return by_unit or by_name


def classify_food_rows() -> pd.DataFrame:
    """Join enrichment and assign only directly computable HPF clusters."""
    foods = pd.read_csv(
        ROOT / "intake" / "intake_foods.csv",
        parse_dates=["date"],
    )
    food_dict = pd.read_csv(
        ROOT / "intake" / "mfp_food_dict.csv",
        low_memory=False,
    )
    data = foods.merge(
        food_dict,
        left_on="food",
        right_on="diary_food_name",
        how="left",
        suffixes=("_diary", "_dict"),
        validate="many_to_one",
    )

    numeric_columns = [
        "calories_diary",
        "calories_dict",
        "fat_g_diary",
        "fat_g_dict",
        "carbs_g_diary",
        "carbs_g_dict",
        "protein_g_diary",
        "protein_g_dict",
        "sugars_g",
        "sugar_g",
        "fiber_g_diary",
        "fiber_g_dict",
        "sodium_mg_diary",
        "sodium_mg_dict",
        "atwater_drift",
        "match_tier",
    ]
    for column in numeric_columns:
        data[column] = pd.to_numeric(data[column], errors="coerce")

    data["calories"] = data["calories_diary"].clip(lower=0)
    nutrient_pairs = {
        "fat_g_used": ("fat_g_diary", "fat_g_dict"),
        "carbs_g_used": ("carbs_g_diary", "carbs_g_dict"),
        "protein_g_used": ("protein_g_diary", "protein_g_dict"),
        "sugar_g_used": ("sugars_g", "sugar_g"),
        "fiber_g_used": ("fiber_g_diary", "fiber_g_dict"),
        "sodium_mg_used": ("sodium_mg_diary", "sodium_mg_dict"),
    }
    for output, (diary_column, enrichment_column) in nutrient_pairs.items():
        data[output] = data[diary_column].combine_first(
            data[enrichment_column]
        )
        data[f"{output}_source"] = np.select(
            [
                data[diary_column].notna(),
                data[enrichment_column].notna(),
            ],
            ["diary", "enrichment"],
            default="missing",
        )

    data["resolved"] = data["status"].str.startswith(
        RESOLVED_PREFIX,
        na=False,
    )
    data["liquid_excluded"] = data.apply(is_liquid, axis=1)
    mass = data.apply(consumed_grams, axis=1, result_type="expand")
    data["consumed_grams"] = mass[0]
    data["mass_source"] = mass[1]

    denominator = data["calories"].where(data["calories"] > 0)
    data["fat_kcal_pct"] = 900.0 * data["fat_g_used"] / denominator
    data["sugar_kcal_pct"] = 400.0 * data["sugar_g_used"] / denominator
    starch_g = (
        data["carbs_g_used"] - data["fiber_g_used"] - data["sugar_g_used"]
    )
    data["starch_g"] = starch_g
    data["starch_kcal_pct"] = 400.0 * starch_g / denominator
    data["sodium_weight_pct"] = (
        data["sodium_mg_used"] / (10.0 * data["consumed_grams"])
    )
    macro_kcal = (
        9.0 * data["fat_g_used"]
        + 4.0 * data["carbs_g_used"]
        + 4.0 * data["protein_g_used"]
    )
    data["selected_atwater_drift"] = (
        (macro_kcal - denominator).abs() / denominator
    )

    macro_values_valid = (
        ~data["liquid_excluded"]
        & data["calories"].gt(0)
        & data["fat_g_used"].ge(0)
        & data["sugar_g_used"].ge(0)
        & data["fat_kcal_pct"].between(0, 105)
        & data["sugar_kcal_pct"].between(0, 105)
    )
    data["fs_assessable"] = macro_values_valid
    data["fs"] = (
        data["fs_assessable"]
        & data["fat_kcal_pct"].gt(FS_FAT_KCAL_PCT)
        & data["sugar_kcal_pct"].gt(FS_SUGAR_KCAL_PCT)
    )
    data["fs_quality_assessable"] = (
        data["fs_assessable"]
        & data["selected_atwater_drift"].le(ATWATER_DRIFT_MAX)
    )
    data["fs_quality"] = data["fs"] & data["fs_quality_assessable"]

    sodium_base = (
        ~data["liquid_excluded"]
        & data["calories"].gt(0)
        & data["consumed_grams"].gt(0)
        & data["sodium_mg_used"].ge(0)
        & data["fat_kcal_pct"].between(0, 105)
        & data["sodium_weight_pct"].between(0, 10)
    )
    data["fsod_assessable"] = sodium_base
    data["fsod"] = (
        data["fsod_assessable"]
        & data["fat_kcal_pct"].gt(FSOD_FAT_KCAL_PCT)
        & data["sodium_weight_pct"].ge(FSOD_SODIUM_WEIGHT_PCT)
    )
    data["csod_assessable"] = (
        sodium_base
        & data["starch_g"].ge(0)
        & data["starch_kcal_pct"].between(0, 105)
    )
    data["csod"] = (
        data["csod_assessable"]
        & data["starch_kcal_pct"].gt(CSOD_STARCH_KCAL_PCT)
        & data["sodium_weight_pct"].ge(CSOD_SODIUM_WEIGHT_PCT)
    )
    data["any_hpf_confirmed"] = data["fs"] | data["fsod"] | data["csod"]
    data["solid_calories"] = data["calories"].where(
        ~data["liquid_excluded"],
        0.0,
    )
    return data


def weighted_calories(
    frame: pd.DataFrame,
    condition: str | pd.Series,
) -> float:
    if isinstance(condition, str):
        mask = frame[condition]
    else:
        mask = condition
    return float(frame.loc[mask, "solid_calories"].sum())


def build_daily_table(food_rows: pd.DataFrame) -> pd.DataFrame:
    """Aggregate exposure/coverage and attach direct daily outcomes."""
    records: list[dict[str, object]] = []
    for date, frame in food_rows.groupby("date", sort=True):
        solid_kcal = float(frame["solid_calories"].sum())
        fs_assessed_kcal = weighted_calories(frame, "fs_assessable")
        fs_quality_assessed_kcal = weighted_calories(
            frame,
            "fs_quality_assessable",
        )
        sodium_assessed_kcal = weighted_calories(
            frame,
            frame["fsod_assessable"] | frame["csod_assessable"],
        )
        source_energy = (
            frame.groupby("source_diary", dropna=False)["calories"]
            .sum()
            .to_dict()
        )
        records.append(
            {
                "date": date,
                "solid_calories": solid_kcal,
                "liquid_calories_excluded": float(
                    frame.loc[frame["liquid_excluded"], "calories"].sum()
                ),
                "fs_assessed_calories": fs_assessed_kcal,
                "fs_calories": weighted_calories(frame, "fs"),
                "fs_coverage": (
                    fs_assessed_kcal / solid_kcal if solid_kcal > 0 else np.nan
                ),
                "fs_share_assessed": (
                    weighted_calories(frame, "fs") / fs_assessed_kcal
                    if fs_assessed_kcal > 0
                    else np.nan
                ),
                "fs_share_total_solid": (
                    weighted_calories(frame, "fs") / solid_kcal
                    if solid_kcal > 0
                    else np.nan
                ),
                "fs_quality_assessed_calories": fs_quality_assessed_kcal,
                "fs_quality_calories": weighted_calories(frame, "fs_quality"),
                "fs_quality_coverage": (
                    fs_quality_assessed_kcal / solid_kcal
                    if solid_kcal > 0
                    else np.nan
                ),
                "fs_quality_share_assessed": (
                    weighted_calories(frame, "fs_quality")
                    / fs_quality_assessed_kcal
                    if fs_quality_assessed_kcal > 0
                    else np.nan
                ),
                "sodium_cluster_assessed_calories": sodium_assessed_kcal,
                "sodium_cluster_coverage": (
                    sodium_assessed_kcal / solid_kcal
                    if solid_kcal > 0
                    else np.nan
                ),
                "fsod_calories": weighted_calories(frame, "fsod"),
                "csod_calories": weighted_calories(frame, "csod"),
                "any_hpf_confirmed_calories": weighted_calories(
                    frame,
                    "any_hpf_confirmed",
                ),
                "any_hpf_confirmed_share_total_solid": (
                    weighted_calories(frame, "any_hpf_confirmed") / solid_kcal
                    if solid_kcal > 0
                    else np.nan
                ),
                "food_rows": len(frame),
                "resolved_rows": int(frame["resolved"].sum()),
                "approx_unverified_calories": float(
                    frame.loc[
                        frame["status"].eq("resolved_approx"),
                        "solid_calories",
                    ].sum()
                ),
                "source_html_calories": float(source_energy.get("html", 0.0)),
                "source_mhtml_calories": float(source_energy.get("mhtml", 0.0)),
                "source_oxps_calories": float(source_energy.get("oxps", 0.0)),
            }
        )
    exposure = pd.DataFrame(records)

    intake = pd.read_csv(
        ROOT / "intake" / "intake_daily.csv",
        parse_dates=["date"],
    )[["date", "calories"]]
    daily = (
        intake.merge(exposure, on="date", how="left", validate="one_to_one")
        .sort_values("date")
        .reset_index(drop=True)
    )
    travel_dates = parse_travel_dates()
    daily["travel_today"] = daily["date"].isin(travel_dates).astype(int)
    daily["travel_tomorrow"] = (
        daily["date"] + pd.Timedelta(days=1)
    ).isin(travel_dates).astype(int)
    daily["current_weekday"] = daily["date"].dt.dayofweek
    daily["target_weekday"] = (
        daily["date"] + pd.Timedelta(days=1)
    ).dt.dayofweek
    daily["day_of_year"] = daily["date"].dt.dayofyear
    daily["season_sin"] = np.sin(
        2 * np.pi * daily["day_of_year"] / 365.25
    )
    daily["season_cos"] = np.cos(
        2 * np.pi * daily["day_of_year"] / 365.25
    )
    daily["time_years"] = (
        (daily["date"] - ANALYSIS_START).dt.days / 365.25
    )
    daily["prior_28d_mean"] = (
        daily["calories"].shift(1).rolling(28, min_periods=14).mean()
    )
    daily["same_day_high_intake"] = (
        daily["calories"] >= HIGH_INTAKE_KCAL
    ).astype(float)
    daily["next_day_calories"] = daily["calories"].shift(-1)
    daily["next_day_high_intake"] = (
        daily["next_day_calories"] >= HIGH_INTAKE_KCAL
    ).astype(float)
    daily.loc[
        daily["next_day_calories"].isna(),
        "next_day_high_intake",
    ] = np.nan

    daily["split"] = "outside"
    daily.loc[
        daily["date"].between(ANALYSIS_START, DISCOVERY_END),
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

    discovery = daily[
        daily["split"].eq("discovery")
        & daily["calories"].ge(MIN_CURRENT_KCAL)
        & daily["fs_coverage"].ge(PRIMARY_COVERAGE_MIN)
    ]
    locked_cut = float(discovery["fs_share_assessed"].quantile(0.75))
    daily["high_fs_exposure_locked"] = (
        daily["fs_share_assessed"] >= locked_cut
    ).astype(float)
    daily.attrs["locked_fs_cut"] = locked_cut
    return daily


def add_column(
    columns: list[np.ndarray],
    names: list[str],
    name: str,
    values: pd.Series | np.ndarray,
) -> None:
    array = np.asarray(values, dtype=float)
    if name != "intercept" and (
        len(array) == 0 or np.nanstd(array) < 1e-12
    ):
        return
    columns.append(array)
    names.append(name)


def design_matrix(
    frame: pd.DataFrame,
    outcome: str,
    exposure: str,
    coverage: str,
) -> tuple[np.ndarray, list[str]]:
    """Build fixed, interpretable covariates for one split."""
    is_next_day = outcome.startswith("next_day")
    weekday_column = "target_weekday" if is_next_day else "current_weekday"

    columns: list[np.ndarray] = []
    names: list[str] = []
    add_column(columns, names, "intercept", np.ones(len(frame)))
    add_column(
        columns,
        names,
        "exposure_per_10pp",
        frame[exposure] * 10.0,
    )
    if is_next_day:
        add_column(
            columns,
            names,
            "current_calories_per_1000",
            frame["calories"] / 1000.0,
        )
    add_column(
        columns,
        names,
        "prior_28d_calories_per_1000",
        frame["prior_28d_mean"] / 1000.0,
    )
    add_column(
        columns,
        names,
        "coverage_per_10pp",
        frame[coverage] * 10.0,
    )
    add_column(columns, names, "time_years", frame["time_years"])
    add_column(columns, names, "season_sin", frame["season_sin"])
    add_column(columns, names, "season_cos", frame["season_cos"])
    add_column(columns, names, "travel_today", frame["travel_today"])
    if is_next_day:
        add_column(
            columns,
            names,
            "travel_tomorrow",
            frame["travel_tomorrow"],
        )
    for weekday in range(1, 7):
        add_column(
            columns,
            names,
            f"weekday_{weekday}",
            frame[weekday_column].eq(weekday).astype(float),
        )
    return np.column_stack(columns), names


def newey_west_covariance(
    design: np.ndarray,
    residuals: np.ndarray,
    lags: int,
) -> np.ndarray:
    """Deterministic Newey-West covariance for date-ordered observations."""
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


def partial_correlation(
    design: np.ndarray,
    names: list[str],
    outcome: np.ndarray,
) -> float:
    exposure_index = names.index("exposure_per_10pp")
    controls = np.delete(design, exposure_index, axis=1)
    exposure = design[:, exposure_index]
    outcome_residual = outcome - controls @ np.linalg.lstsq(
        controls,
        outcome,
        rcond=None,
    )[0]
    exposure_residual = exposure - controls @ np.linalg.lstsq(
        controls,
        exposure,
        rcond=None,
    )[0]
    if np.std(outcome_residual) == 0 or np.std(exposure_residual) == 0:
        return np.nan
    return float(np.corrcoef(outcome_residual, exposure_residual)[0, 1])


def fit_models(daily: pd.DataFrame) -> pd.DataFrame:
    """Fit prespecified continuous-exposure OLS/LPM models per ordered split."""
    specifications = [
        (
            "primary_all_context",
            "fs_share_assessed",
            "fs_coverage",
            PRIMARY_COVERAGE_MIN,
            False,
        ),
        (
            "home_only",
            "fs_share_assessed",
            "fs_coverage",
            PRIMARY_COVERAGE_MIN,
            True,
        ),
        (
            "coverage_90pct",
            "fs_share_assessed",
            "fs_coverage",
            HIGH_COVERAGE_MIN,
            False,
        ),
        (
            "atwater_quality",
            "fs_quality_share_assessed",
            "fs_quality_coverage",
            QUALITY_COVERAGE_MIN,
            False,
        ),
    ]
    outcomes = [
        "calories",
        "same_day_high_intake",
        "next_day_calories",
        "next_day_high_intake",
    ]
    rows: list[dict[str, object]] = []
    for sample, exposure, coverage, coverage_minimum, home_only in specifications:
        for split in ["discovery", "validation", "confirmation"]:
            for outcome in outcomes:
                required = [
                    exposure,
                    coverage,
                    outcome,
                    "prior_28d_mean",
                    "time_years",
                    "season_sin",
                    "season_cos",
                ]
                frame = daily[
                    daily["split"].eq(split)
                    & daily["calories"].ge(MIN_CURRENT_KCAL)
                    & daily[coverage].ge(coverage_minimum)
                ].dropna(subset=required)
                if home_only:
                    frame = frame[
                        frame["travel_today"].eq(0)
                        & frame["travel_tomorrow"].eq(0)
                    ]
                frame = frame.sort_values("date")
                if len(frame) < 50 or frame[exposure].nunique() < 3:
                    continue
                design, names = design_matrix(
                    frame,
                    outcome,
                    exposure,
                    coverage,
                )
                y = frame[outcome].to_numpy(dtype=float)
                beta = np.linalg.lstsq(design, y, rcond=None)[0]
                residuals = y - design @ beta
                covariance = newey_west_covariance(
                    design,
                    residuals,
                    HAC_LAGS,
                )
                exposure_index = names.index("exposure_per_10pp")
                effect = float(beta[exposure_index])
                standard_error = float(
                    np.sqrt(max(covariance[exposure_index, exposure_index], 0.0))
                )
                rows.append(
                    {
                        "sample": sample,
                        "split": split,
                        "outcome": outcome,
                        "n": len(frame),
                        "start": frame["date"].min().date(),
                        "end": frame["date"].max().date(),
                        "effect_per_10pp": effect,
                        "hac_se": standard_error,
                        "ci95_low": effect - 1.96 * standard_error,
                        "ci95_high": effect + 1.96 * standard_error,
                        "partial_r": partial_correlation(design, names, y),
                        "exposure_mean": float(frame[exposure].mean()),
                        "exposure_sd": float(frame[exposure].std(ddof=0)),
                        "coverage_mean": float(frame[coverage].mean()),
                        "high_intake_rate": float(
                            frame[
                                "next_day_high_intake"
                                if outcome.startswith("next_day")
                                else "same_day_high_intake"
                            ].mean()
                        ),
                        "controls": ",".join(
                            name
                            for name in names
                            if name not in {"intercept", "exposure_per_10pp"}
                        ),
                    }
                )
    return pd.DataFrame(rows)


def build_split_summary(
    daily: pd.DataFrame,
    locked_cut: float,
) -> pd.DataFrame:
    """Summarize coverage and a discovery-locked high-exposure contrast."""
    rows: list[dict[str, object]] = []
    for split in ["discovery", "validation", "confirmation"]:
        frame = daily[daily["split"].eq(split)].copy()
        eligible = frame[
            frame["calories"].ge(MIN_CURRENT_KCAL)
            & frame["fs_coverage"].ge(PRIMARY_COVERAGE_MIN)
        ]
        high = eligible[eligible["fs_share_assessed"].ge(locked_cut)]
        lower = eligible[eligible["fs_share_assessed"].lt(locked_cut)]
        rows.append(
            {
                "split": split,
                "start": frame["date"].min().date(),
                "end": frame["date"].max().date(),
                "calendar_days": len(frame),
                "eligible_days": len(eligible),
                "mean_fs_coverage": float(frame["fs_coverage"].mean()),
                "calorie_weighted_fs_coverage": float(
                    frame["fs_assessed_calories"].sum()
                    / frame["solid_calories"].sum()
                ),
                "mean_quality_coverage": float(
                    frame["fs_quality_coverage"].mean()
                ),
                "mean_sodium_cluster_coverage": float(
                    frame["sodium_cluster_coverage"].mean()
                ),
                "calorie_weighted_sodium_cluster_coverage": float(
                    frame["sodium_cluster_assessed_calories"].sum()
                    / frame["solid_calories"].sum()
                ),
                "mean_fs_share_assessed_eligible": float(
                    eligible["fs_share_assessed"].mean()
                ),
                "locked_high_fs_cut": locked_cut,
                "locked_high_days": len(high),
                "lower_days": len(lower),
                "same_day_calorie_delta_high_minus_lower": float(
                    high["calories"].mean() - lower["calories"].mean()
                ),
                "next_day_calorie_delta_high_minus_lower": float(
                    high["next_day_calories"].mean()
                    - lower["next_day_calories"].mean()
                ),
                "same_day_high_rate_high": float(
                    high["same_day_high_intake"].mean()
                ),
                "same_day_high_rate_lower": float(
                    lower["same_day_high_intake"].mean()
                ),
                "next_day_high_rate_high": float(
                    high["next_day_high_intake"].mean()
                ),
                "next_day_high_rate_lower": float(
                    lower["next_day_high_intake"].mean()
                ),
                "travel_days": int(frame["travel_today"].sum()),
            }
        )
    return pd.DataFrame(rows)


def build_food_summary(food_rows: pd.DataFrame) -> pd.DataFrame:
    """Aggregate the largest calorie sources and retain classification flags."""
    modern = food_rows[
        food_rows["date"].between(ANALYSIS_START, CONFIRMATION_END)
        & ~food_rows["liquid_excluded"]
    ].copy()
    summary = (
        modern.groupby("food", as_index=False)
        .agg(
            rows=("date", "size"),
            first_date=("date", "min"),
            last_date=("date", "max"),
            calories=("calories", "sum"),
            status=("status", "first"),
            match_tier=("match_tier", "first"),
            atwater_drift=("atwater_drift", "first"),
            selected_atwater_drift=("selected_atwater_drift", "first"),
            fat_source=("fat_g_used_source", "first"),
            sugar_source=("sugar_g_used_source", "first"),
            sodium_source=("sodium_mg_used_source", "first"),
            fs_assessable=("fs_assessable", "max"),
            fs=("fs", "max"),
            fs_quality_assessable=("fs_quality_assessable", "max"),
            fsod_assessable=("fsod_assessable", "max"),
            fsod=("fsod", "max"),
            csod_assessable=("csod_assessable", "max"),
            csod=("csod", "max"),
            fat_kcal_pct=("fat_kcal_pct", "first"),
            sugar_kcal_pct=("sugar_kcal_pct", "first"),
            starch_kcal_pct=("starch_kcal_pct", "first"),
            sodium_weight_pct=("sodium_weight_pct", "first"),
            mass_source=("mass_source", "first"),
        )
        .sort_values(["calories", "food"], ascending=[False, True])
        .reset_index(drop=True)
    )
    summary["first_date"] = summary["first_date"].dt.date
    summary["last_date"] = summary["last_date"].dt.date
    return summary


def print_report(
    daily: pd.DataFrame,
    results: pd.DataFrame,
    splits: pd.DataFrame,
    foods: pd.DataFrame,
) -> None:
    print("FAZZINO-STYLE HYPER-PALATABLE-FOOD AUDIT")
    print("=" * 112)
    print(
        f"Discovery {ANALYSIS_START.date()}–{DISCOVERY_END.date()} | "
        f"validation 2022 | confirmation 2023-01-01–{CONFIRMATION_END.date()}"
    )
    print(
        "Primary exposure: calorie share of assessable solid foods meeting "
        "fat+sugar (>20% kcal each)."
    )
    print(
        "No P1/P2/P3/P4, modeled TDEE/composition, or modeled drug level is read."
    )

    print("\nCoverage and discovery-locked high-exposure contrasts")
    columns = [
        "split",
        "calendar_days",
        "eligible_days",
        "calorie_weighted_fs_coverage",
        "calorie_weighted_sodium_cluster_coverage",
        "mean_fs_share_assessed_eligible",
        "locked_high_days",
        "same_day_calorie_delta_high_minus_lower",
        "next_day_calorie_delta_high_minus_lower",
        "next_day_high_rate_high",
        "next_day_high_rate_lower",
    ]
    print(splits[columns].round(3).to_string(index=False))

    print("\nPrimary adjusted effects per +10 percentage points FS exposure")
    primary = results[results["sample"].eq("primary_all_context")].copy()
    print(
        primary[
            [
                "outcome",
                "split",
                "n",
                "effect_per_10pp",
                "ci95_low",
                "ci95_high",
                "partial_r",
            ]
        ].round(3).to_string(index=False)
    )

    print("\nHome-only next-day sensitivity")
    home = results[
        results["sample"].eq("home_only")
        & results["outcome"].isin(
            ["next_day_calories", "next_day_high_intake"]
        )
    ]
    print(
        home[
            [
                "outcome",
                "split",
                "n",
                "effect_per_10pp",
                "ci95_low",
                "ci95_high",
                "partial_r",
            ]
        ].round(3).to_string(index=False)
    )

    print("\nLargest FS-classified calorie sources in the analysis period")
    top = foods[foods["fs"]].head(20)
    print(
        top[
            [
                "food",
                "rows",
                "calories",
                "fat_kcal_pct",
                "sugar_kcal_pct",
                "status",
                "selected_atwater_drift",
            ]
        ].round(3).to_string(index=False)
    )

    period = daily[daily["split"].ne("outside")]
    print("\nEnrichment cautions")
    print(
        "  Mean daily FS calorie coverage: "
        f"{period['fs_coverage'].mean():.1%}; "
        "calorie-weighted: "
        f"{period['fs_assessed_calories'].sum() / period['solid_calories'].sum():.1%}."
    )
    print(
        "  Calorie-weighted exact-mass coverage for sodium clusters: "
        f"{period['sodium_cluster_assessed_calories'].sum() / period['solid_calories'].sum():.1%}."
    )
    print(
        "  Approximate unverified matches, missing sugar, recipe aggregation, "
        "and MFP database errors remain possible."
    )


def main() -> None:
    food_rows = classify_food_rows()
    daily = build_daily_table(food_rows)
    locked_cut = float(daily.attrs["locked_fs_cut"])
    results = fit_models(daily)
    splits = build_split_summary(daily, locked_cut)
    foods = build_food_summary(food_rows)

    daily.to_csv(ROOT / "analysis" / "EJ_hpf_daily.csv", index=False)
    results.to_csv(ROOT / "analysis" / "EJ_hpf_results.csv", index=False)
    splits.to_csv(
        ROOT / "analysis" / "EJ_hpf_split_summary.csv",
        index=False,
    )
    foods.to_csv(ROOT / "analysis" / "EJ_hpf_food_summary.csv", index=False)
    print_report(daily, results, splits, foods)


if __name__ == "__main__":
    main()
