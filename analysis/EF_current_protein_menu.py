"""Describe current protein distribution and calorie-neutral swap candidates.

The analysis covers the final 90 intake days and uses the foods already logged
by the subject. It does not infer protein quality or prescribe a target.

Run:

    uv run --with numpy --with pandas python analysis/EF_current_protein_menu.py

Outputs:
    analysis/EF_current_food_sources.csv
    analysis/EF_current_protein_by_meal.csv
    analysis/EF_current_swap_candidates.csv
    analysis/EF_current_menu_summary.csv
"""

from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
END_DATE = pd.Timestamp("2026-02-28")
WINDOW_DAYS = 90
START_DATE = END_DATE - pd.Timedelta(days=WINDOW_DAYS - 1)
MIN_DAYS_FOR_MENU_ITEM = 4
LOW_PROTEIN_DENSITY = 3.0
HIGH_PROTEIN_DENSITY = 7.0


def normalize_food_name(value: str) -> str:
    text = str(value).lower().replace("’", "'").strip()
    text = re.sub(r"\([^)]*\)", " ", text)
    text = re.sub(r",\s*\d+(?:\.\d+)?\s+.*$", "", text)
    text = re.sub(r"\b\d+(?:\.\d+)?\s*(?:g|kg|ml|oz|cup|cups|tbsp|tsp)\b", "", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip(" ,-")


def build_food_sources(foods: pd.DataFrame) -> pd.DataFrame:
    frame = foods.copy()
    frame["food_key"] = frame["food"].map(normalize_food_name)
    frame["calories"] = pd.to_numeric(frame["calories"], errors="coerce").fillna(0)
    frame["protein_g"] = pd.to_numeric(frame["protein_g"], errors="coerce").fillna(0)
    grouped = (
        frame.groupby("food_key")
        .agg(
            days=("date", "nunique"),
            entries=("date", "size"),
            total_calories=("calories", "sum"),
            total_protein_g=("protein_g", "sum"),
            representative_food=("food", lambda values: values.value_counts().index[0]),
        )
        .reset_index()
    )
    grouped["protein_g_per_100_kcal"] = np.where(
        grouped["total_calories"] > 0,
        100 * grouped["total_protein_g"] / grouped["total_calories"],
        np.nan,
    )
    grouped["share_of_window_calories"] = (
        grouped["total_calories"] / frame["calories"].sum()
    )
    grouped["share_of_window_protein"] = (
        grouped["total_protein_g"] / frame["protein_g"].sum()
    )
    return grouped.sort_values("total_calories", ascending=False)


def build_meal_summary(foods: pd.DataFrame) -> pd.DataFrame:
    daily_meal = (
        foods.groupby(["date", "meal"], as_index=False)
        .agg(calories=("calories", "sum"), protein_g=("protein_g", "sum"))
    )
    rows = []
    for meal, frame in daily_meal.groupby("meal"):
        rows.append(
            {
                "meal": meal,
                "days_logged": frame["date"].nunique(),
                "mean_calories_on_logged_days": frame["calories"].mean(),
                "mean_protein_g_on_logged_days": frame["protein_g"].mean(),
                "median_protein_g_on_logged_days": frame["protein_g"].median(),
                "days_at_least_20g_protein": int((frame["protein_g"] >= 20).sum()),
                "share_days_at_least_20g_protein": (
                    frame["protein_g"] >= 20
                ).mean(),
            }
        )
    return pd.DataFrame(rows).sort_values(
        "mean_protein_g_on_logged_days",
        ascending=False,
    )


def build_swap_candidates(sources: pd.DataFrame) -> pd.DataFrame:
    eligible = sources[sources["days"] >= MIN_DAYS_FOR_MENU_ITEM].copy()
    low = eligible[
        eligible["protein_g_per_100_kcal"] < LOW_PROTEIN_DENSITY
    ].nlargest(30, "total_calories")
    low["candidate_type"] = "recurring low-protein calorie source"

    high = eligible[
        eligible["protein_g_per_100_kcal"] >= HIGH_PROTEIN_DENSITY
    ].nlargest(30, "total_protein_g")
    high["candidate_type"] = "existing protein-dense familiar food"

    return pd.concat([low, high], ignore_index=True)[
        [
            "candidate_type",
            "food_key",
            "representative_food",
            "days",
            "total_calories",
            "total_protein_g",
            "protein_g_per_100_kcal",
            "share_of_window_calories",
            "share_of_window_protein",
        ]
    ]


def print_report(
    sources: pd.DataFrame,
    meals: pd.DataFrame,
    candidates: pd.DataFrame,
    foods: pd.DataFrame,
) -> None:
    daily = foods.groupby("date").agg(
        calories=("calories", "sum"),
        protein_g=("protein_g", "sum"),
    )
    print("CURRENT PROTEIN / MENU AUDIT")
    print("=" * 120)
    print(
        f"Window: {START_DATE.date()} to {END_DATE.date()} | "
        f"{len(daily)} days | mean {daily['calories'].mean():.0f} kcal, "
        f"{daily['protein_g'].mean():.1f} g protein"
    )

    print("\nMeal-level protein distribution")
    print(meals.round(2).to_string(index=False))

    print("\nLargest recurring low-protein calorie sources")
    low = candidates[
        candidates["candidate_type"] == "recurring low-protein calorie source"
    ]
    print(
        low[
            [
                "representative_food",
                "days",
                "total_calories",
                "total_protein_g",
                "protein_g_per_100_kcal",
            ]
        ].head(20).round(2).to_string(index=False)
    )

    print("\nExisting familiar protein-dense foods")
    high = candidates[
        candidates["candidate_type"] == "existing protein-dense familiar food"
    ]
    print(
        high[
            [
                "representative_food",
                "days",
                "total_calories",
                "total_protein_g",
                "protein_g_per_100_kcal",
            ]
        ].head(20).round(2).to_string(index=False)
    )

    top_protein_share = (
        sources.nlargest(10, "total_protein_g")["share_of_window_protein"].sum()
    )
    print(
        f"\nTop 10 food clusters provide {top_protein_share:.1%} of all protein. "
        "Use familiar foods for substitutions to avoid adding menu novelty."
    )


def build_window_summary(
    sources: pd.DataFrame,
    candidates: pd.DataFrame,
    foods: pd.DataFrame,
) -> pd.DataFrame:
    total_calories = foods["calories"].sum()
    total_protein = foods["protein_g"].sum()
    all_low = sources[
        (sources["days"] >= MIN_DAYS_FOR_MENU_ITEM)
        & (sources["protein_g_per_100_kcal"] < LOW_PROTEIN_DENSITY)
    ]
    top_ten_low = (
        candidates[
            candidates["candidate_type"]
            == "recurring low-protein calorie source"
        ]
        .nlargest(10, "total_calories")
    )
    return pd.DataFrame(
        [
            {
                "group": "full 90-day window",
                "calories": total_calories,
                "calories_per_day": total_calories / WINDOW_DAYS,
                "calorie_share": 1.0,
                "protein_g": total_protein,
                "protein_g_per_day": total_protein / WINDOW_DAYS,
            },
            {
                "group": "all recurring foods below 3g protein per 100 kcal",
                "calories": all_low["total_calories"].sum(),
                "calories_per_day": all_low["total_calories"].sum() / WINDOW_DAYS,
                "calorie_share": all_low["total_calories"].sum() / total_calories,
                "protein_g": all_low["total_protein_g"].sum(),
                "protein_g_per_day": all_low["total_protein_g"].sum()
                / WINDOW_DAYS,
            },
            {
                "group": "top 10 recurring low-protein calorie sources",
                "calories": top_ten_low["total_calories"].sum(),
                "calories_per_day": top_ten_low["total_calories"].sum()
                / WINDOW_DAYS,
                "calorie_share": top_ten_low["total_calories"].sum()
                / total_calories,
                "protein_g": top_ten_low["total_protein_g"].sum(),
                "protein_g_per_day": top_ten_low["total_protein_g"].sum()
                / WINDOW_DAYS,
            },
        ]
    )


def main() -> None:
    foods = pd.read_csv(
        ROOT / "intake" / "intake_foods.csv",
        parse_dates=["date"],
    )
    foods = foods[foods["date"].between(START_DATE, END_DATE)].copy()
    sources = build_food_sources(foods)
    meals = build_meal_summary(foods)
    candidates = build_swap_candidates(sources)
    summary = build_window_summary(sources, candidates, foods)

    sources.to_csv(
        ROOT / "analysis" / "EF_current_food_sources.csv",
        index=False,
    )
    meals.to_csv(
        ROOT / "analysis" / "EF_current_protein_by_meal.csv",
        index=False,
    )
    candidates.to_csv(
        ROOT / "analysis" / "EF_current_swap_candidates.csv",
        index=False,
    )
    summary.to_csv(
        ROOT / "analysis" / "EF_current_menu_summary.csv",
        index=False,
    )
    print_report(sources, meals, candidates, foods)
    print("\nWindow-level leverage")
    print(summary.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
