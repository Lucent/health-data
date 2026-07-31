#!/usr/bin/env python3
"""Build a deterministic, observation-first clinical and nutrition snapshot.

This script does not diagnose disease or impute missing clinical labs.  It
summarizes direct repository observations and explicitly labels database-
enriched nutrient estimates.
"""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
CURRENT_WINDOW_START = pd.Timestamp("2025-12-01")
CURRENT_WINDOW_END = pd.Timestamp("2026-02-28")

YEAR_OUTPUT = ROOT / "analysis" / "EH_nutrition_year_summary.csv"
SOURCE_OUTPUT = ROOT / "analysis" / "EH_top_saturated_fat_sources.csv"
SNAPSHOT_OUTPUT = ROOT / "analysis" / "EH_clinical_snapshot.csv"


def _number(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _latest_food_dictionary() -> pd.DataFrame:
    food_dict = pd.read_csv(ROOT / "intake" / "mfp_food_dict.csv", low_memory=False)
    food_dict["fetched_at"] = pd.to_datetime(
        food_dict["fetched_at"], errors="coerce", utc=True
    )
    return (
        food_dict.sort_values("fetched_at", na_position="first")
        .drop_duplicates("diary_food_name", keep="last")
        .rename(columns={"calories": "dictionary_calories"})
    )


def _nutrition_outputs() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float]]:
    daily = pd.read_csv(ROOT / "intake" / "intake_daily.csv")
    foods = pd.read_csv(ROOT / "intake" / "intake_foods.csv")
    food_dict = _latest_food_dictionary()

    daily["date"] = pd.to_datetime(daily["date"])
    foods["date"] = pd.to_datetime(foods["date"])
    foods = foods.rename(columns={"food": "diary_food_name", "calories": "diary_calories"})

    enrichment_columns = [
        "diary_food_name",
        "status",
        "saturated_fat_g",
        "potassium_mg",
        "fetched_at",
    ]
    enriched = foods.merge(
        food_dict[enrichment_columns], on="diary_food_name", how="left", validate="many_to_one"
    )
    enriched["diary_calories"] = _number(enriched["diary_calories"])
    enriched["saturated_fat_g"] = _number(enriched["saturated_fat_g"])
    enriched["potassium_mg"] = _number(enriched["potassium_mg"])
    enriched["year"] = enriched["date"].dt.year
    daily["year"] = daily["date"].dt.year

    direct_year = (
        daily.groupby("year", as_index=False)
        .agg(
            days=("date", "size"),
            mean_calories=("calories", "mean"),
            mean_sodium_mg=("sodium_mg", "mean"),
            mean_fiber_g=("fiber_g", "mean"),
            mean_protein_g=("protein_g", "mean"),
        )
    )
    known_sat_year = (
        enriched.groupby("year", as_index=False)
        .agg(
            known_saturated_fat_g=("saturated_fat_g", "sum"),
            saturated_fat_covered_calories=(
                "diary_calories",
                lambda values: values[enriched.loc[values.index, "saturated_fat_g"].notna()].sum(),
            ),
            food_calories=("diary_calories", "sum"),
        )
    )
    annual = direct_year.merge(known_sat_year, on="year", how="left")
    annual["known_saturated_fat_g_per_day"] = (
        annual["known_saturated_fat_g"] / annual["days"]
    )
    annual["known_saturated_fat_pct_energy"] = (
        annual["known_saturated_fat_g_per_day"] * 9 / annual["mean_calories"] * 100
    )
    annual["saturated_fat_kcal_coverage_pct"] = (
        annual["saturated_fat_covered_calories"] / annual["food_calories"] * 100
    )
    annual = annual[
        [
            "year",
            "days",
            "mean_calories",
            "mean_sodium_mg",
            "mean_fiber_g",
            "mean_protein_g",
            "known_saturated_fat_g_per_day",
            "known_saturated_fat_pct_energy",
            "saturated_fat_kcal_coverage_pct",
        ]
    ].round(2)
    annual.to_csv(YEAR_OUTPUT, index=False)

    current_daily = daily[
        daily["date"].between(CURRENT_WINDOW_START, CURRENT_WINDOW_END)
    ].copy()
    current_foods = enriched[
        enriched["date"].between(CURRENT_WINDOW_START, CURRENT_WINDOW_END)
    ].copy()
    current_days = int(current_daily["date"].nunique())
    known_sat_g = float(current_foods["saturated_fat_g"].sum())
    total_food_kcal = float(current_foods["diary_calories"].sum())
    covered_food_kcal = float(
        current_foods.loc[current_foods["saturated_fat_g"].notna(), "diary_calories"].sum()
    )

    sources = (
        current_foods[current_foods["saturated_fat_g"].notna()]
        .groupby("diary_food_name", as_index=False)
        .agg(
            days=("date", "nunique"),
            entries=("date", "size"),
            saturated_fat_g=("saturated_fat_g", "sum"),
            calories=("diary_calories", "sum"),
        )
    )
    sources["saturated_fat_g_per_day"] = sources["saturated_fat_g"] / current_days
    sources["share_of_known_saturated_fat_pct"] = (
        sources["saturated_fat_g"] / known_sat_g * 100
    )
    sources = sources.sort_values("saturated_fat_g", ascending=False).round(2)
    sources.to_csv(SOURCE_OUTPUT, index=False)

    current = {
        "days": current_days,
        "mean_calories": float(current_daily["calories"].mean()),
        "mean_sodium_mg": float(current_daily["sodium_mg"].mean()),
        "mean_fiber_g": float(current_daily["fiber_g"].mean()),
        "mean_protein_g": float(current_daily["protein_g"].mean()),
        "known_saturated_fat_g_per_day": known_sat_g / current_days,
        "known_saturated_fat_pct_energy": (
            known_sat_g / current_days * 9 / current_daily["calories"].mean() * 100
        ),
        "saturated_fat_kcal_coverage_pct": covered_food_kcal / total_food_kcal * 100,
        "known_potassium_mg_per_day": float(current_foods["potassium_mg"].sum())
        / current_days,
        "potassium_kcal_coverage_pct": float(
            current_foods.loc[
                current_foods["potassium_mg"].notna(), "diary_calories"
            ].sum()
        )
        / total_food_kcal
        * 100,
    }
    return annual, sources, current


def _clinical_snapshot(current: dict[str, float]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    def add(
        metric: str,
        value: float,
        unit: str,
        date_start: str,
        date_end: str,
        source: str,
        evidence_note: str,
    ) -> None:
        rows.append(
            {
                "metric": metric,
                "value": round(float(value), 3),
                "unit": unit,
                "date_start": date_start,
                "date_end": date_end,
                "source": source,
                "evidence_note": evidence_note,
            }
        )

    window_start = CURRENT_WINDOW_START.date().isoformat()
    window_end = CURRENT_WINDOW_END.date().isoformat()
    direct_daily_metrics = [
        ("daily_calories_mean", "mean_calories", "kcal/day"),
        ("daily_sodium_mean", "mean_sodium_mg", "mg/day"),
        ("daily_fiber_mean", "mean_fiber_g", "g/day"),
        ("daily_protein_mean", "mean_protein_g", "g/day"),
    ]
    for metric, key, unit in direct_daily_metrics:
        add(
            metric,
            current[key],
            unit,
            window_start,
            window_end,
            "intake/intake_daily.csv",
            "direct diary total; recording accuracy is not independently verified",
        )
    for metric, key, unit in [
        ("known_saturated_fat_mean", "known_saturated_fat_g_per_day", "g/day"),
        ("known_saturated_fat_energy", "known_saturated_fat_pct_energy", "% kcal"),
        (
            "saturated_fat_calorie_coverage",
            "saturated_fat_kcal_coverage_pct",
            "% diary kcal",
        ),
        ("known_potassium_mean", "known_potassium_mg_per_day", "mg/day"),
        ("potassium_calorie_coverage", "potassium_kcal_coverage_pct", "% diary kcal"),
    ]:
        add(
            metric,
            current[key],
            unit,
            window_start,
            window_end,
            "intake/intake_foods.csv + intake/mfp_food_dict.csv",
            "database-enriched estimate; missing foods are not zero and values are not biomarkers",
        )

    composition = pd.read_csv(ROOT / "composition" / "composition.csv")
    composition["date"] = pd.to_datetime(composition["date"])
    h30 = composition[composition["device"].eq("H30")].sort_values("date")
    for label, row in [("h30_first", h30.iloc[0]), ("h30_latest", h30.iloc[-1])]:
        date = row["date"].date().isoformat()
        for metric, column, unit in [
            ("weight", "weight_lbs", "lb"),
            ("bmi", "bmi", "kg/m2"),
            ("fat_mass", "fat_mass_lbs", "lb"),
            ("fat_pct", "fat_pct", "%"),
            ("lean_mass", "lean_mass_lbs", "lb"),
            ("skeletal_muscle_mass", "smm_lbs", "lb"),
            ("visceral_fat_level", "visceral_fat_level", "device level"),
            ("waist_hip_ratio", "waist_hip_ratio", "ratio"),
        ]:
            if pd.notna(row[column]):
                add(
                    f"{label}_{metric}",
                    row[column],
                    unit,
                    date,
                    date,
                    "composition/composition.csv",
                    "direct device reading; bioimpedance-derived composition is hydration-sensitive",
                )

    rmr = pd.read_csv(ROOT / "RMR" / "rmr.csv")
    rmr["date"] = pd.to_datetime(rmr["date"])
    latest_rmr = rmr.sort_values("date").iloc[-1]
    latest_rmr_date = latest_rmr["date"].date().isoformat()
    add(
        "latest_measured_rmr",
        latest_rmr["rmr_kcal"],
        "kcal/day",
        latest_rmr_date,
        latest_rmr_date,
        "RMR/rmr.csv",
        f"direct {latest_rmr['device']} test; single-test uncertainty applies",
    )
    comparison_rmr = rmr[rmr["date"].isin(pd.to_datetime(["2024-09-17", "2024-09-19"]))]
    if len(comparison_rmr) == 2:
        add(
            "2024_same_week_rmr_mean",
            comparison_rmr["rmr_kcal"].mean(),
            "kcal/day",
            "2024-09-17",
            "2024-09-19",
            "RMR/rmr.csv",
            "mean of two direct cosmed_fitmate tests",
        )

    fatty = pd.read_csv(ROOT / "fatty-acids" / "fatty_acids.csv")
    fatty["date"] = pd.to_datetime(fatty["date"])
    latest_fatty_date = fatty["date"].max()
    latest_fatty = fatty[fatty["date"].eq(latest_fatty_date)]
    for label, analyte in [
        ("omega_3_index", "Omega-3 Index (RBC)"),
        ("epa", "EPA (20:5n3)"),
        ("dha", "DHA (22:6n3)"),
        ("arachidonic_acid", "Arachidonic Acid (AA 20:4n6)"),
        ("aa_epa_ratio", "AA:EPA ratio"),
        ("omega_6_omega_3_ratio", "Omega-6:Omega-3 ratio"),
    ]:
        match = latest_fatty[latest_fatty["fatty_acid"].eq(analyte)]
        if not match.empty:
            row = match.iloc[0]
            date = latest_fatty_date.date().isoformat()
            add(
                f"latest_{label}",
                row["value_pct"],
                row["unit"],
                date,
                date,
                "fatty-acids/fatty_acids.csv",
                "manually extracted lab value; source PDF is absent from the repository",
            )

    heart_rate = pd.read_csv(ROOT / "analysis" / "BX_tirz_heart_rate_summary.csv")
    for metric in ["night_median_hr_pre_tirz", "night_median_hr_on_tirz"]:
        match = heart_rate[heart_rate["metric"].eq(metric)]
        if not match.empty:
            add(
                metric,
                match.iloc[0]["median"],
                "bpm",
                "",
                "",
                "analysis/BX_tirz_heart_rate_summary.csv",
                "derived nightly proxy; underlying external Samsung archive is not in this repository",
            )

    temperature = pd.read_csv(ROOT / "temperature" / "temperature.csv")
    temperature["date"] = pd.to_datetime(temperature["date"])
    add(
        "temperature_reading_count",
        len(temperature),
        "readings",
        temperature["date"].min().date().isoformat(),
        temperature["date"].max().date().isoformat(),
        "temperature/temperature.csv",
        "direct logged values; measurement site, device, and standardized protocol are not recorded",
    )

    medicine = pd.read_csv(ROOT / "drugs" / "medicine.csv")
    medicine["date"] = pd.to_datetime(medicine["date"])
    click_rows = medicine[_number(medicine["clicks"]).notna()].sort_values("date")
    add(
        "medicine_log_row_count",
        len(medicine),
        "records",
        medicine["date"].min().date().isoformat(),
        medicine["date"].max().date().isoformat(),
        "drugs/medicine.csv",
        "observed records only; current regimen is not established by the last log entry",
    )
    if len(click_rows) >= 2:
        gaps = click_rows["date"].diff().dt.days
        add(
            "click_dose_three_day_intervals",
            (gaps == 3).sum(),
            "intervals",
            click_rows["date"].min().date().isoformat(),
            click_rows["date"].max().date().isoformat(),
            "drugs/medicine.csv",
            "observed record only; click-count dosing is not operationalized or endorsed",
        )

    snapshot = pd.DataFrame(rows)
    snapshot.to_csv(SNAPSHOT_OUTPUT, index=False)
    return snapshot


def main() -> None:
    annual, sources, current = _nutrition_outputs()
    snapshot = _clinical_snapshot(current)
    print(f"Wrote {YEAR_OUTPUT.relative_to(ROOT)} ({len(annual)} rows)")
    print(f"Wrote {SOURCE_OUTPUT.relative_to(ROOT)} ({len(sources)} rows)")
    print(f"Wrote {SNAPSHOT_OUTPUT.relative_to(ROOT)} ({len(snapshot)} rows)")
    print(
        "Current-window known saturated fat: "
        f"{current['known_saturated_fat_g_per_day']:.1f} g/day, "
        f"{current['known_saturated_fat_pct_energy']:.1f}% kcal, "
        f"{current['saturated_fat_kcal_coverage_pct']:.1f}% calorie coverage"
    )


if __name__ == "__main__":
    main()
