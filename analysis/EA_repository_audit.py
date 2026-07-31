"""Audit repository coverage, continuity, staleness, and nutrient enrichment.

This is intentionally descriptive. It does not modify source or derived data.
Run from any directory with:

    uv run --with pandas python analysis/EA_repository_audit.py
"""

from pathlib import Path
import re

import pandas as pd


ROOT = Path(__file__).resolve().parent.parent

DATED_FILES = {
    "daily intake": ("intake/intake_daily.csv", "date"),
    "food items": ("intake/intake_foods.csv", "date"),
    "weight": ("weight/weight.csv", "date"),
    "composition": ("composition/composition.csv", "date"),
    "RMR": ("RMR/rmr.csv", "date"),
    "steps": ("steps-sleep/steps.csv", "date"),
    "sleep": ("steps-sleep/sleep.csv", "date"),
    "cardio exercise": ("steps-sleep/exercises.csv", "date"),
    "strength": ("workout/strength.csv", "date"),
    "medicine": ("drugs/medicine.csv", "date"),
    "tirzepatide PK": ("drugs/tirzepatide.csv", "date"),
    "temperature": ("temperature/temperature.csv", "date"),
    "fatty-acid labs": ("fatty-acids/fatty_acids.csv", "date"),
    "baths": ("baths/baths.csv", "date"),
    "Kalman daily": ("analysis/P4_kalman_daily.csv", "date"),
}


def date_summary(label: str, relative_path: str, date_column: str) -> dict:
    path = ROOT / relative_path
    data = pd.read_csv(path)
    dates = pd.to_datetime(data[date_column], errors="coerce")
    return {
        "label": label,
        "path": relative_path,
        "rows": len(data),
        "first": dates.min(),
        "last": dates.max(),
        "bad_dates": int(dates.isna().sum()),
        "duplicate_rows": int(data.duplicated().sum()),
        "duplicate_dates": int(dates.duplicated().sum()),
    }


def latest_dated_filename(directory: Path, suffix: str) -> pd.Timestamp:
    dates = []
    for path in directory.glob(f"*{suffix}"):
        match = re.search(r"((?:19|20)\d{2})[-_]?(\d{2})[-_]?(\d{2})", path.stem)
        if not match:
            continue
        dates.append(pd.Timestamp("-".join(match.groups())))
    return max(dates) if dates else pd.NaT


def latest_inbody_measurement(directory: Path) -> pd.Timestamp:
    latest = pd.NaT
    for path in directory.glob("InBody-*.csv"):
        data = pd.read_csv(path, encoding="utf-8-sig", usecols=["date"])
        dates = pd.to_datetime(
            data["date"].astype(str).str.slice(0, 8),
            format="%Y%m%d",
            errors="coerce",
        )
        candidate = dates.max()
        if pd.isna(latest) or (pd.notna(candidate) and candidate > latest):
            latest = candidate
    return latest


def fmt_date(value: pd.Timestamp) -> str:
    if pd.isna(value):
        return "n/a"
    return value.strftime("%Y-%m-%d")


def print_source_coverage() -> pd.DataFrame:
    summaries = pd.DataFrame(
        date_summary(label, path, column)
        for label, (path, column) in DATED_FILES.items()
    )
    print("SOURCE COVERAGE")
    print("=" * 112)
    print(
        f"{'source':<20} {'rows':>8} {'first':>12} {'last':>12} "
        f"{'bad dates':>10} {'dup rows':>10} {'dup dates':>10}"
    )
    for row in summaries.itertuples():
        print(
            f"{row.label:<20} {row.rows:>8,} {fmt_date(row.first):>12} "
            f"{fmt_date(row.last):>12} {row.bad_dates:>10,} "
            f"{row.duplicate_rows:>10,} {row.duplicate_dates:>10,}"
        )
    return summaries


def print_intake_continuity() -> None:
    intake = pd.read_csv(ROOT / "intake/intake_daily.csv", parse_dates=["date"])
    foods = pd.read_csv(ROOT / "intake/intake_foods.csv", parse_dates=["date"])
    expected = pd.date_range(intake["date"].min(), intake["date"].max(), freq="D")
    missing = expected.difference(pd.DatetimeIndex(intake["date"]))

    item_totals = foods.groupby("date", as_index=False)["calories"].sum()
    compare = intake[["date", "calories"]].merge(
        item_totals,
        on="date",
        how="outer",
        suffixes=("_daily", "_items"),
        indicator=True,
    )
    compare["absolute_gap"] = (
        compare["calories_daily"] - compare["calories_items"]
    ).abs()

    print("\nINTAKE CONTINUITY AND DAILY/ITEM RECONCILIATION")
    print("=" * 112)
    print(f"Calendar days in range: {len(expected):,}")
    print(f"Rows in intake_daily: {len(intake):,}")
    print(f"Missing calendar dates: {len(missing):,}")
    if len(missing):
        print("Missing dates:", ", ".join(date.strftime("%Y-%m-%d") for date in missing))
    print(f"Days at or below 10 logged kcal: {(intake['calories'] <= 10).sum():,}")
    print(f"Days at or below 100 logged kcal: {(intake['calories'] <= 100).sum():,}")
    print(f"Food-item rows: {len(foods):,}")
    print(f"Unique food-item dates: {foods['date'].nunique():,}")
    print(f"Dates missing from either side of item/daily merge: {(compare['_merge'] != 'both').sum():,}")
    print(f"Median absolute daily-vs-item calorie gap: {compare['absolute_gap'].median():.1f}")
    print(f"Days with >1 kcal gap: {(compare['absolute_gap'] > 1).sum():,}")
    print(f"Days with >50 kcal gap: {(compare['absolute_gap'] > 50).sum():,}")


def print_artifact_staleness(summaries: pd.DataFrame) -> None:
    last = summaries.set_index("label")["last"]
    latest_workout_pdf = latest_dated_filename(ROOT / "workout", ".pdf")
    latest_inbody_measurement_date = latest_inbody_measurement(
        ROOT / "composition"
    )

    checks = [
        (
            "strength.csv vs workout PDFs",
            last["strength"],
            latest_workout_pdf,
        ),
        (
            "composition.csv vs InBody measurements",
            last["composition"],
            latest_inbody_measurement_date,
        ),
        (
            "Kalman output vs intake",
            last["Kalman daily"],
            last["daily intake"],
        ),
        (
            "Kalman output vs weight",
            last["Kalman daily"],
            last["weight"],
        ),
        (
            "tirzepatide PK vs medicine log",
            last["tirzepatide PK"],
            last["medicine"],
        ),
    ]

    print("\nDERIVED-ARTIFACT STALENESS CHECK")
    print("=" * 112)
    print(f"{'check':<48} {'derived last':>14} {'source last':>14} {'lag days':>10}")
    for label, derived_last, source_last in checks:
        if pd.isna(source_last) or pd.isna(derived_last):
            lag_text = "n/a"
        else:
            lag_text = str(
                (source_last.normalize() - derived_last.normalize()).days
            )
        print(
            f"{label:<48} {fmt_date(derived_last):>14} "
            f"{fmt_date(source_last):>14} {lag_text:>10}"
        )


def print_enrichment_coverage() -> None:
    foods = pd.read_csv(ROOT / "intake/intake_foods.csv", parse_dates=["date"])
    enriched = pd.read_csv(ROOT / "intake/mfp_food_dict.csv")

    duplicate_keys = enriched["diary_food_name"].duplicated(keep=False)
    one_per_key = (
        enriched.sort_values(["diary_food_name", "fetched_at"])
        .drop_duplicates("diary_food_name", keep="last")
    )
    joined = foods.merge(
        one_per_key,
        left_on="food",
        right_on="diary_food_name",
        how="left",
        suffixes=("_diary", "_enriched"),
    )

    joined["diary_calories"] = pd.to_numeric(
        joined["calories_diary"], errors="coerce"
    ).fillna(0)
    joined["resolved_mapping"] = joined["status"].fillna("").str.startswith(
        "resolved"
    )
    total_calories = joined["diary_calories"].sum()

    print("\nMFP ENRICHMENT COVERAGE")
    print("=" * 112)
    print(f"Rows in mfp_food_dict.csv: {len(enriched):,}")
    print(f"Unique diary_food_name keys: {enriched['diary_food_name'].nunique():,}")
    print(f"Rows whose key is duplicated: {duplicate_keys.sum():,}")
    print(
        "Diary item rows linked to a dictionary record: "
        f"{joined['diary_food_name'].notna().sum():,} / {len(joined):,}"
    )
    print(
        "Diary rows with a resolved MFP mapping: "
        f"{joined['resolved_mapping'].mean():.1%}"
    )
    print(
        "Diary calories with a resolved MFP mapping: "
        f"{joined.loc[joined['resolved_mapping'], 'diary_calories'].sum() / total_calories:.1%}"
    )

    for nutrient in [
        "saturated_fat_g",
        "mono_fat_g",
        "poly_fat_g",
        "trans_fat_g",
        "potassium_mg",
        "added_sugars_g",
    ]:
        covered_rows = joined[nutrient].notna()
        covered_calories = joined.loc[covered_rows, "diary_calories"].sum()
        print(
            f"{nutrient:<20} rows={covered_rows.mean():>6.1%} "
            f"calories={covered_calories / total_calories:>6.1%}"
        )

    yearly = (
        joined.assign(year=joined["date"].dt.year)
        .groupby("year")
        .apply(
            lambda group: pd.Series(
                {
                    "rows": len(group),
                    "resolved_share": group["resolved_mapping"].mean(),
                    "sat_fat_row_share": group["saturated_fat_g"].notna().mean(),
                    "poly_fat_row_share": group["poly_fat_g"].notna().mean(),
                }
            ),
            include_groups=False,
        )
    )
    print("\nYearly enriched-row coverage:")
    print(
        f"{'year':<6} {'rows':>7} {'resolved':>10} "
        f"{'sat fat':>10} {'poly fat':>10}"
    )
    for year, row in yearly.iterrows():
        print(
            f"{year:<6} {int(row['rows']):>7,} "
            f"{row['resolved_share']:>9.1%} "
            f"{row['sat_fat_row_share']:>9.1%} "
            f"{row['poly_fat_row_share']:>9.1%}"
        )

    print("\nEnrichment status counts:")
    print(
        one_per_key["status"]
        .fillna("(missing)")
        .value_counts(dropna=False)
        .to_string()
    )


def main() -> None:
    summaries = print_source_coverage()
    print_intake_continuity()
    print_artifact_staleness(summaries)
    print_enrichment_coverage()


if __name__ == "__main__":
    main()
