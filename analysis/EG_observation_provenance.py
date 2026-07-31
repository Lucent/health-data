"""Build a daily observation/provenance table without latent-state estimates.

The table makes it explicit which values are observed, which source/device
regime applies, and how far each date is from a direct measurement.

Run:

    uv run --with numpy --with pandas python analysis/EG_observation_provenance.py

Output:
    analysis/EG_observation_provenance_daily.csv
"""

from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
VALID_SLEEP_MIN_HOURS = 2
VALID_SLEEP_MAX_HOURS = 11


def nearest_distance_days(
    target_dates: pd.Series,
    observed_dates: pd.Series,
) -> np.ndarray:
    target = target_dates.to_numpy(dtype="datetime64[D]").astype("int64")
    observed = np.sort(
        observed_dates.dropna()
        .to_numpy(dtype="datetime64[D]")
        .astype("int64")
    )
    if len(observed) == 0:
        return np.full(len(target), np.nan)
    insertion = np.searchsorted(observed, target)
    left_index = np.clip(insertion - 1, 0, len(observed) - 1)
    right_index = np.clip(insertion, 0, len(observed) - 1)
    left_distance = np.abs(target - observed[left_index])
    right_distance = np.abs(target - observed[right_index])
    return np.minimum(left_distance, right_distance).astype(float)


def parse_travel_dates() -> set[pd.Timestamp]:
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


def measurement_tier(days: float) -> str:
    if pd.isna(days):
        return "no measurement"
    if days == 0:
        return "same day"
    if days <= 7:
        return "within 7 days"
    if days <= 30:
        return "within 30 days"
    if days <= 90:
        return "within 90 days"
    return "more than 90 days"


def main() -> None:
    intake = pd.read_csv(
        ROOT / "intake" / "intake_daily.csv",
        parse_dates=["date"],
    )
    foods = pd.read_csv(
        ROOT / "intake" / "intake_foods.csv",
        parse_dates=["date"],
        usecols=["date", "source"],
    )
    weights = pd.read_csv(
        ROOT / "weight" / "weight.csv",
        parse_dates=["date"],
    )
    composition = pd.read_csv(
        ROOT / "composition" / "composition.csv",
        parse_dates=["date"],
    )
    rmr = pd.read_csv(ROOT / "RMR" / "rmr.csv", parse_dates=["date"])
    steps = pd.read_csv(
        ROOT / "steps-sleep" / "steps.csv",
        parse_dates=["date"],
    )
    samsung_steps = pd.read_csv(
        ROOT / "steps-sleep" / "steps_samsung.csv",
        parse_dates=["date"],
    )
    sleep = pd.read_csv(
        ROOT / "steps-sleep" / "sleep.csv",
        parse_dates=["date"],
    )
    strength = pd.read_csv(
        ROOT / "workout" / "strength.csv",
        parse_dates=["date"],
    )
    medicine = pd.read_csv(
        ROOT / "drugs" / "medicine.csv",
        parse_dates=["date"],
    )
    temperature = pd.read_csv(
        ROOT / "temperature" / "temperature.csv",
        parse_dates=["date"],
    )

    food_sources = (
        foods.groupby("date")
        .agg(
            intake_item_rows=("source", "size"),
            intake_source=("source", lambda values: "+".join(sorted(set(values)))),
        )
        .reset_index()
    )
    daily = intake.merge(food_sources, on="date", how="left")

    daily = daily.merge(
        weights[["date", "weight_lbs"]].rename(
            columns={"weight_lbs": "observed_weight_lbs"}
        ),
        on="date",
        how="left",
    )
    daily["weight_observed"] = daily["observed_weight_lbs"].notna()
    daily["days_to_nearest_weight"] = nearest_distance_days(
        daily["date"],
        weights["date"],
    )
    daily["weight_observation_tier"] = daily["days_to_nearest_weight"].map(
        measurement_tier
    )

    composition_columns = [
        "date",
        "device",
        "fat_mass_lbs",
        "lean_mass_lbs",
    ]
    daily = daily.merge(
        composition[composition_columns].rename(
            columns={
                "device": "composition_device",
                "fat_mass_lbs": "observed_fat_mass_lbs",
                "lean_mass_lbs": "observed_lean_mass_lbs",
            }
        ),
        on="date",
        how="left",
    )
    daily["composition_observed"] = daily["composition_device"].notna()
    daily["days_to_nearest_composition"] = nearest_distance_days(
        daily["date"],
        composition["date"],
    )
    daily["composition_observation_tier"] = daily[
        "days_to_nearest_composition"
    ].map(measurement_tier)

    daily = daily.merge(
        rmr[["date", "rmr_kcal"]].rename(
            columns={"rmr_kcal": "observed_rmr_kcal"}
        ),
        on="date",
        how="left",
    )
    daily["rmr_observed"] = daily["observed_rmr_kcal"].notna()
    daily["days_to_nearest_rmr"] = nearest_distance_days(
        daily["date"],
        rmr["date"],
    )

    daily = daily.merge(
        steps[["date", "steps"]],
        on="date",
        how="left",
    )
    samsung_dates = set(samsung_steps["date"])
    daily["steps_source"] = "missing"
    daily.loc[daily["steps"].notna(), "steps_source"] = "estimated backfill"
    daily.loc[daily["date"].isin(samsung_dates), "steps_source"] = "Samsung"

    sleep_one = sleep.sort_values("date").drop_duplicates("date", keep="last")
    daily = daily.merge(
        sleep_one[["date", "sleep_hours"]],
        on="date",
        how="left",
    )
    daily["sleep_observed"] = daily["sleep_hours"].notna()
    daily["sleep_plausible_2_to_11h"] = daily["sleep_hours"].between(
        VALID_SLEEP_MIN_HOURS,
        VALID_SLEEP_MAX_HOURS,
    )

    strength_dates = set(strength["date"])
    daily["strength_session"] = daily["date"].isin(strength_dates)

    tirz = medicine[
        medicine["drug"].str.contains("tirzepatide|zepbound", case=False, na=False)
    ][["date", "dose_mg", "clicks"]].rename(
        columns={
            "dose_mg": "injected_tirzepatide_mg",
            "clicks": "injection_clicks",
        }
    )
    daily = daily.merge(tirz, on="date", how="left")
    daily["tirzepatide_injection_observed"] = daily[
        "injected_tirzepatide_mg"
    ].notna()
    daily["tirzepatide_era"] = daily["date"] >= medicine.loc[
        medicine["drug"].str.contains(
            "tirzepatide|zepbound",
            case=False,
            na=False,
        ),
        "date",
    ].min()

    daily["travel"] = daily["date"].isin(parse_travel_dates())
    daily["temperature_readings"] = daily["date"].map(
        temperature.groupby("date").size()
    ).fillna(0).astype(int)

    daily.to_csv(
        ROOT / "analysis" / "EG_observation_provenance_daily.csv",
        index=False,
    )

    print("OBSERVATION / PROVENANCE TABLE")
    print("=" * 108)
    print(f"Rows: {len(daily):,}, {daily['date'].min().date()} to {daily['date'].max().date()}")
    print("\nWeight proximity")
    print(daily["weight_observation_tier"].value_counts().to_string())
    print("\nComposition proximity")
    print(daily["composition_observation_tier"].value_counts().to_string())
    print("\nStep source")
    print(daily["steps_source"].value_counts().to_string())
    print(
        f"\nSleep observed: {daily['sleep_observed'].sum():,}; "
        f"plausible 2–11 h: {daily['sleep_plausible_2_to_11h'].sum():,}"
    )
    print(
        f"Direct composition days: {daily['composition_observed'].sum():,}; "
        f"direct RMR days: {daily['rmr_observed'].sum():,}; "
        f"direct injection days inside intake range: "
        f"{daily['tirzepatide_injection_observed'].sum():,}"
    )


if __name__ == "__main__":
    main()
