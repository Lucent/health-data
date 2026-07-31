"""Audit named diet epochs using only direct intake, food, and scale observations.

This is intentionally independent of P1/P2/P3/P4, modeled body composition,
modeled expenditure, and modeled medication exposure. It answers narrower,
defensible questions:

1. What was logged before, during, and after each annotated calendar epoch?
2. Which annotations overlap, and how many days are exclusive to one label?
3. For focal natural experiments, how many days satisfy a transparent,
   food-log or macro-based confirmation rule?
4. What happened between exact observed scale measurements, without filling
   missing dates or interpreting short-term scale changes as fat change?

Run without arguments:
    uv run --with pandas python analysis/EJ_diet_epoch_direct_audit.py
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
COMPARISON_WINDOW_CAP_DAYS = 28
BINGE_CALORIES = 3000
LOW_INTAKE_CALORIES = 1200
POTATO_MIN_LOGGED_CALORIES = 500

FOCAL_LABELS = [
    "travel_binge_1",
    "high_protein_experiment",
    "keto_phase_1",
    "keto_phase_2",
    "enlightened_era",
    "travel_binge_2",
    "post_travel_restriction",
    "weekend_fasting",
    "potato_diet_1",
    "potato_diet_2",
    "low_protein_experiment",
    "potato_diet_3",
    "potato_diet_4",
]

OUTPUT_WINDOW = ROOT / "analysis" / "EJ_epoch_window_summary.csv"
OUTPUT_CONTRAST = ROOT / "analysis" / "EJ_epoch_contrasts.csv"
OUTPUT_OVERLAP = ROOT / "analysis" / "EJ_epoch_overlaps.csv"
OUTPUT_EXPOSURE = ROOT / "analysis" / "EJ_confirmed_exposure_summary.csv"


def load_sources() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load only minimally transformed observation tables."""
    daily = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    foods = pd.read_csv(ROOT / "intake" / "intake_foods.csv", parse_dates=["date"])
    weights = pd.read_csv(ROOT / "weight" / "weight.csv", parse_dates=["date"])
    epochs = pd.read_csv(ROOT / "intake" / "diet_epochs.csv", parse_dates=["start", "end"])

    daily = daily.sort_values("date").reset_index(drop=True)
    foods = foods.sort_values(["date", "meal", "food"]).reset_index(drop=True)
    weights = weights.sort_values("date").reset_index(drop=True)
    epochs = epochs.sort_values(["start", "end", "label"]).reset_index(drop=True)

    if daily["date"].duplicated().any():
        raise ValueError("intake_daily.csv must have one row per date")
    if weights["date"].duplicated().any():
        raise ValueError("weight.csv must have one row per date")
    if epochs["label"].duplicated().any():
        raise ValueError("diet epoch labels must be unique")
    if (epochs["end"] < epochs["start"]).any():
        raise ValueError("diet epoch end precedes start")

    numeric_cols = [
        "calories",
        "carbs_g",
        "fat_g",
        "protein_g",
        "fiber_g",
        "sugars_g",
        "sodium_mg",
    ]
    for col in numeric_cols:
        daily[col] = pd.to_numeric(daily[col], errors="coerce")
    foods["calories"] = pd.to_numeric(foods["calories"], errors="coerce").fillna(0)
    weights["weight_lbs"] = pd.to_numeric(weights["weight_lbs"], errors="coerce")
    weights = weights.dropna(subset=["weight_lbs"])

    daily["carb_pct_cal"] = 400 * daily["carbs_g"] / daily["calories"]
    daily["fat_pct_cal"] = 900 * daily["fat_g"] / daily["calories"]
    daily["protein_pct_cal"] = 400 * daily["protein_g"] / daily["calories"]
    daily["high_intake_day"] = daily["calories"] >= BINGE_CALORIES
    daily["low_intake_day"] = daily["calories"] <= LOW_INTAKE_CALORIES

    return daily, foods, weights, epochs


def safe_round(value: float, digits: int = 1) -> float | None:
    """Round a scalar while preserving missingness as an empty CSV cell."""
    if pd.isna(value):
        return None
    return round(float(value), digits)


def summarize_observed_window(
    daily: pd.DataFrame,
    weights: pd.DataFrame,
    start: pd.Timestamp,
    end: pd.Timestamp,
) -> dict[str, object]:
    """Summarize exact calendar rows and exact observed scale measurements."""
    expected_days = max((end - start).days + 1, 0)
    window = daily[daily["date"].between(start, end)].copy()
    scale = weights[weights["date"].between(start, end)].copy()

    row: dict[str, object] = {
        "window_start": start.strftime("%Y-%m-%d"),
        "window_end": end.strftime("%Y-%m-%d"),
        "expected_calendar_days": expected_days,
        "intake_days": int(len(window)),
        "intake_completeness": safe_round(len(window) / expected_days if expected_days else None, 4),
        "mean_calories": safe_round(window["calories"].mean()),
        "median_calories": safe_round(window["calories"].median()),
        "mean_carbs_g": safe_round(window["carbs_g"].mean()),
        "mean_fat_g": safe_round(window["fat_g"].mean()),
        "mean_protein_g": safe_round(window["protein_g"].mean()),
        "mean_fiber_g": safe_round(window["fiber_g"].mean()),
        "mean_carb_pct_cal": safe_round(window["carb_pct_cal"].mean()),
        "mean_fat_pct_cal": safe_round(window["fat_pct_cal"].mean()),
        "mean_protein_pct_cal": safe_round(window["protein_pct_cal"].mean()),
        "high_intake_days": int(window["high_intake_day"].sum()) if len(window) else 0,
        "high_intake_rate": safe_round(window["high_intake_day"].mean(), 4),
        "low_intake_days": int(window["low_intake_day"].sum()) if len(window) else 0,
        "low_intake_rate": safe_round(window["low_intake_day"].mean(), 4),
        "observed_weight_days": int(len(scale)),
        "first_observed_weight_date": None,
        "first_observed_weight_lbs": None,
        "last_observed_weight_date": None,
        "last_observed_weight_lbs": None,
        "observed_weight_span_days": None,
        "first_to_last_observed_weight_change_lbs": None,
    }
    if len(scale):
        first = scale.iloc[0]
        last = scale.iloc[-1]
        row.update(
            {
                "first_observed_weight_date": first["date"].strftime("%Y-%m-%d"),
                "first_observed_weight_lbs": safe_round(first["weight_lbs"], 2),
                "last_observed_weight_date": last["date"].strftime("%Y-%m-%d"),
                "last_observed_weight_lbs": safe_round(last["weight_lbs"], 2),
                "observed_weight_span_days": int((last["date"] - first["date"]).days),
                "first_to_last_observed_weight_change_lbs": (
                    safe_round(last["weight_lbs"] - first["weight_lbs"], 2) if len(scale) >= 2 else None
                ),
            }
        )
    return row


def build_window_summary(
    daily: pd.DataFrame,
    weights: pd.DataFrame,
    epochs: pd.DataFrame,
) -> pd.DataFrame:
    """Create capped pre/full-epoch/capped-post summaries for every label."""
    rows: list[dict[str, object]] = []
    for epoch in epochs.itertuples(index=False):
        duration = (epoch.end - epoch.start).days + 1
        flank_days = min(duration, COMPARISON_WINDOW_CAP_DAYS)
        phases = [
            ("pre", epoch.start - pd.Timedelta(days=flank_days), epoch.start - pd.Timedelta(days=1)),
            ("epoch", epoch.start, epoch.end),
            ("post", epoch.end + pd.Timedelta(days=1), epoch.end + pd.Timedelta(days=flank_days)),
        ]
        for phase, start, end in phases:
            row = {
                "label": epoch.label,
                "phase": phase,
                "epoch_start": epoch.start.strftime("%Y-%m-%d"),
                "epoch_end": epoch.end.strftime("%Y-%m-%d"),
                "epoch_calendar_days": duration,
                "flank_window_days": flank_days,
            }
            row.update(summarize_observed_window(daily, weights, start, end))
            rows.append(row)
    return pd.DataFrame(rows)


def build_overlap_summary(epochs: pd.DataFrame) -> pd.DataFrame:
    """List every pair of annotations sharing at least one calendar day."""
    rows: list[dict[str, object]] = []
    for left, right in combinations(epochs.itertuples(index=False), 2):
        start = max(left.start, right.start)
        end = min(left.end, right.end)
        if end < start:
            continue
        rows.append(
            {
                "label_1": left.label,
                "label_2": right.label,
                "overlap_start": start.strftime("%Y-%m-%d"),
                "overlap_end": end.strftime("%Y-%m-%d"),
                "overlap_calendar_days": int((end - start).days + 1),
            }
        )
    return pd.DataFrame(rows)


def overlap_context(
    epoch: pd.Series,
    epochs: pd.DataFrame,
    daily: pd.DataFrame,
) -> dict[str, object]:
    """Quantify other labels active within an epoch and summarize exclusive days."""
    target_days = pd.date_range(epoch["start"], epoch["end"], freq="D")
    other_labels: list[str] = []
    overlapped_dates: set[pd.Timestamp] = set()
    for other in epochs.itertuples(index=False):
        if other.label == epoch["label"]:
            continue
        start = max(epoch["start"], other.start)
        end = min(epoch["end"], other.end)
        if end >= start:
            other_labels.append(other.label)
            overlapped_dates.update(pd.date_range(start, end, freq="D"))

    exclusive_dates = set(target_days) - overlapped_dates
    exclusive = daily[daily["date"].isin(exclusive_dates)]
    return {
        "overlapping_labels": ";".join(sorted(other_labels)),
        "overlapped_calendar_days": len(overlapped_dates),
        "exclusive_calendar_days": len(exclusive_dates),
        "exclusive_intake_days": int(len(exclusive)),
        "exclusive_mean_calories": safe_round(exclusive["calories"].mean()),
        "exclusive_mean_protein_g": safe_round(exclusive["protein_g"].mean()),
        "exclusive_mean_carb_pct_cal": safe_round(exclusive["carb_pct_cal"].mean()),
        "exclusive_mean_fat_pct_cal": safe_round(exclusive["fat_pct_cal"].mean()),
    }


def build_contrasts(
    windows: pd.DataFrame,
    epochs: pd.DataFrame,
    daily: pd.DataFrame,
) -> pd.DataFrame:
    """Contrast during-epoch observations against transparent calendar flanks."""
    rows: list[dict[str, object]] = []
    measures = [
        "mean_calories",
        "mean_protein_g",
        "mean_fiber_g",
        "mean_carb_pct_cal",
        "mean_fat_pct_cal",
        "mean_protein_pct_cal",
        "high_intake_rate",
        "low_intake_rate",
    ]
    for epoch in epochs.itertuples(index=False):
        group = windows[windows["label"] == epoch.label].set_index("phase")
        row: dict[str, object] = {
            "label": epoch.label,
            "epoch_start": epoch.start.strftime("%Y-%m-%d"),
            "epoch_end": epoch.end.strftime("%Y-%m-%d"),
            "epoch_calendar_days": int((epoch.end - epoch.start).days + 1),
        }
        for measure in measures:
            pre = group.at["pre", measure]
            during = group.at["epoch", measure]
            post = group.at["post", measure]
            row[f"pre_{measure}"] = pre
            row[f"epoch_{measure}"] = during
            row[f"post_{measure}"] = post
            row[f"epoch_minus_pre_{measure}"] = safe_round(during - pre, 4)
            row[f"post_minus_epoch_{measure}"] = safe_round(post - during, 4)

        epoch_weight = group.loc["epoch"]
        row.update(
            {
                "observed_weight_days_during": epoch_weight["observed_weight_days"],
                "first_observed_weight_date_during": epoch_weight["first_observed_weight_date"],
                "first_observed_weight_lbs_during": epoch_weight["first_observed_weight_lbs"],
                "last_observed_weight_date_during": epoch_weight["last_observed_weight_date"],
                "last_observed_weight_lbs_during": epoch_weight["last_observed_weight_lbs"],
                "observed_weight_span_days_during": epoch_weight["observed_weight_span_days"],
                "first_to_last_observed_weight_change_lbs_during": epoch_weight[
                    "first_to_last_observed_weight_change_lbs"
                ],
            }
        )
        epoch_series = pd.Series(
            {"label": epoch.label, "start": epoch.start, "end": epoch.end},
            dtype=object,
        )
        row.update(overlap_context(epoch_series, epochs, daily))
        rows.append(row)
    return pd.DataFrame(rows)


def food_calories_by_date(foods: pd.DataFrame, pattern: str, column: str) -> pd.DataFrame:
    """Aggregate matching logged food calories by calendar date."""
    matches = foods[foods["food"].str.contains(pattern, case=False, na=False, regex=True)]
    return matches.groupby("date", as_index=False)["calories"].sum().rename(columns={"calories": column})


def confirmation_rule(
    label: str,
    frame: pd.DataFrame,
) -> tuple[pd.Series, str]:
    """Return a prespecified, auditable exposure-confirmation mask."""
    if label.startswith("potato_diet_"):
        return (
            frame["potato_food_calories"] >= POTATO_MIN_LOGGED_CALORIES,
            f"logged potato/yukon/russet foods >= {POTATO_MIN_LOGGED_CALORIES} kcal",
        )
    if label == "weekend_fasting":
        return (
            (frame["calories"] < 200) & (frame["date"].dt.dayofweek >= 5),
            "Saturday or Sunday with <200 logged kcal",
        )
    if label == "low_protein_experiment":
        return frame["protein_pct_cal"] <= 10, "protein <=10% of logged calories"
    if label == "high_protein_experiment":
        return frame["protein_pct_cal"] >= 20, "protein >=20% of logged calories"
    if label.startswith("keto_phase_"):
        return (
            (frame["carb_pct_cal"] <= 25) & (frame["fat_pct_cal"] >= 50),
            "carbohydrate <=25% and fat >=50% of logged calories",
        )
    if label == "enlightened_era":
        return frame["enlightened_food_calories"] > 0, "any logged Enlightened-brand calories"
    if label.startswith("travel_binge_"):
        return frame["calories"] >= BINGE_CALORIES, f"logged intake >= {BINGE_CALORIES} kcal"
    if label == "post_travel_restriction":
        return frame["calories"] <= 1800, "logged intake <=1800 kcal"
    raise KeyError(f"No confirmation rule for {label}")


def build_confirmed_exposure_summary(
    daily: pd.DataFrame,
    foods: pd.DataFrame,
    weights: pd.DataFrame,
    epochs: pd.DataFrame,
) -> pd.DataFrame:
    """Compare confirmed versus unconfirmed days inside focal annotations."""
    potato = food_calories_by_date(foods, r"potato|yukon|russet", "potato_food_calories")
    enlightened = food_calories_by_date(foods, r"enlightened", "enlightened_food_calories")
    weight_map = weights.set_index("date")["weight_lbs"]

    augmented = daily.merge(potato, on="date", how="left").merge(enlightened, on="date", how="left")
    augmented[["potato_food_calories", "enlightened_food_calories"]] = augmented[
        ["potato_food_calories", "enlightened_food_calories"]
    ].fillna(0)
    augmented["observed_weight_lbs"] = augmented["date"].map(weight_map)
    augmented["next_date"] = augmented["date"] + pd.Timedelta(days=1)
    augmented["next_day_observed_weight_lbs"] = augmented["next_date"].map(weight_map)
    augmented["next_calendar_day_weight_change_lbs"] = (
        augmented["next_day_observed_weight_lbs"] - augmented["observed_weight_lbs"]
    )

    rows: list[dict[str, object]] = []
    epoch_index = epochs.set_index("label")
    for label in FOCAL_LABELS:
        epoch = epoch_index.loc[label]
        frame = augmented[augmented["date"].between(epoch["start"], epoch["end"])].copy()
        confirmed, criterion = confirmation_rule(label, frame)
        frame["confirmed"] = confirmed.fillna(False)

        for status, subset in frame.groupby("confirmed", dropna=False):
            rows.append(
                {
                    "label": label,
                    "criterion": criterion,
                    "group": "confirmed_exposure" if bool(status) else "other_labeled_days",
                    "days": int(len(subset)),
                    "share_of_labeled_days": safe_round(len(subset) / len(frame) if len(frame) else None, 4),
                    "mean_calories": safe_round(subset["calories"].mean()),
                    "median_calories": safe_round(subset["calories"].median()),
                    "mean_carbs_g": safe_round(subset["carbs_g"].mean()),
                    "mean_fat_g": safe_round(subset["fat_g"].mean()),
                    "mean_protein_g": safe_round(subset["protein_g"].mean()),
                    "mean_fiber_g": safe_round(subset["fiber_g"].mean()),
                    "mean_carb_pct_cal": safe_round(subset["carb_pct_cal"].mean()),
                    "mean_fat_pct_cal": safe_round(subset["fat_pct_cal"].mean()),
                    "mean_protein_pct_cal": safe_round(subset["protein_pct_cal"].mean()),
                    "high_intake_rate": safe_round(subset["high_intake_day"].mean(), 4),
                    "exact_next_day_weight_pairs": int(
                        subset["next_calendar_day_weight_change_lbs"].notna().sum()
                    ),
                    "mean_next_calendar_day_weight_change_lbs": safe_round(
                        subset["next_calendar_day_weight_change_lbs"].mean(), 3
                    ),
                    "median_next_calendar_day_weight_change_lbs": safe_round(
                        subset["next_calendar_day_weight_change_lbs"].median(), 3
                    ),
                    "total_potato_food_calories": safe_round(subset["potato_food_calories"].sum(), 0),
                    "total_enlightened_food_calories": safe_round(
                        subset["enlightened_food_calories"].sum(), 0
                    ),
                }
            )
        if frame["confirmed"].all():
            rows.append(
                {
                    "label": label,
                    "criterion": criterion,
                    "group": "other_labeled_days",
                    "days": 0,
                    "share_of_labeled_days": 0.0,
                }
            )
    return pd.DataFrame(rows)


def print_key_results(
    contrasts: pd.DataFrame,
    exposure: pd.DataFrame,
    overlaps: pd.DataFrame,
) -> None:
    """Print a compact, reproducible console synopsis."""
    focal = contrasts[contrasts["label"].isin(FOCAL_LABELS)].copy()
    cols = [
        "label",
        "pre_mean_calories",
        "epoch_mean_calories",
        "epoch_minus_pre_mean_calories",
        "post_mean_calories",
        "first_to_last_observed_weight_change_lbs_during",
        "observed_weight_days_during",
        "overlapping_labels",
    ]
    print("\n=== Direct diet-epoch audit ===")
    print(focal[cols].to_string(index=False))

    confirmed = exposure[exposure["group"] == "confirmed_exposure"]
    exp_cols = [
        "label",
        "days",
        "share_of_labeled_days",
        "mean_calories",
        "mean_protein_g",
        "exact_next_day_weight_pairs",
        "mean_next_calendar_day_weight_change_lbs",
    ]
    print("\nConfirmed exposure-day summaries:")
    print(confirmed[exp_cols].to_string(index=False))

    print("\nOverlapping annotations:")
    print(overlaps.to_string(index=False) if len(overlaps) else "none")
    print(f"\nWrote {OUTPUT_WINDOW.relative_to(ROOT)}")
    print(f"Wrote {OUTPUT_CONTRAST.relative_to(ROOT)}")
    print(f"Wrote {OUTPUT_OVERLAP.relative_to(ROOT)}")
    print(f"Wrote {OUTPUT_EXPOSURE.relative_to(ROOT)}")


def main() -> None:
    daily, foods, weights, epochs = load_sources()
    windows = build_window_summary(daily, weights, epochs)
    overlaps = build_overlap_summary(epochs)
    contrasts = build_contrasts(windows, epochs, daily)
    exposure = build_confirmed_exposure_summary(daily, foods, weights, epochs)

    windows.to_csv(OUTPUT_WINDOW, index=False)
    contrasts.to_csv(OUTPUT_CONTRAST, index=False)
    overlaps.to_csv(OUTPUT_OVERLAP, index=False)
    exposure.to_csv(OUTPUT_EXPOSURE, index=False)
    print_key_results(contrasts, exposure, overlaps)


if __name__ == "__main__":
    main()
