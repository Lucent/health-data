#!/usr/bin/env python3
"""Summarize Samsung exercise types and daylight walking regimes.

This uses session-level exercise exports to separate:
- confirmed exercise type codes
- strict paired daylight walks
- borderline paired daylight walks
- single daylight walks

The goal is to test whether deliberate sunlight walk structure carries signal
that bare daily step counts do not.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
EXERCISES_PATH = ROOT / "steps-sleep" / "exercises_samsung.csv"
STEPS_PATH = ROOT / "steps-sleep" / "steps.csv"
INTAKE_PATH = ROOT / "intake" / "intake_daily.csv"
KALMAN_PATH = ROOT / "analysis" / "P4_kalman_daily.csv"
COMPOSITION_PATH = ROOT / "analysis" / "P2_daily_composition.csv"

TYPE_SUMMARY_PATH = ROOT / "analysis" / "V_exercise_type_summary.csv"
WALK_DAYS_PATH = ROOT / "analysis" / "V_daylight_walk_regime_days.csv"
WALK_SUMMARY_PATH = ROOT / "analysis" / "V_daylight_walk_regime_summary.csv"
WALK_CONTRAST_PATH = ROOT / "analysis" / "V_daylight_walk_regime_contrast.csv"


# Walk regime detection parameters
MIN_YEAR = 2020
HOUR_START = 12       # noon
HOUR_END = 19         # 7pm
MIN_SESSION_MIN = 20  # minimum per-leg duration
MAX_SESSION_MIN = 45
PAIRED_BREAK_MIN = (5, 40)   # break between legs
SINGLE_MIN_DURATION = 20

TYPE_LABELS = {
    1001: ("walking", "high"),
    1002: ("running", "high"),
    11007: ("bike", "high"),
    13001: ("hiking", "high"),
    15003: ("indoor_bike", "high"),
}


def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    exercises = pd.read_csv(EXERCISES_PATH, parse_dates=["start_time", "end_time"])
    steps = pd.read_csv(STEPS_PATH, parse_dates=["date"])
    return exercises, steps


def build_type_summary(exercises: pd.DataFrame) -> pd.DataFrame:
    df = exercises.copy()
    for col in ["exercise_type", "count", "distance", "calorie", "duration_min"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["pace_m_per_min"] = df["distance"] / df["duration_min"]
    df["cal_per_min"] = df["calorie"] / df["duration_min"]
    df["year"] = df["start_time"].dt.year

    rows = []
    for exercise_type, group in df.groupby("exercise_type", dropna=False):
        if pd.isna(exercise_type):
            label, confidence = ("unknown", "low")
            exercise_type_value = ""
        else:
            label, confidence = TYPE_LABELS.get(int(exercise_type), ("unknown", "low"))
            exercise_type_value = int(exercise_type)
        rows.append(
            {
                "exercise_type": exercise_type_value,
                "exercise_label": label,
                "label_confidence": confidence,
                "sessions": int(len(group)),
                "first_date": group["start_time"].min().date().isoformat(),
                "last_date": group["start_time"].max().date().isoformat(),
                "recent_sessions_2024_plus": int((group["year"] >= 2024).sum()),
                "source_types": ",".join(sorted(set(group["source_type"].dropna().astype(str)))),
                "mean_duration_min": round(group["duration_min"].mean(), 2),
                "median_duration_min": round(group["duration_min"].median(), 2),
                "mean_distance_m": round(group["distance"].mean(), 2),
                "mean_steps": round(group["count"].mean(), 2),
                "nonnull_step_fraction": round(group["count"].notna().mean(), 4),
                "mean_pace_m_per_min": round(group["pace_m_per_min"].mean(), 2),
                "mean_cal_per_min": round(group["cal_per_min"].mean(), 2),
            }
        )
    out = pd.DataFrame(rows).sort_values(["sessions", "exercise_type"], ascending=[False, True])
    out.to_csv(TYPE_SUMMARY_PATH, index=False)
    return out


def classify_daylight_walk_regimes(exercises: pd.DataFrame) -> pd.DataFrame:
    walk = exercises.copy()
    for col in ["exercise_type", "count", "distance", "duration_min"]:
        walk[col] = pd.to_numeric(walk[col], errors="coerce")
    walk = walk[walk["exercise_type"] == 1001].copy()
    walk["date"] = walk["start_time"].dt.normalize()
    walk["year"] = walk["start_time"].dt.year
    walk["month"] = walk["start_time"].dt.month
    walk["hour"] = walk["start_time"].dt.hour

    base = (
        (walk["year"] >= MIN_YEAR)
        & (walk["hour"] >= HOUR_START)
        & (walk["hour"] < HOUR_END)
    )
    candidates = walk[base].sort_values(["date", "start_time"]).copy()

    rows = []
    for date, group in candidates.groupby("date"):
        group = group.reset_index(drop=True)
        regime = None
        chosen = None

        for i in range(len(group) - 1):
            break_min = (group.loc[i + 1, "start_time"] - group.loc[i, "end_time"]).total_seconds() / 60
            if (
                group.loc[i, "duration_min"] >= MIN_SESSION_MIN
                and group.loc[i, "duration_min"] <= MAX_SESSION_MIN
                and group.loc[i + 1, "duration_min"] >= MIN_SESSION_MIN
                and group.loc[i + 1, "duration_min"] <= MAX_SESSION_MIN
                and PAIRED_BREAK_MIN[0] <= break_min <= PAIRED_BREAK_MIN[1]
            ):
                regime = "paired_daylight_walk"
                chosen = {
                    "first_duration_min": round(group.loc[i, "duration_min"], 2),
                    "second_duration_min": round(group.loc[i + 1, "duration_min"], 2),
                    "break_min": round(break_min, 2),
                    "regime_duration_min": round(group.loc[i : i + 1, "duration_min"].sum(), 2),
                    "regime_distance_m": round(group.loc[i : i + 1, "distance"].sum(), 2),
                    "regime_steps": round(group.loc[i : i + 1, "count"].fillna(0).sum(), 2),
                    "source_types": ",".join(sorted(set(group.loc[i : i + 1, "source_type"].astype(str)))),
                }
                break

        if regime is None:
            for i in range(len(group) - 1):
                break_min = (group.loc[i + 1, "start_time"] - group.loc[i, "end_time"]).total_seconds() / 60
                if (
                    group.loc[i, "duration_min"] >= 20
                    and group.loc[i, "duration_min"] <= 45
                    and group.loc[i + 1, "duration_min"] >= 10
                    and group.loc[i + 1, "duration_min"] <= 45
                    and 0 <= break_min <= 90
                ):
                    regime = "borderline_paired_daylight_walk"
                    chosen = {
                        "first_duration_min": round(group.loc[i, "duration_min"], 2),
                        "second_duration_min": round(group.loc[i + 1, "duration_min"], 2),
                        "break_min": round(break_min, 2),
                        "regime_duration_min": round(group.loc[i : i + 1, "duration_min"].sum(), 2),
                        "regime_distance_m": round(group.loc[i : i + 1, "distance"].sum(), 2),
                        "regime_steps": round(group.loc[i : i + 1, "count"].fillna(0).sum(), 2),
                        "source_types": ",".join(sorted(set(group.loc[i : i + 1, "source_type"].astype(str)))),
                    }
                    break

        if regime is None:
            singles = group[(group["duration_min"] >= SINGLE_MIN_DURATION) & (group["duration_min"] <= MAX_SESSION_MIN)]
            if not singles.empty:
                first = singles.iloc[0]
                regime = "single_daylight_walk"
                chosen = {
                    "first_duration_min": round(first["duration_min"], 2),
                    "second_duration_min": "",
                    "break_min": "",
                    "regime_duration_min": round(first["duration_min"], 2),
                    "regime_distance_m": round(first["distance"], 2),
                    "regime_steps": round(0 if pd.isna(first["count"]) else first["count"], 2),
                    "source_types": str(first["source_type"]),
                }

        if regime is not None:
            rows.append(
                {
                    "date": date.date().isoformat(),
                    "year": int(group.loc[0, "year"]),
                    "walk_regime": regime,
                    "first_start": group.loc[0, "start_time"].strftime("%Y-%m-%d %H:%M:%S"),
                    **chosen,
                }
            )

    out = pd.DataFrame(rows).sort_values(["date", "walk_regime"])
    out.to_csv(WALK_DAYS_PATH, index=False)
    return out


def summarize_walk_regimes(walk_days: pd.DataFrame, steps: pd.DataFrame) -> pd.DataFrame:
    steps = steps.copy()
    steps["date"] = steps["date"].dt.strftime("%Y-%m-%d")
    steps["steps"] = pd.to_numeric(steps["steps"], errors="coerce")

    merged = walk_days.merge(steps[["date", "steps"]], on="date", how="left")
    rows = []
    for regime, group in merged.groupby("walk_regime"):
        rows.append(
            {
                "walk_regime": regime,
                "days": int(len(group)),
                "mean_regime_duration_min": round(group["regime_duration_min"].mean(), 2),
                "median_regime_duration_min": round(group["regime_duration_min"].median(), 2),
                "mean_regime_distance_m": round(group["regime_distance_m"].mean(), 2),
                "median_regime_distance_m": round(group["regime_distance_m"].median(), 2),
                "mean_regime_steps": round(group["regime_steps"].mean(), 2),
                "median_regime_steps": round(group["regime_steps"].median(), 2),
                "mean_total_day_steps": round(group["steps"].mean(), 2),
                "median_total_day_steps": round(group["steps"].median(), 2),
            }
        )
    out = pd.DataFrame(rows).sort_values("walk_regime")
    out.to_csv(WALK_SUMMARY_PATH, index=False)
    return out


def build_walk_contrast(walk_days: pd.DataFrame, steps: pd.DataFrame) -> pd.DataFrame:
    intake = pd.read_csv(INTAKE_PATH, parse_dates=["date"])[["date", "calories"]]
    kalman = pd.read_csv(KALMAN_PATH, parse_dates=["date"])[["date", "tdee"]]
    composition = pd.read_csv(COMPOSITION_PATH, parse_dates=["date"])[["date", "expected_rmr"]]

    daily = (
        steps.merge(intake, on="date", how="left")
        .merge(kalman, on="date", how="left")
        .merge(composition, on="date", how="left")
    )
    daily["date"] = daily["date"].dt.strftime("%Y-%m-%d")
    daily["steps"] = pd.to_numeric(daily["steps"], errors="coerce")
    daily["tdee_rmr_ratio"] = daily["tdee"] / daily["expected_rmr"]
    daily = daily.merge(walk_days[["date", "walk_regime"]], on="date", how="left")

    for horizon in [14]:
        daily[f"future{horizon}_calories"] = daily["calories"].shift(-1).rolling(
            horizon, min_periods=7
        ).mean()
        daily[f"future{horizon}_ratio"] = daily["tdee_rmr_ratio"].shift(-1).rolling(
            horizon, min_periods=7
        ).mean()

    rows = []
    for regime in [
        "paired_daylight_walk",
        "borderline_paired_daylight_walk",
        "single_daylight_walk",
    ]:
        group = daily[daily["walk_regime"] == regime]
        rows.append(
            {
                "group": regime,
                "days": int(len(group)),
                "mean_steps": round(group["steps"].mean(), 2),
                "mean_same_day_calories": round(group["calories"].mean(), 2),
                "mean_future14_calories": round(group["future14_calories"].mean(), 2),
                "mean_future14_tdee_rmr_ratio": round(group["future14_ratio"].mean(), 4),
                "mean_walk_regime_duration_min": round(
                    walk_days.loc[walk_days["walk_regime"] == regime, "regime_duration_min"].mean(), 2
                ),
            }
        )

    for regime in ["paired_daylight_walk", "single_daylight_walk"]:
        target = daily[daily["walk_regime"] == regime][["steps", "future14_calories", "future14_ratio"]].dropna()
        others = (
            daily[daily["walk_regime"] != regime][["steps", "future14_calories", "future14_ratio"]]
            .dropna()
            .sort_values("steps")
            .reset_index(drop=True)
        )
        matched_rows = []
        for _, row in target.iterrows():
            if others.empty:
                break
            idx = (others["steps"] - row["steps"]).abs().idxmin()
            match = others.loc[idx]
            matched_rows.append(
                {
                    "step_distance_to_match": abs(row["steps"] - match["steps"]),
                    "future14_calories_diff": row["future14_calories"] - match["future14_calories"],
                    "future14_tdee_rmr_ratio_diff": row["future14_ratio"] - match["future14_ratio"],
                }
            )
            others = others.drop(idx).reset_index(drop=True)
        matched = pd.DataFrame(matched_rows)
        rows.append(
            {
                "group": f"{regime}_step_matched_minus_other_days",
                "days": int(len(matched)),
                "mean_steps": round(target["steps"].mean(), 2),
                "mean_same_day_calories": "",
                "mean_future14_calories": round(matched["future14_calories_diff"].mean(), 2),
                "mean_future14_tdee_rmr_ratio": round(matched["future14_tdee_rmr_ratio_diff"].mean(), 4),
                "mean_walk_regime_duration_min": round(
                    walk_days.loc[walk_days["walk_regime"] == regime, "regime_duration_min"].mean(), 2
                ),
                "step_distance_to_match": round(matched["step_distance_to_match"].mean(), 2),
            }
        )

    out = pd.DataFrame(rows)
    out.to_csv(WALK_CONTRAST_PATH, index=False)
    return out


def main() -> None:
    exercises, steps = load_data()
    type_summary = build_type_summary(exercises)
    walk_days = classify_daylight_walk_regimes(exercises)
    walk_summary = summarize_walk_regimes(walk_days, steps)
    walk_contrast = build_walk_contrast(walk_days, steps)

    walking = type_summary[type_summary["exercise_label"] == "walking"].iloc[0]
    paired = walk_summary[walk_summary["walk_regime"] == "paired_daylight_walk"].iloc[0]
    single = walk_summary[walk_summary["walk_regime"] == "single_daylight_walk"].iloc[0]
    paired_match = walk_contrast[
        walk_contrast["group"] == "paired_daylight_walk_step_matched_minus_other_days"
    ].iloc[0]
    single_match = walk_contrast[
        walk_contrast["group"] == "single_daylight_walk_step_matched_minus_other_days"
    ].iloc[0]

    print(
        f"Type 1001 walking sessions: {int(walking['sessions'])}, "
        f"mean {walking['mean_duration_min']:.2f} min, mean {walking['mean_distance_m']:.0f} m, "
        f"mean {walking['mean_steps']:.0f} steps."
    )
    print(
        f"Paired daylight walks: {int(paired['days'])} days, "
        f"median {paired['median_regime_duration_min']:.2f} min, "
        f"median {paired['median_regime_steps']:.0f} regime steps, "
        f"median {paired['median_total_day_steps']:.0f} total day steps."
    )
    print(
        f"Single daylight walks: {int(single['days'])} days, "
        f"median {single['median_regime_duration_min']:.2f} min, "
        f"median {single['median_regime_steps']:.0f} regime steps, "
        f"median {single['median_total_day_steps']:.0f} total day steps."
    )
    print(
        f"Step-matched paired daylight walks vs other days: "
        f"future14 calories {paired_match['mean_future14_calories']:+.2f}, "
        f"future14 TDEE/RMR {paired_match['mean_future14_tdee_rmr_ratio']:+.4f}."
    )
    print(
        f"Step-matched single daylight walks vs other days: "
        f"future14 calories {single_match['mean_future14_calories']:+.2f}, "
        f"future14 TDEE/RMR {single_match['mean_future14_tdee_rmr_ratio']:+.4f}."
    )


if __name__ == "__main__":
    main()
