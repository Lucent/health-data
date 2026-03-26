"""Temperature by trajectory phase at matched fat mass.

This is a partial identifiability check rather than a definitive causal test.
Temperature starts in Dec 2023, so pre-tirzepatide overlap mostly covers
stable/rising high-fat-mass days, while on-drug overlap mostly covers falling
days during weight loss.

Outputs:
    analysis/temperature_daily.csv
    analysis/temperature_phase_overlap.csv
    analysis/temperature_phase_band_summary.csv
    analysis/temperature_phase_regression.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

TREND_WINDOW_DAYS = 90
TREND_THRESHOLD_LBS = 3.0


def load_data():
    temp = pd.read_csv(ROOT / "temperature" / "temperature.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    temp = temp[(temp["temp_f"] > 95) & (temp["temp_f"] < 101)].copy()
    temp["date"] = temp["date"].dt.floor("D")
    daily_temp = (
        temp.groupby("date")
        .agg(temp_mean=("temp_f", "mean"), temp_median=("temp_f", "median"), n_readings=("temp_f", "size"))
        .reset_index()
    )

    daily = daily_temp.merge(kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="left")
    daily = daily.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    daily = daily.merge(tirz[["date", "effective_level", "dose_mg"]], on="date", how="left")
    daily["effective_level"] = daily["effective_level"].fillna(0)
    daily["on_tirz"] = daily["effective_level"] > 0
    return daily.sort_values("date").reset_index(drop=True)


def build_frame(daily):
    d = daily.copy()
    d["fat_delta_window"] = d["fat_mass_lbs"].diff(TREND_WINDOW_DAYS)
    d["phase"] = np.where(
        d["fat_delta_window"] <= -TREND_THRESHOLD_LBS,
        "falling",
        np.where(d["fat_delta_window"] >= TREND_THRESHOLD_LBS, "rising", "stable"),
    )
    d["phase_falling"] = (d["phase"] == "falling").astype(int)
    d["phase_rising"] = (d["phase"] == "rising").astype(int)
    d["tdee_rmr_ratio"] = d["tdee"] / d["expected_rmr"]
    return d.dropna(subset=["temp_mean", "fat_mass_lbs", "tdee", "expected_rmr"])


def overlap_summary(daily):
    rows = []
    grouped = daily.groupby(["phase", "on_tirz"], dropna=False)
    for (phase, on_tirz), grp in grouped:
        rows.append(
            {
                "phase": phase,
                "on_tirz": bool(on_tirz),
                "n_days": len(grp),
                "mean_fat_mass_lbs": round(grp["fat_mass_lbs"].mean(), 2),
                "mean_temp_f": round(grp["temp_mean"].mean(), 4),
                "mean_tdee_rmr_ratio": round(grp["tdee_rmr_ratio"].mean(), 4),
            }
        )
    return pd.DataFrame(rows).sort_values(["on_tirz", "phase"])


def matched_band_summary(daily):
    rows = []
    for low, high in [(60, 70), (65, 75), (70, 80), (75, 85)]:
        band = daily[(daily["fat_mass_lbs"] >= low) & (daily["fat_mass_lbs"] < high)].copy()
        if band.empty:
            continue
        for (phase, on_tirz), grp in band.groupby(["phase", "on_tirz"], dropna=False):
            rows.append(
                {
                    "fat_band_lbs": f"{low}-{high}",
                    "phase": phase,
                    "on_tirz": bool(on_tirz),
                    "n_days": len(grp),
                    "mean_temp_f": round(grp["temp_mean"].mean(), 4),
                    "mean_tdee_rmr_ratio": round(grp["tdee_rmr_ratio"].mean(), 4),
                }
            )
    return pd.DataFrame(rows)


def regression_summary(daily):
    rows = []

    # Overall regression: temperature ~ fat_mass + on_tirz + phase + effective level
    X = np.column_stack(
        [
            np.ones(len(daily)),
            daily["fat_mass_lbs"].values,
            daily["on_tirz"].astype(int).values,
            daily["phase_falling"].values,
            daily["phase_rising"].values,
            daily["effective_level"].values,
        ]
    )
    coef = np.linalg.lstsq(X, daily["temp_mean"].values, rcond=None)[0]
    labels = ["intercept", "fat_mass_lbs", "on_tirz", "phase_falling", "phase_rising", "effective_level"]
    for label, value in zip(labels, coef):
        rows.append({"model": "overall", "term": label, "coef": round(value, 6)})

    # On-drug only
    on = daily[daily["on_tirz"]].copy()
    if len(on) >= 30:
        X = np.column_stack(
            [
                np.ones(len(on)),
                on["fat_mass_lbs"].values,
                on["phase_falling"].values,
                on["phase_rising"].values,
                on["effective_level"].values,
            ]
        )
        coef = np.linalg.lstsq(X, on["temp_mean"].values, rcond=None)[0]
        labels = ["intercept", "fat_mass_lbs", "phase_falling", "phase_rising", "effective_level"]
        for label, value in zip(labels, coef):
            rows.append({"model": "on_tirz_only", "term": label, "coef": round(value, 6)})

    # Pre-drug only
    pre = daily[~daily["on_tirz"]].copy()
    if len(pre) >= 30:
        X = np.column_stack(
            [
                np.ones(len(pre)),
                pre["fat_mass_lbs"].values,
                pre["phase_falling"].values,
                pre["phase_rising"].values,
            ]
        )
        coef = np.linalg.lstsq(X, pre["temp_mean"].values, rcond=None)[0]
        labels = ["intercept", "fat_mass_lbs", "phase_falling", "phase_rising"]
        for label, value in zip(labels, coef):
            rows.append({"model": "pre_tirz_only", "term": label, "coef": round(value, 6)})

    return pd.DataFrame(rows)


def save_outputs(temp_df, overlap_df, band_df, reg_df):
    temp_df.to_csv(ROOT / "analysis" / "temperature_daily.csv", index=False)
    overlap_df.to_csv(ROOT / "analysis" / "temperature_phase_overlap.csv", index=False)
    band_df.to_csv(ROOT / "analysis" / "temperature_phase_band_summary.csv", index=False)
    reg_df.to_csv(ROOT / "analysis" / "temperature_phase_regression.csv", index=False)


def print_report(overlap_df, band_df, reg_df):
    print("\n=== Temperature by Phase (Partial Identifiability) ===")
    print("Overlap by phase and tirzepatide status:")
    for _, row in overlap_df.iterrows():
        cohort = "on_tirz" if row["on_tirz"] else "pre_tirz"
        print(
            f"  {cohort:>8}  {row['phase']:>7}: n={int(row['n_days'])}  "
            f"temp={row['mean_temp_f']:.3f}F  FM={row['mean_fat_mass_lbs']:.1f}  "
            f"TDEE/RMR={row['mean_tdee_rmr_ratio']:.3f}"
        )

    print("\nMatched fat-mass bands:")
    for _, row in band_df.iterrows():
        cohort = "on_tirz" if row["on_tirz"] else "pre_tirz"
        print(
            f"  FM {row['fat_band_lbs']:>5}  {cohort:>8}  {row['phase']:>7}: "
            f"temp={row['mean_temp_f']:.3f}F  ratio={row['mean_tdee_rmr_ratio']:.3f}  n={int(row['n_days'])}"
        )

    print("\nRegression coefficients:")
    for _, row in reg_df.iterrows():
        print(f"  {row['model']:>12}  {row['term']:>14}: {row['coef']:+.6f}")

    print("\nInterpretation: temperature-by-phase is only partially identifiable because")
    print("pre-drug overlap has almost no falling-phase days, while on-drug overlap")
    print("is dominated by falling-phase weight loss. Treat coefficients as descriptive.")


def main():
    raw = load_data()
    daily = build_frame(raw)
    overlap_df = overlap_summary(daily)
    band_df = matched_band_summary(daily)
    reg_df = regression_summary(daily)
    save_outputs(raw, overlap_df, band_df, reg_df)
    print_report(overlap_df, band_df, reg_df)


if __name__ == "__main__":
    main()
