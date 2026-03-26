"""Test whether the 'dead zone' is really a phase-dependent expenditure pattern.

The original claim was a single intake band (about 2000-2500 kcal) where
metabolism adapts to maintain. This script checks a stronger alternative:
the intake->TDEE slope depends on trajectory phase, and the alleged dead zone
may be branch-dependent rather than universal.

Outputs:
    analysis/deadzone_phase_bin_summary.csv
    analysis/deadzone_phase_regression.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

TREND_WINDOW_DAYS = 90
TREND_THRESHOLD_LBS = 3.0
DEADZONE_LOW = 2000
DEADZONE_HIGH = 2500


def load_data():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])

    daily = intake.merge(kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="left")
    daily = daily.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    daily = daily.merge(tirz[["date", "effective_level"]], on="date", how="left")
    daily["effective_level"] = daily["effective_level"].fillna(0)
    return daily.sort_values("date").reset_index(drop=True)


def build_frame(daily):
    d = daily[daily["effective_level"] == 0].copy()
    d["fat_delta_window"] = d["fat_mass_lbs"].diff(TREND_WINDOW_DAYS)
    d["phase"] = np.where(
        d["fat_delta_window"] <= -TREND_THRESHOLD_LBS,
        "falling",
        np.where(d["fat_delta_window"] >= TREND_THRESHOLD_LBS, "rising", "stable"),
    )
    d["tdee_minus_rmr"] = d["tdee"] - d["expected_rmr"]
    d["is_deadzone"] = (d["calories"] >= DEADZONE_LOW) & (d["calories"] < DEADZONE_HIGH)
    d["intake_bin"] = pd.cut(
        d["calories"],
        bins=[0, 1600, 2000, 2500, 3000, 10000],
        labels=["<1600", "1600-1999", "2000-2499", "2500-2999", "3000+"],
        right=False,
    )
    return d.dropna(subset=["fat_mass_lbs", "tdee", "expected_rmr", "calories"])


def bin_summary(daily):
    rows = []
    for phase in ["falling", "stable", "rising"]:
        grp = daily[daily["phase"] == phase]
        for intake_bin, sub in grp.groupby("intake_bin", observed=False):
            if len(sub) < 30:
                continue
            rows.append(
                {
                    "phase": phase,
                    "intake_bin": str(intake_bin),
                    "n_days": len(sub),
                    "mean_calories": round(sub["calories"].mean(), 1),
                    "mean_tdee": round(sub["tdee"].mean(), 1),
                    "mean_tdee_minus_rmr": round(sub["tdee_minus_rmr"].mean(), 1),
                    "mean_fat_mass_lbs": round(sub["fat_mass_lbs"].mean(), 1),
                }
            )
    return pd.DataFrame(rows)


def phase_regressions(daily):
    rows = []
    subsets = {
        "all": daily["calories"] >= 0,
        "deadzone": daily["is_deadzone"],
        "outside_deadzone": ~daily["is_deadzone"],
    }

    for subset_name, mask in subsets.items():
        subset = daily[mask]
        for phase in ["falling", "stable", "rising"]:
            grp = subset[subset["phase"] == phase]
            if len(grp) < 50:
                continue
            X = np.column_stack([np.ones(len(grp)), grp["fat_mass_lbs"].values, grp["calories"].values])
            coef = np.linalg.lstsq(X, grp["tdee"].values, rcond=None)[0]
            rows.append(
                {
                    "subset": subset_name,
                    "phase": phase,
                    "n_days": len(grp),
                    "coef_fat_mass": round(coef[1], 4),
                    "coef_calories": round(coef[2], 5),
                }
            )
    return pd.DataFrame(rows)


def save_outputs(bin_df, reg_df):
    bin_df.to_csv(ROOT / "analysis" / "I_deadzone_phase_bin_summary.csv", index=False)
    reg_df.to_csv(ROOT / "analysis" / "I_deadzone_phase_regression.csv", index=False)


def print_report(bin_df, reg_df):
    print("\n=== Dead Zone vs Phase-Dependent Expenditure ===")
    print("Mean TDEE by intake bin and phase:")
    for _, row in bin_df.iterrows():
        print(
            f"{row['phase']:>7}  {row['intake_bin']:>9}: TDEE {row['mean_tdee']:.0f}  "
            f"TDEE-RMR {row['mean_tdee_minus_rmr']:.0f}  FM {row['mean_fat_mass_lbs']:.1f}  n={int(row['n_days'])}"
        )

    print("\nRegression slopes (TDEE ~ fat_mass + calories):")
    for _, row in reg_df.iterrows():
        print(
            f"{row['subset']:>15}  {row['phase']:>7}: coef_calories={row['coef_calories']:+.5f}  "
            f"coef_fat_mass={row['coef_fat_mass']:+.4f}  n={int(row['n_days'])}"
        )


def main():
    daily = build_frame(load_data())
    bin_df = bin_summary(daily)
    reg_df = phase_regressions(daily)
    save_outputs(bin_df, reg_df)
    print_report(bin_df, reg_df)


if __name__ == "__main__":
    main()
