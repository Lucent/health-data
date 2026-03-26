"""Describe binge risk as a function of trajectory phase and set-point distance.

This asks two different questions:
1. Descriptively, how does binge rate vary across falling/stable/rising phases?
2. Predictively, does phase add out-of-sample value beyond distance or yesterday's
   calories in the pre-tirzepatide era?

Outputs:
    analysis/phase_binge_summary.csv
    analysis/phase_binge_distance_bins.csv
    analysis/phase_binge_model_auc.csv
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from C_binge_analysis import roc_auc_score_np, fit_logistic_regression, predict_logistic_regression, standardize_train_test

TREND_WINDOW_DAYS = 90
TREND_THRESHOLD_LBS = 3.0


def load_data():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    daily = intake.merge(
        kalman[["date", "fat_mass_lbs", "fat_mass_lbs_filtered"]],
        on="date",
        how="left",
    )
    daily = daily.merge(tirz[["date", "effective_level"]], on="date", how="left")
    daily["effective_level"] = daily["effective_level"].fillna(0)
    return daily.sort_values("date").reset_index(drop=True)


def build_frame(daily):
    d = daily[daily["effective_level"] == 0].copy()
    d["binge"] = (d["calories"] > 2800).astype(int)
    d["prev_cal"] = d["calories"].shift(1)
    d["fm_pred"] = d["fat_mass_lbs_filtered"].fillna(d["fat_mass_lbs"])
    d["fat_delta_window"] = d["fat_mass_lbs"].diff(TREND_WINDOW_DAYS)
    d["phase"] = np.where(
        d["fat_delta_window"] <= -TREND_THRESHOLD_LBS,
        "falling",
        np.where(d["fat_delta_window"] >= TREND_THRESHOLD_LBS, "rising", "stable"),
    )
    d["phase_falling"] = (d["phase"] == "falling").astype(int)
    d["phase_rising"] = (d["phase"] == "rising").astype(int)
    prev_fm = d["fm_pred"].shift(1)
    d["sp_180d"] = prev_fm.rolling(180, min_periods=90).mean()
    d["dist_180d"] = prev_fm - d["sp_180d"]
    return d


def phase_summary(daily):
    rows = []
    for phase, grp in daily.groupby("phase"):
        rows.append(
            {
                "phase": phase,
                "n_days": len(grp),
                "mean_distance_lbs": round(grp["dist_180d"].mean(), 2),
                "mean_calories": round(grp["calories"].mean(), 1),
                "binge_rate": round(grp["binge"].mean(), 4),
            }
        )
    return pd.DataFrame(rows).sort_values("phase")


def distance_bin_summary(daily):
    rows = []
    for phase in ["falling", "stable", "rising"]:
        grp = daily[daily["phase"] == phase]
        for lo, hi in [(-10, -5), (-5, 0), (0, 5), (5, 10)]:
            sub = grp[(grp["dist_180d"] >= lo) & (grp["dist_180d"] < hi)]
            if len(sub) < 40:
                continue
            rows.append(
                {
                    "phase": phase,
                    "distance_bin": f"{lo} to {hi}",
                    "n_days": len(sub),
                    "mean_calories": round(sub["calories"].mean(), 1),
                    "binge_rate": round(sub["binge"].mean(), 4),
                }
            )
    return pd.DataFrame(rows)


def auc_summary(daily):
    valid = daily.dropna(subset=["dist_180d", "prev_cal"]).copy()
    models = [
        ("distance_only", ["dist_180d"]),
        ("phase_only", ["phase_falling", "phase_rising"]),
        ("distance_plus_phase", ["dist_180d", "phase_falling", "phase_rising"]),
        ("prev_cal_only", ["prev_cal"]),
        ("prev_cal_plus_distance", ["prev_cal", "dist_180d"]),
        ("prev_cal_plus_phase", ["prev_cal", "phase_falling", "phase_rising"]),
        (
            "prev_cal_plus_distance_plus_phase",
            ["prev_cal", "dist_180d", "phase_falling", "phase_rising"],
        ),
    ]

    rows = []
    for label, cols in models:
        result = ba.evaluate_predictor(valid, cols, label)
        if result is None:
            continue
        rows.append(
            {
                "model": label,
                "auc": round(float(result["auc"]), 4),
                "n_days": result["n"],
                "n_binge": result["n_binge"],
            }
        )
    return pd.DataFrame(rows)


def save_outputs(summary_df, bins_df, auc_df):
    summary_df.to_csv(ROOT / "analysis" / "H_phase_binge_summary.csv", index=False)
    bins_df.to_csv(ROOT / "analysis" / "H_phase_binge_distance_bins.csv", index=False)
    auc_df.to_csv(ROOT / "analysis" / "H_phase_binge_model_auc.csv", index=False)


def print_report(summary_df, bins_df, auc_df):
    print("\n=== Phase-Aware Binge Landscape (Pre-tirzepatide) ===")
    for _, row in summary_df.iterrows():
        print(
            f"{row['phase']:>7}: binge rate {row['binge_rate']*100:5.1f}%  "
            f"mean calories {row['mean_calories']:.0f}  mean dist {row['mean_distance_lbs']:+.1f} lbs"
        )

    print("\nDistance bins within phase:")
    for _, row in bins_df.iterrows():
        print(
            f"{row['phase']:>7}  {row['distance_bin']:>9}: binge {row['binge_rate']*100:5.1f}%  "
            f"calories {row['mean_calories']:.0f}  n={int(row['n_days'])}"
        )

    print("\nWalk-forward AUC:")
    for _, row in auc_df.iterrows():
        print(f"  {row['model']:30s}  AUC={row['auc']:.4f}  n={int(row['n_days'])}")


def main():
    daily = build_frame(load_data())
    summary_df = phase_summary(daily.dropna(subset=["dist_180d"]))
    bins_df = distance_bin_summary(daily.dropna(subset=["dist_180d"]))
    auc_df = auc_summary(daily)
    save_outputs(summary_df, bins_df, auc_df)
    print_report(summary_df, bins_df, auc_df)


if __name__ == "__main__":
    main()
