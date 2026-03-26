"""Test expenditure hysteresis at matched fat mass in the pre-tirzepatide era.

Hypothesis:
At the same fat mass, TDEE differs depending on whether that fat mass was
reached from above (falling phase) or below (rising phase).

This script is retrospective: it uses the RTS-smoothed Kalman states
(`fat_mass_lbs`, `tdee`) rather than the causal filtered columns.

Outputs:
    analysis/tdee_hysteresis_phase_summary.csv
    analysis/tdee_hysteresis_band_summary.csv
    analysis/tdee_hysteresis_regression.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

TREND_WINDOW_DAYS = 90
TREND_THRESHOLD_LBS = 3.0
FAT_BANDS = [(25, 45), (45, 65), (65, 85)]


def load_data():
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])

    daily = kalman.merge(tirz[["date", "effective_level"]], on="date", how="left")
    daily = daily.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    daily["effective_level"] = daily["effective_level"].fillna(0)
    return daily.sort_values("date").reset_index(drop=True)


def classify_phase(delta_lbs, threshold=TREND_THRESHOLD_LBS):
    if pd.isna(delta_lbs):
        return np.nan
    if delta_lbs <= -threshold:
        return "falling"
    if delta_lbs >= threshold:
        return "rising"
    return "stable"


def build_pre_tirz_frame(daily):
    d = daily[daily["effective_level"] == 0].copy()
    d["fat_delta_window"] = d["fat_mass_lbs"].diff(TREND_WINDOW_DAYS)
    d["phase"] = d["fat_delta_window"].map(classify_phase)
    d["tdee_minus_rmr"] = d["tdee"] - d["expected_rmr"]
    d["tdee_rmr_ratio"] = d["tdee"] / d["expected_rmr"]
    return d.dropna(subset=["fat_mass_lbs", "tdee", "expected_rmr", "phase"])


def phase_summary(daily):
    rows = []
    for phase, grp in daily.groupby("phase"):
        rows.append(
            {
                "phase": phase,
                "n_days": len(grp),
                "mean_fat_mass_lbs": round(grp["fat_mass_lbs"].mean(), 2),
                "mean_tdee": round(grp["tdee"].mean(), 1),
                "mean_expected_rmr": round(grp["expected_rmr"].mean(), 1),
                "mean_tdee_minus_rmr": round(grp["tdee_minus_rmr"].mean(), 1),
                "mean_tdee_rmr_ratio": round(grp["tdee_rmr_ratio"].mean(), 4),
            }
        )
    return pd.DataFrame(rows).sort_values("phase")


def matched_band_summary(daily):
    rows = []
    for low, high in FAT_BANDS:
        band = daily[
            (daily["fat_mass_lbs"] >= low)
            & (daily["fat_mass_lbs"] < high)
            & (daily["phase"].isin(["rising", "falling"]))
        ].copy()
        if band.empty:
            continue
        rising = band[band["phase"] == "rising"]
        falling = band[band["phase"] == "falling"]
        if len(rising) < 20 or len(falling) < 20:
            continue
        rows.append(
            {
                "fat_band_lbs": f"{low}-{high}",
                "mean_fat_mass_lbs": round(band["fat_mass_lbs"].mean(), 2),
                "rising_days": len(rising),
                "falling_days": len(falling),
                "rising_mean_tdee": round(rising["tdee"].mean(), 1),
                "falling_mean_tdee": round(falling["tdee"].mean(), 1),
                "falling_minus_rising_tdee": round(falling["tdee"].mean() - rising["tdee"].mean(), 1),
                "rising_mean_ratio": round(rising["tdee_rmr_ratio"].mean(), 4),
                "falling_mean_ratio": round(falling["tdee_rmr_ratio"].mean(), 4),
                "falling_minus_rising_ratio": round(
                    falling["tdee_rmr_ratio"].mean() - rising["tdee_rmr_ratio"].mean(), 4
                ),
                "rising_mean_tdee_minus_rmr": round(rising["tdee_minus_rmr"].mean(), 1),
                "falling_mean_tdee_minus_rmr": round(falling["tdee_minus_rmr"].mean(), 1),
                "falling_minus_rising_tdee_minus_rmr": round(
                    falling["tdee_minus_rmr"].mean() - rising["tdee_minus_rmr"].mean(), 1
                ),
            }
        )
    return pd.DataFrame(rows)


def fit_simple_regression(daily):
    """OLS: TDEE ~ intercept + fat_mass + phase indicators."""
    X = np.column_stack(
        [
            np.ones(len(daily)),
            daily["fat_mass_lbs"].values,
            (daily["phase"] == "falling").astype(float).values,
            (daily["phase"] == "rising").astype(float).values,
        ]
    )
    y = daily["tdee"].values
    coef = np.linalg.lstsq(X, y, rcond=None)[0]
    labels = ["intercept", "fat_mass_lbs", "falling_vs_stable", "rising_vs_stable"]
    out = pd.DataFrame({"term": labels, "coef": np.round(coef, 4)})
    out["interpretation"] = [
        "baseline stable-phase TDEE at zero fat mass (not directly meaningful)",
        "cal/day increase in TDEE per additional lb fat mass, holding phase fixed",
        "phase offset for falling vs stable at matched fat mass",
        "phase offset for rising vs stable at matched fat mass",
    ]
    return out


def print_report(phase_df, band_df, reg_df):
    phase_map = {row["phase"]: row for _, row in phase_df.iterrows()}
    rising = phase_map.get("rising")
    falling = phase_map.get("falling")
    stable = phase_map.get("stable")

    print("\n=== TDEE Hysteresis (Pre-tirzepatide, retrospective Kalman states) ===")
    if rising is not None and falling is not None:
        print(
            f"Overall matched-era comparison: rising mean TDEE {rising['mean_tdee']:.0f} vs "
            f"falling {falling['mean_tdee']:.0f} cal/day"
        )
        print(
            f"TDEE-RMR gap: rising {rising['mean_tdee_minus_rmr']:.0f} vs "
            f"falling {falling['mean_tdee_minus_rmr']:.0f} cal/day"
        )
        print(
            f"TDEE/RMR ratio: rising {rising['mean_tdee_rmr_ratio']:.3f} vs "
            f"falling {falling['mean_tdee_rmr_ratio']:.3f}"
        )

    if not band_df.empty:
        print("\nMatched fat-mass bands:")
        for _, row in band_df.iterrows():
            print(
                f"  FM {row['fat_band_lbs']} lbs: rising {row['rising_mean_tdee']:.0f} vs "
                f"falling {row['falling_mean_tdee']:.0f} cal/day "
                f"(diff {row['falling_minus_rising_tdee']:+.0f})"
            )

    coeff = dict(zip(reg_df["term"], reg_df["coef"]))
    print("\nSimple regression:")
    print(
        f"  TDEE = {coeff['intercept']:.1f} + {coeff['fat_mass_lbs']:.2f}*fat_mass "
        f"+ {coeff['falling_vs_stable']:.1f}*falling + {coeff['rising_vs_stable']:.1f}*rising"
    )


def save_outputs(phase_df, band_df, reg_df):
    phase_df.to_csv(ROOT / "analysis" / "K_tdee_hysteresis_phase_summary.csv", index=False)
    band_df.to_csv(ROOT / "analysis" / "K_tdee_hysteresis_band_summary.csv", index=False)
    reg_df.to_csv(ROOT / "analysis" / "K_tdee_hysteresis_regression.csv", index=False)


def main():
    daily = build_pre_tirz_frame(load_data())
    phase_df = phase_summary(daily)
    band_df = matched_band_summary(daily)
    reg_df = fit_simple_regression(daily)
    save_outputs(phase_df, band_df, reg_df)
    print_report(phase_df, band_df, reg_df)


if __name__ == "__main__":
    main()
