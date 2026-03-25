"""Diagnostic plots for evaluating glycogen smoothing and weight interpolation.

Generates a multi-panel figure showing:
1. Full weight trajectory: observed vs interpolated vs smoothed
2. Zoomed windows showing model fit quality
3. Derived TDEE over time with Mifflin-St Jeor reference
4. Cumulative energy balance residual (model drift detector)
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

sys.path.insert(0, str(ROOT / "analysis"))
from interpolate_weight import CALORIMETRY, CAL_PER_LB

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec


def load():
    daily = pd.read_csv(ROOT / "analysis" / "daily_weight.csv", parse_dates=["date"])
    return daily


def fig1_full_trajectory(daily):
    """Full 15-year weight trajectory: observed, interpolated, smoothed."""
    fig, ax = plt.subplots(figsize=(18, 6))

    obs = daily[daily["observed_weight_lbs"].notna()]
    ax.scatter(obs["date"], obs["observed_weight_lbs"], s=3, color="tab:blue",
               alpha=0.3, zorder=3, label="Observed (scale)")

    ax.plot(daily["date"], daily["interpolated_weight_lbs"],
            color="tab:green", lw=0.4, alpha=0.5, label="Interpolated (simulated scale)")

    # 7-day rolling smoothed for readability
    sm7 = daily["smoothed_weight_lbs"].rolling(7, center=True, min_periods=1).mean()
    ax.plot(daily["date"], sm7, color="tab:red", lw=1.2,
            label="Smoothed fat mass (7-day avg)")

    ax.set_ylabel("Weight (lbs)")
    ax.set_title("Full Weight Trajectory: Observed vs Model")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

    plt.tight_layout()
    out = ROOT / "analysis" / "plot_trajectory.png"
    plt.savefig(out, dpi=150)
    print(f"Saved {out}")
    plt.close()


def fig2_zoom_windows(daily):
    """Zoomed views of 6 interesting periods showing model fit."""
    windows = [
        ("2019-10-10", "2019-12-05", "Oct-Nov 2019: Weekend fasts"),
        ("2022-05-01", "2022-07-15", "May-Jul 2022: Dense weighing, weight loss"),
        ("2017-08-01", "2018-12-31", "2017-2018: 399-day gap (interpolated)"),
        ("2024-09-01", "2025-03-01", "Sep 2024-Feb 2025: Tirzepatide start"),
        ("2013-06-01", "2015-06-01", "2013-2015: 517-day gap (lowest weight)"),
        ("2020-01-01", "2020-12-31", "2020: COVID year"),
    ]

    fig, axes = plt.subplots(3, 2, figsize=(18, 14))
    axes = axes.flatten()

    for idx, (start, end, title) in enumerate(windows):
        ax = axes[idx]
        mask = (daily["date"] >= start) & (daily["date"] <= end)
        sub = daily[mask]
        obs = sub[sub["observed_weight_lbs"].notna()]

        ax.scatter(obs["date"], obs["observed_weight_lbs"], s=12, color="tab:blue",
                   zorder=4, label="Observed", alpha=0.7)
        ax.plot(sub["date"], sub["interpolated_weight_lbs"],
                color="tab:green", lw=0.8, alpha=0.6, label="Interpolated")
        ax.plot(sub["date"], sub["smoothed_weight_lbs"],
                color="tab:red", lw=1.2, label="Smoothed (fat mass)")

        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7, loc="best")
        ax.grid(True, alpha=0.2)
        ax.set_ylabel("lbs", fontsize=8)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b '%y"))

    plt.suptitle("Model Fit: Zoomed Windows", fontsize=13)
    plt.tight_layout()
    out = ROOT / "analysis" / "plot_zoom_windows.png"
    plt.savefig(out, dpi=150)
    print(f"Saved {out}")
    plt.close()


def fig3_tdee(daily):
    """TDEE over time: derived vs Mifflin-St Jeor expected vs intake."""
    fig, axes = plt.subplots(2, 1, figsize=(18, 10), sharex=True,
                             gridspec_kw={"height_ratios": [2, 1]})

    # --- Top: TDEE and intake ---
    ax = axes[0]

    # Load composition-aware expected RMR
    comp_path = ROOT / "analysis" / "daily_composition.csv"
    if comp_path.exists():
        comp = pd.read_csv(comp_path, parse_dates=["date"])
        rmr_ref = daily.merge(comp[["date", "expected_rmr"]], on="date", how="left")["expected_rmr"]
    else:
        # MSJ fallback
        height_cm = 5 * 30.48 + 11.75 * 2.54
        rmr_ref = daily.apply(
            lambda r: 10 * (r["smoothed_weight_lbs"] * 0.4536) + 6.25 * height_cm
            - 5 * ((r["date"] - pd.Timestamp("1982-10-15")).days / 365.25) + 5
            if pd.notna(r["smoothed_weight_lbs"]) else np.nan, axis=1
        )
    rmr_30 = rmr_ref.rolling(30, center=True, min_periods=10).mean()

    tdee_30 = daily["tdee"].rolling(30, center=True, min_periods=10).mean()
    intake_30 = daily["calories"].rolling(30, center=True, min_periods=10).mean()

    ax.plot(daily["date"], intake_30, color="tab:green", lw=1, alpha=0.6,
            label="30-day intake")
    ax.plot(daily["date"], tdee_30, color="tab:purple", lw=1.5,
            label="30-day derived TDEE")
    ax.plot(daily["date"], rmr_30, color="tab:orange", lw=1, ls="--",
            label="30-day expected RMR (composition-aware)")

    # Fill between intake and TDEE: green when surplus, red when deficit
    ax.fill_between(daily["date"], intake_30, tdee_30,
                    where=intake_30 > tdee_30, color="tab:red", alpha=0.1,
                    label="Surplus (gaining)")
    ax.fill_between(daily["date"], intake_30, tdee_30,
                    where=intake_30 <= tdee_30, color="tab:blue", alpha=0.1,
                    label="Deficit (losing)")

    # Calorimetry anchors
    for year, rmr in CALORIMETRY.items():
        ax.plot(pd.Timestamp(f"{year}-06-01"), rmr, "D", color="black",
                ms=8, zorder=5)
        ax.annotate(f"  Measured RMR {year}: {rmr}",
                    xy=(pd.Timestamp(f"{year}-06-01"), rmr),
                    fontsize=7, va="bottom")

    ax.set_ylabel("Calories/day")
    ax.set_ylim(1200, 3500)
    ax.legend(fontsize=8, ncol=3, loc="upper right")
    ax.set_title("Derived TDEE vs Intake vs Expected RMR")
    ax.grid(True, alpha=0.2)

    # --- Bottom: TDEE / MSJ ratio (metabolic efficiency) ---
    ax = axes[1]
    ratio = tdee_30 / rmr_30
    ax.plot(daily["date"], ratio, color="tab:purple", lw=1.2)
    ax.axhline(y=1.0, color="gray", ls="--", alpha=0.5, label="TDEE = RMR")
    ax.axhline(y=1.2, color="tab:orange", ls=":", alpha=0.5,
               label="Typical sedentary (1.2×)")
    ax.fill_between(daily["date"], ratio, 1.2,
                    where=ratio < 1.2, color="tab:blue", alpha=0.15)
    ax.fill_between(daily["date"], ratio, 1.2,
                    where=ratio >= 1.2, color="tab:red", alpha=0.15)

    ax.set_ylabel("TDEE / MSJ ratio")
    ax.set_ylim(0.7, 1.8)
    ax.legend(fontsize=8)
    ax.set_title("Metabolic Efficiency: TDEE relative to composition-aware RMR "
                 "(below 1.2 = adaptation/undercounting)")
    ax.grid(True, alpha=0.2)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

    plt.tight_layout()
    out = ROOT / "analysis" / "plot_tdee.png"
    plt.savefig(out, dpi=150)
    print(f"Saved {out}")
    plt.close()


def fig4_residual(daily):
    """Cumulative energy balance residual: detects model drift.

    If the model is correct, the cumulative (intake - TDEE) should track
    weight change * 3500. Any drift indicates the model or data is off.
    """
    fig, axes = plt.subplots(2, 1, figsize=(18, 8), sharex=True)

    # Compute cumulative predicted vs actual weight change
    valid = daily[daily["tdee"].notna() & daily["smoothed_weight_lbs"].notna()].copy()
    valid["daily_surplus"] = valid["calories"] - valid["tdee"]
    valid["cum_surplus_lbs"] = valid["daily_surplus"].cumsum() / CAL_PER_LB
    valid["actual_change"] = valid["smoothed_weight_lbs"] - valid["smoothed_weight_lbs"].iloc[0]

    ax = axes[0]
    ax.plot(valid["date"], valid["cum_surplus_lbs"], color="tab:green", lw=1,
            label="Cumulative surplus (intake - TDEE) / 3500")
    ax.plot(valid["date"], valid["actual_change"], color="tab:red", lw=1,
            label="Actual smoothed weight change")
    ax.set_ylabel("Pounds")
    ax.legend(fontsize=9)
    ax.set_title("Cumulative Energy Balance: Does the model track reality?")
    ax.grid(True, alpha=0.2)

    # Residual (should be near zero if model is right)
    ax = axes[1]
    residual = valid["cum_surplus_lbs"] - valid["actual_change"]
    ax.plot(valid["date"], residual, color="tab:purple", lw=0.8)
    ax.axhline(y=0, color="gray", ls="--")
    ax.fill_between(valid["date"], 0, residual, alpha=0.2, color="tab:purple")
    ax.set_ylabel("Residual (lbs)")
    ax.set_title("Cumulative Residual (drift = model error or data bias)")
    ax.grid(True, alpha=0.2)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

    plt.tight_layout()
    out = ROOT / "analysis" / "plot_residual.png"
    plt.savefig(out, dpi=150)
    print(f"Saved {out}")
    plt.close()


def main():
    print("Loading daily_weight.csv...")
    daily = load()
    print(f"  {len(daily)} days, {daily['observed_weight_lbs'].notna().sum()} observed\n")

    fig1_full_trajectory(daily)
    fig2_zoom_windows(daily)
    fig3_tdee(daily)
    fig4_residual(daily)

    print("\nAll plots saved to analysis/")


if __name__ == "__main__":
    main()
