"""Derive activity calories from daily steps.

Converts step counts to estimated activity energy expenditure (AEE) using:
    AEE = steps × stride_length_m × weight_kg × efficiency

Then decomposes derived TDEE into:
    TDEE = RMR + TEF + AEE + NEAT + residual

Where:
- RMR: from composition-aware model (daily_composition.csv)
- TEF: thermic effect of food, ~10% of intake
- AEE: from steps (this script)
- NEAT: non-exercise activity thermogenesis (fidgeting, posture, etc.)
- Residual: unexplained (metabolic adaptation, undercounting, or model error)

Inputs:
    steps-sleep/steps.csv — 4,275 days of step counts
    analysis/daily_weight.csv — daily weight + derived TDEE
    analysis/daily_composition.csv — daily expected RMR

Outputs:
    analysis/daily_activity.csv — steps, activity_cal, tdee_decomposition
    Prints calibration report
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

# Walking energy cost parameters
# Gross cost of walking ≈ 0.5 kcal/kg/km (net, above resting)
# Average stride length for 5'11.75" male ≈ 0.78m
# But we calibrate against the data rather than assuming

STRIDE_M_DEFAULT = 0.78  # meters per step (5'11.75" male)
WALK_COST_DEFAULT = 0.5  # kcal per kg per km (net above resting)


def load_data():
    steps = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])
    daily = pd.read_csv(ROOT / "analysis" / "daily_weight.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "daily_composition.csv", parse_dates=["date"])
    return steps, daily, comp


def compute_activity_cal(steps_arr, weight_lbs_arr, distance_arr,
                         cal_per_step_per_lb):
    """Compute activity calories from steps.

    Uses a simple linear model: activity_cal = steps * weight_lbs * k
    where k is calibrated against the data.
    """
    return steps_arr * weight_lbs_arr * cal_per_step_per_lb


def calibrate(merged):
    """Find cal_per_step_per_lb that best explains the TDEE-RMR-TEF gap.

    TDEE = RMR + TEF + AEE + NEAT
    AEE = steps * weight * k
    NEAT ≈ constant

    We fit k and NEAT simultaneously:
    TDEE - RMR - TEF = k * steps * weight + NEAT

    Using least squares on days with both steps and reliable TDEE.
    """
    valid = merged.dropna(subset=["tdee", "expected_rmr", "steps",
                                   "smoothed_weight_lbs"]).copy()
    # Exclude extreme TDEE windows (< 1200 or > 3200)
    valid = valid[(valid["tdee"] > 1200) & (valid["tdee"] < 3200)]

    # TEF ≈ 10% of intake
    valid["tef"] = valid["calories"] * 0.10
    valid["gap"] = valid["tdee"] - valid["expected_rmr"] - valid["tef"]

    # gap = k * steps * weight + NEAT
    X = np.column_stack([
        valid["steps"].values * valid["smoothed_weight_lbs"].values,
        np.ones(len(valid)),
    ])
    y = valid["gap"].values

    coeffs, residuals, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    k, neat = coeffs

    predicted = X @ coeffs
    residual = y - predicted
    rmse = np.sqrt(np.mean(residual ** 2))

    print(f"\n=== Activity Calibration ===")
    print(f"  Fitted on {len(valid)} days (TDEE between 1200-3200)")
    print(f"  cal_per_step_per_lb = {k:.6f}")
    print(f"  NEAT (intercept)    = {neat:.0f} cal/day")
    print(f"  RMSE                = {rmse:.0f} cal/day")
    print(f"")

    # Sanity check: what does this give for typical values?
    for step_count in [2000, 4000, 6000, 8000, 10000]:
        aee = k * step_count * 210
        print(f"  {step_count:5d} steps at 210 lbs → {aee:.0f} cal AEE")

    # Check: physical units
    # k * steps * weight_lbs = kcal
    # k has units of kcal / (step * lb)
    # Convert to kcal/kg/km for comparison:
    # k * steps * weight_lbs = k * (distance_km / stride_km) * (weight_kg / 0.4536)
    stride_km = STRIDE_M_DEFAULT / 1000
    kcal_per_kg_per_km = k / stride_km * 0.4536
    print(f"\n  Equivalent: {kcal_per_kg_per_km:.2f} kcal/kg/km "
          f"(literature: 0.4-0.6 for walking)")

    return k, neat


def build_daily_activity(steps, daily, comp, k, neat):
    """Build complete daily activity table."""
    # Start with daily_weight as base
    result = daily[["date", "calories", "smoothed_weight_lbs", "tdee"]].copy()
    result = result.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    result = result.merge(steps[["date", "steps", "distance"]], on="date", how="left")

    # Activity calories
    has_steps = result["steps"].notna() & result["smoothed_weight_lbs"].notna()
    result["activity_cal"] = np.where(
        has_steps,
        result["steps"] * result["smoothed_weight_lbs"] * k,
        np.nan
    )

    # TEF
    result["tef_cal"] = result["calories"] * 0.10

    # NEAT (constant)
    result["neat_cal"] = neat

    # Residual: TDEE - RMR - TEF - AEE - NEAT
    result["residual_cal"] = (result["tdee"]
                              - result["expected_rmr"]
                              - result["tef_cal"]
                              - result["activity_cal"]
                              - result["neat_cal"])

    # Round
    for col in ["activity_cal", "tef_cal", "neat_cal", "residual_cal"]:
        result[col] = result[col].round(0)

    return result


def report(result):
    """Print decomposition report."""
    valid = result.dropna(subset=["tdee", "expected_rmr", "activity_cal"]).copy()
    valid = valid[(valid["tdee"] > 1200) & (valid["tdee"] < 3200)]

    print(f"\n=== TDEE Decomposition ({len(valid)} days) ===")
    print(f"{'Component':<15} {'Mean':>7} {'Median':>7} {'P5':>7} {'P95':>7}")
    for col, label in [("expected_rmr", "RMR"),
                       ("tef_cal", "TEF"),
                       ("activity_cal", "Activity"),
                       ("neat_cal", "NEAT"),
                       ("residual_cal", "Residual"),
                       ("tdee", "TDEE (total)")]:
        vals = valid[col].dropna()
        print(f"{label:<15} {vals.mean():7.0f} {vals.median():7.0f} "
              f"{vals.quantile(0.05):7.0f} {vals.quantile(0.95):7.0f}")

    # Does the residual correlate with anything?
    print(f"\n=== Residual correlations ===")
    for col, label in [("steps", "Steps"), ("calories", "Intake"),
                       ("smoothed_weight_lbs", "Weight")]:
        r = valid[["residual_cal", col]].dropna()
        if len(r) > 10:
            corr = r.corr().iloc[0, 1]
            print(f"  Residual vs {label}: r={corr:.3f}")

    # By year
    print(f"\n=== Activity calories by year ===")
    print(f"{'Year':>5} {'Steps':>6} {'AEE':>5} {'TDEE':>5} {'RMR':>5} {'Resid':>6} {'n':>5}")
    for year in range(2014, 2027):
        mask = valid["date"].dt.year == year
        if mask.sum() > 10:
            yr = valid[mask]
            print(f"{year:5d} {yr['steps'].median():6.0f} {yr['activity_cal'].median():5.0f} "
                  f"{yr['tdee'].median():5.0f} {yr['expected_rmr'].median():5.0f} "
                  f"{yr['residual_cal'].median():6.0f} {mask.sum():5d}")


def save_output(result):
    cols = ["date", "steps", "activity_cal", "tef_cal", "neat_cal",
            "expected_rmr", "tdee", "residual_cal"]
    out = result[cols].copy()
    out["date"] = out["date"].dt.strftime("%Y-%m-%d")
    path = ROOT / "analysis" / "daily_activity.csv"
    out.to_csv(path, index=False)
    step_days = out["steps"].notna().sum()
    print(f"\nWrote {len(out)} rows ({step_days} with steps) to {path}")


def plot_activity(result):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
    except ImportError:
        print("matplotlib not available")
        return

    valid = result.dropna(subset=["activity_cal", "tdee", "expected_rmr"]).copy()
    valid = valid[(valid["tdee"] > 1200) & (valid["tdee"] < 3200)]

    fig, axes = plt.subplots(3, 1, figsize=(18, 12), sharex=True)

    # Panel 1: TDEE decomposition (stacked)
    ax = axes[0]
    r30 = lambda s: s.rolling(30, center=True, min_periods=10).mean()
    ax.fill_between(valid["date"], 0, r30(valid["expected_rmr"]),
                    alpha=0.3, color="tab:blue", label="RMR")
    ax.fill_between(valid["date"], r30(valid["expected_rmr"]),
                    r30(valid["expected_rmr"] + valid["tef_cal"]),
                    alpha=0.3, color="tab:green", label="TEF")
    ax.fill_between(valid["date"],
                    r30(valid["expected_rmr"] + valid["tef_cal"]),
                    r30(valid["expected_rmr"] + valid["tef_cal"] + valid["activity_cal"]),
                    alpha=0.3, color="tab:orange", label="Activity (steps)")
    ax.plot(valid["date"], r30(valid["tdee"]),
            color="tab:purple", lw=1.5, label="Derived TDEE (30d)")
    ax.set_ylabel("Calories/day")
    ax.set_ylim(1000, 3500)
    ax.legend(fontsize=8, ncol=4)
    ax.set_title("TDEE Decomposition: RMR + TEF + Activity from Steps")
    ax.grid(True, alpha=0.2)

    # Panel 2: Steps
    ax = axes[1]
    ax.bar(valid["date"], valid["steps"], width=1, color="tab:orange", alpha=0.3)
    ax.plot(valid["date"], r30(valid["steps"]),
            color="tab:orange", lw=1.5, label="30-day avg steps")
    ax.set_ylabel("Steps/day")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2)

    # Panel 3: Residual
    ax = axes[2]
    resid_30 = r30(valid["residual_cal"])
    ax.plot(valid["date"], resid_30, color="tab:purple", lw=1.2)
    ax.axhline(y=0, color="gray", ls="--")
    ax.fill_between(valid["date"], 0, resid_30, alpha=0.2, color="tab:purple")
    ax.set_ylabel("Residual (cal/day)")
    ax.set_title("Residual: TDEE - RMR - TEF - Activity - NEAT "
                 "(positive = unexplained burn, negative = undercounting)")
    ax.grid(True, alpha=0.2)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

    plt.tight_layout()
    out = ROOT / "analysis" / "plot_activity.png"
    plt.savefig(out, dpi=150)
    print(f"Saved {out}")
    plt.close()


def main():
    print("Loading data...")
    steps, daily, comp = load_data()
    print(f"  Steps: {len(steps)} days ({steps['date'].min().strftime('%Y-%m-%d')} "
          f"to {steps['date'].max().strftime('%Y-%m-%d')})")

    # Merge for calibration
    merged = daily.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    merged = merged.merge(steps[["date", "steps", "distance"]], on="date", how="left")

    print("Calibrating activity model...")
    k, neat = calibrate(merged)

    print("Building daily activity table...")
    result = build_daily_activity(steps, daily, comp, k, neat)

    report(result)
    save_output(result)
    plot_activity(result)


if __name__ == "__main__":
    main()
