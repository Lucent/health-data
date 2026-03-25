"""RMR-based weight interpolation.

Derives TDEE from windows of >=7 days between weigh-ins. Snaps to every
observed weight. Fills gaps by simulating day-by-day energy balance.

For short windows: constant TDEE = (intake - weight_change * 3500) / n_days
For long windows (>60 days): TDEE = k * Mifflin-St Jeor(weight, age), with k
fitted via binary search to match the endpoint.

Inputs:
    intake/intake_daily.csv
    weight/weight.csv

Outputs:
    analysis/daily_weight.csv — complete 5,429-day series
    analysis/interpolation_diagnostic.png
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analysis"))
from glycogen_smooth import (
    simulate_glycogen, glycogen_to_water_lbs,
    G_MAX, CARB_REF, RATE_UP, RATE_DOWN, WATER_RATIO, GRAMS_PER_LB,
)

CAL_PER_LB = 3500
LONG_GAP_THRESHOLD = 60
MIN_TDEE_WINDOW = 7  # minimum days for TDEE derivation

CALORIMETRY = {"2011": 2415, "2012": 1956, "2016": 1700}

# Composition-aware RMR (loaded at runtime from daily_composition.csv)
_COMP_RMR = None  # dict: date_idx -> expected_rmr


def load_data():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    weight = pd.read_csv(ROOT / "weight" / "weight.csv", parse_dates=["date"])
    weight = weight[["date", "weight_lbs"]].dropna(subset=["weight_lbs"])
    return intake, weight


def load_composition_rmr():
    """Load composition-aware expected RMR from rmr_model.py output."""
    global _COMP_RMR
    path = ROOT / "analysis" / "daily_composition.csv"
    if not path.exists():
        print(f"  WARNING: {path} not found, falling back to MSJ")
        return False
    comp = pd.read_csv(path, parse_dates=["date"])
    _COMP_RMR = dict(zip(comp["date"], comp["expected_rmr"]))
    print(f"  Loaded {len(_COMP_RMR)} days of composition-aware RMR")
    return True


def expected_rmr(weight_lbs, date):
    """Get expected RMR: composition-aware if available, else MSJ fallback."""
    if _COMP_RMR is not None:
        rmr = _COMP_RMR.get(date)
        if rmr is not None and not np.isnan(rmr):
            return rmr
    # MSJ fallback
    height_cm = 5 * 30.48 + 11.75 * 2.54
    kg = weight_lbs * 0.453592
    age = (date - pd.Timestamp("1982-10-15")).days / 365.25
    return 10 * kg + 6.25 * height_cm - 5 * age + 5


def build_daily(intake, weight):
    date_range = pd.date_range(intake["date"].min(), intake["date"].max(), freq="D")
    daily = pd.DataFrame({"date": date_range})
    daily = daily.merge(intake[["date", "calories", "carbs_g", "sodium_mg"]], on="date", how="left")
    daily["carbs_g"] = daily["carbs_g"].fillna(0)
    daily["calories"] = daily["calories"].fillna(0)
    daily["sodium_mg"] = daily["sodium_mg"].fillna(daily["sodium_mg"].median())
    daily = daily.merge(weight, on="date", how="left")
    return daily


def compute_water_corrections(daily):
    """Compute glycogen + sodium water-weight corrections.

    Matches glycogen_smooth.py's apply_model() exactly.
    """
    from glycogen_smooth import SODIUM_K

    # Glycogen
    glycogen = simulate_glycogen(daily["carbs_g"].values)
    glycogen_lagged = np.roll(glycogen, 1)
    glycogen_lagged[0] = G_MAX * 0.67
    median_glycogen = np.median(glycogen)
    water_median = glycogen_to_water_lbs(median_glycogen)
    water_lagged = glycogen_to_water_lbs(glycogen_lagged)
    glycogen_correction = water_median - water_lagged

    # Sodium
    sodium = daily["sodium_mg"].values
    sodium_lagged = np.roll(sodium, 1)
    sodium_lagged[0] = np.nanmedian(sodium)
    median_sodium = np.nanmedian(sodium)
    sodium_correction = (median_sodium - sodium_lagged) * SODIUM_K

    daily = daily.copy()
    daily["glycogen_g"] = np.round(glycogen_lagged, 1)
    daily["glycogen_correction_lbs"] = np.round(glycogen_correction + sodium_correction, 4)
    return daily


def find_tdee_windows(daily):
    """Find windows for TDEE derivation.

    Groups consecutive weigh-ins into windows of >= MIN_TDEE_WINDOW days.
    Each window yields one TDEE estimate. Returns list of
    (start_idx, end_idx, n_days) where start and end are observed weight indices.
    """
    obs_idx = daily.index[daily["weight_lbs"].notna()].tolist()
    if len(obs_idx) < 2:
        return []

    windows = []
    i = 0
    while i < len(obs_idx) - 1:
        s = obs_idx[i]
        # Extend window until we have at least MIN_TDEE_WINDOW days
        j = i + 1
        while j < len(obs_idx):
            e = obs_idx[j]
            n_days = (daily.loc[e, "date"] - daily.loc[s, "date"]).days
            if n_days >= MIN_TDEE_WINDOW:
                windows.append({"start_idx": s, "end_idx": e, "n_days": n_days,
                                "start_date": daily.loc[s, "date"],
                                "end_date": daily.loc[e, "date"],
                                "start_weight": daily.loc[s, "weight_lbs"],
                                "end_weight": daily.loc[e, "weight_lbs"]})
                i = j
                break
            j += 1
        else:
            # Ran out of observations without reaching MIN_TDEE_WINDOW
            # Use whatever we have
            e = obs_idx[-1]
            n_days = (daily.loc[e, "date"] - daily.loc[s, "date"]).days
            if n_days > 0:
                windows.append({"start_idx": s, "end_idx": e, "n_days": n_days,
                                "start_date": daily.loc[s, "date"],
                                "end_date": daily.loc[e, "date"],
                                "start_weight": daily.loc[s, "weight_lbs"],
                                "end_weight": daily.loc[e, "weight_lbs"]})
            break

    return windows


def derive_tdee_const(daily, w):
    """Derive constant TDEE for a window.

    Weight on day e reflects food eaten on days s through e-1 (you weigh
    before eating on day e). So total intake is days s..e-1, not s..e.
    """
    s, e = w["start_idx"], w["end_idx"]
    sw = daily.loc[s, "weight_lbs"] + daily.loc[s, "glycogen_correction_lbs"]
    ew = daily.loc[e, "weight_lbs"] + daily.loc[e, "glycogen_correction_lbs"]
    weight_change = ew - sw
    total_intake = daily.loc[s:e - 1, "calories"].sum()
    return (total_intake - weight_change * CAL_PER_LB) / w["n_days"]


def derive_tdee_long(daily, w):
    """Derive k factor for long windows: TDEE = k * expected_RMR."""
    s, e = w["start_idx"], w["end_idx"]
    start_smooth = daily.loc[s, "weight_lbs"] + daily.loc[s, "glycogen_correction_lbs"]
    target_smooth = daily.loc[e, "weight_lbs"] + daily.loc[e, "glycogen_correction_lbs"]

    def simulate_k(k):
        wt = start_smooth
        for idx in range(s, e):
            tdee = k * expected_rmr(wt, daily.loc[idx, "date"])
            wt += (daily.loc[idx, "calories"] - tdee) / CAL_PER_LB
        return wt

    k_lo, k_hi = 0.5, 2.0
    for _ in range(60):
        k_mid = (k_lo + k_hi) / 2
        if simulate_k(k_mid) > target_smooth:
            k_lo = k_mid
        else:
            k_hi = k_mid
    return (k_lo + k_hi) / 2


def interpolate_all(daily, windows):
    """Fill every day with interpolated weight and TDEE.

    Strategy:
    - Derive TDEE per window (constant for short, MSJ-scaled for long)
    - Between each pair of observed weigh-ins within a window, simulate
      day-by-day using the window's TDEE, snapping at each observation.
    - The simulation resets at every observed weight.
    """
    n = len(daily)
    smoothed = np.full(n, np.nan)
    interpolated = np.full(n, np.nan)
    tdee_arr = np.full(n, np.nan)
    window_id = np.full(n, -1, dtype=int)

    # Pre-compute all obs indices for quick lookup
    obs_mask = daily["weight_lbs"].notna().values

    for wid, w in enumerate(windows):
        s, e = w["start_idx"], w["end_idx"]
        is_long = w["n_days"] > LONG_GAP_THRESHOLD

        if is_long:
            k = derive_tdee_long(daily, w)
            w["k_factor"] = k
        else:
            const_tdee = derive_tdee_const(daily, w)
            w["const_tdee"] = const_tdee

        # Walk through each day in the window, snapping at observations
        current_smooth = (daily.loc[s, "weight_lbs"]
                          + daily.loc[s, "glycogen_correction_lbs"])

        for idx in range(s, e + 1):
            # If this day has an observation, snap to it
            if obs_mask[idx]:
                current_smooth = (daily.loc[idx, "weight_lbs"]
                                  + daily.loc[idx, "glycogen_correction_lbs"])

            # Compute TDEE for this day
            if is_long:
                day_tdee = k * expected_rmr(current_smooth, daily.loc[idx, "date"])
            else:
                day_tdee = const_tdee

            smoothed[idx] = current_smooth
            interpolated[idx] = current_smooth - daily.loc[idx, "glycogen_correction_lbs"]
            tdee_arr[idx] = day_tdee
            window_id[idx] = wid

            # Advance underlying weight (for non-observed days ahead)
            if idx < e:
                surplus = daily.loc[idx, "calories"] - day_tdee
                current_smooth += surplus / CAL_PER_LB

    daily = daily.copy()
    daily["smoothed_weight_lbs"] = np.round(smoothed, 2)
    daily["interpolated_weight_lbs"] = np.round(interpolated, 2)
    daily["tdee"] = np.round(tdee_arr, 0)
    daily["window_id"] = window_id
    return daily, windows


def validate(daily, windows):
    """Validation report."""
    print("\n=== Validation ===\n")

    # Interpolation error at observed points
    obs = daily[daily["weight_lbs"].notna() & daily["interpolated_weight_lbs"].notna()]
    errors = (obs["interpolated_weight_lbs"] - obs["weight_lbs"]).abs()
    print(f"Interpolation error at observed points:")
    print(f"  Max: {errors.max():.2f} lbs  Mean: {errors.mean():.3f}  "
          f"Median: {errors.median():.3f}  >1 lb: {(errors > 1).sum()}")

    # TDEE distribution
    valid_tdee = daily["tdee"].dropna()
    print(f"\nTDEE distribution ({len(valid_tdee)} days):")
    print(f"  Mean: {valid_tdee.mean():.0f}  Median: {valid_tdee.median():.0f}  "
          f"Std: {valid_tdee.std():.0f}")
    print(f"  P5: {valid_tdee.quantile(0.05):.0f}  P95: {valid_tdee.quantile(0.95):.0f}")
    outliers = valid_tdee[(valid_tdee < 1400) | (valid_tdee > 3000)]
    print(f"  Outside 1400-3000: {len(outliers)} ({len(outliers)/len(valid_tdee)*100:.1f}%)")

    # Calorimetry
    print(f"\nCalorimetry anchors:")
    for year, rmr in CALORIMETRY.items():
        mask = daily["date"].dt.year == int(year)
        year_tdee = daily.loc[mask, "tdee"].dropna()
        if len(year_tdee) > 0:
            med = year_tdee.median()
            print(f"  {year}: RMR={rmr}  TDEE={med:.0f}  ratio={med/rmr:.2f}")

    # Window stats
    short = [w for w in windows if w["n_days"] <= LONG_GAP_THRESHOLD]
    long = [w for w in windows if w["n_days"] > LONG_GAP_THRESHOLD]
    print(f"\nWindows: {len(windows)} ({len(short)} short, {len(long)} long)")

    # Coverage
    filled = daily["smoothed_weight_lbs"].notna().sum()
    print(f"Days with interpolated weight: {filled}/{len(daily)} "
          f"({filled/len(daily)*100:.1f}%)")


def save_output(daily):
    cols = ["date", "calories", "carbs_g", "observed_weight_lbs",
            "glycogen_g", "glycogen_correction_lbs",
            "interpolated_weight_lbs", "smoothed_weight_lbs", "tdee",
            "window_id"]
    out = daily.rename(columns={"weight_lbs": "observed_weight_lbs"}).copy()
    out["observed_weight_lbs"] = out["observed_weight_lbs"].round(1)
    out["date"] = out["date"].dt.strftime("%Y-%m-%d")
    out_path = ROOT / "analysis" / "daily_weight.csv"
    out[cols].to_csv(out_path, index=False)
    print(f"\nWrote {len(out)} rows to {out_path}")


def plot_diagnostic(daily):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
    except ImportError:
        print("matplotlib not available")
        return

    fig, axes = plt.subplots(2, 1, figsize=(18, 10), sharex=True)

    # Weight
    ax = axes[0]
    obs = daily[daily["weight_lbs"].notna()]
    ax.scatter(obs["date"], obs["weight_lbs"], s=2, color="tab:blue",
               alpha=0.4, label="Observed", zorder=3)
    ax.plot(daily["date"], daily["smoothed_weight_lbs"],
            color="tab:red", lw=0.5, alpha=0.7, label="Smoothed (fat mass)")
    ax.set_ylabel("Weight (lbs)")
    ax.legend(fontsize=8)
    ax.set_title("Complete Daily Weight Series — Observed + Interpolated")
    ax.grid(True, alpha=0.2)

    # TDEE
    ax = axes[1]
    tdee = daily[["date", "tdee", "calories"]].copy()
    tdee["tdee_30d"] = tdee["tdee"].rolling(30, center=True, min_periods=10).mean()
    intake_30d = tdee["calories"].rolling(30, center=True, min_periods=10).mean()

    ax.plot(tdee["date"], tdee["tdee"], color="gray", lw=0.3, alpha=0.2)
    ax.plot(tdee["date"], tdee["tdee_30d"], color="tab:purple", lw=1.5,
            label="30-day TDEE")
    ax.plot(daily["date"], intake_30d, color="tab:green", lw=0.8, alpha=0.5,
            label="30-day intake")

    for year, rmr in CALORIMETRY.items():
        ax.plot(pd.Timestamp(f"{year}-01-01"), rmr, "D", color="tab:orange",
                ms=8, zorder=5)
        ax.annotate(f"RMR {year}: {rmr}", xy=(pd.Timestamp(f"{year}-03-01"), rmr),
                    fontsize=7, color="tab:orange")

    ax.set_ylabel("Calories/day")
    ax.set_ylim(1000, 3500)
    ax.legend(fontsize=8)
    ax.set_title("Derived TDEE vs Intake (30-day rolling)")
    ax.grid(True, alpha=0.2)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

    plt.tight_layout()
    out = ROOT / "analysis" / "interpolation_diagnostic.png"
    plt.savefig(out, dpi=150)
    print(f"Saved {out}")
    plt.close()


def main():
    print("Loading data...")
    intake, weight = load_data()
    weight = weight[weight["date"] >= intake["date"].min()]
    daily = build_daily(intake, weight)
    print(f"Daily frame: {len(daily)} days, {daily['weight_lbs'].notna().sum()} with weight")

    print("Computing water corrections (glycogen + sodium)...")
    daily = compute_water_corrections(daily)

    print("Loading composition-aware RMR...")
    load_composition_rmr()

    print(f"Finding TDEE windows (min {MIN_TDEE_WINDOW} days)...")
    windows = find_tdee_windows(daily)
    short = sum(1 for w in windows if w["n_days"] <= LONG_GAP_THRESHOLD)
    print(f"  {len(windows)} windows ({short} short, {len(windows)-short} long)")

    print("Interpolating...")
    daily, windows = interpolate_all(daily, windows)

    validate(daily, windows)
    save_output(daily)
    plot_diagnostic(daily)


if __name__ == "__main__":
    main()
