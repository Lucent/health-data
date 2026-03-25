"""Water weight smoothing: glycogen + sodium corrections.

Two independent water-weight mechanisms, both validated by the same criterion
(nutrient→weight correlation drops to ~zero after correction):

1. **Glycogen-water**: glycogen stores seek a target proportional to daily carb
   intake. Each gram of glycogen binds ~3g water. Validated: carb→weight
   partial correlation drops from r=0.27 to r≈0. Variance reduction: 9.6%.

2. **Sodium-water**: dietary sodium causes water retention at ~136ml per gram
   sodium (literature: 130-150ml). Validated: sodium→weight partial correlation
   drops from r=0.175 to r≈-0.015. Variance reduction: 9.3%.

Both corrections are lagged by 1 day (weight measured upon waking BEFORE
eating reflects previous day's intake). Both are centered on median levels
so typical days get ~zero correction.

Inputs:
    intake/intake_daily.csv — daily carbs
    weight/weight.csv — daily weight (sparse)

Outputs:
    analysis/smoothed_weight.csv — observed + smoothed weight with glycogen state
    analysis/glycogen_validation.png — multi-window diagnostic plots
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

# Glycogen parameters (grid search over 840 combinations, 6 validation windows)
G_MAX = 600         # effective glycogen+water capacity (grams)
CARB_REF = 350      # carb intake at which glycogen target = G_MAX
RATE_UP = 0.60      # refill: 60% toward target per day of eating
RATE_DOWN = 0.45    # depletion: 45% toward target per day of restriction
WATER_RATIO = 3.0   # grams water per gram glycogen
GRAMS_PER_LB = 453.592

# Sodium parameters (sweep over 13 values, validated by correlation + variance)
SODIUM_K = 0.00030  # lbs water retained per mg sodium
                    # = 136ml/g Na, literature range 130-150ml/g


def load_data():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    weight = pd.read_csv(ROOT / "weight" / "weight.csv", parse_dates=["date"])
    weight = weight[["date", "weight_lbs"]].dropna(subset=["weight_lbs"])
    return intake, weight


def simulate_glycogen(carbs, g_max=G_MAX, carb_ref=CARB_REF,
                      rate_up=RATE_UP, rate_down=RATE_DOWN):
    """Simulate daily glycogen using linear target-seeking dynamics.

    target(t) = g_max * clamp(carbs(t) / carb_ref, 0, 1)

    With carb_ref=350 and median intake ~233g, most days have target ~400g
    (below max). Only high-carb days (>350g) reach full glycogen. This means
    the model corrects across the entire carb range, not just extreme fasts.
    """
    n = len(carbs)
    glycogen = np.empty(n)
    g = g_max * 0.67  # ~median level
    for i in range(n):
        target = g_max * min(carbs[i] / carb_ref, 1.0)
        if target >= g:
            g += rate_up * (target - g)
        else:
            g += rate_down * (target - g)
        g = np.clip(g, 0, g_max)
        glycogen[i] = g
    return glycogen


def glycogen_to_water_lbs(glycogen_g):
    """Glycogen + bound water mass in lbs."""
    return glycogen_g * (1 + WATER_RATIO) / GRAMS_PER_LB


def build_daily_frame(intake, weight):
    """Day-level DataFrame spanning full intake range, joined with sparse weight."""
    date_range = pd.date_range(intake["date"].min(), intake["date"].max(), freq="D")
    daily = pd.DataFrame({"date": date_range})
    daily = daily.merge(intake[["date", "calories", "carbs_g", "sodium_mg"]], on="date", how="left")
    daily["carbs_g"] = daily["carbs_g"].fillna(0)
    daily["calories"] = daily["calories"].fillna(0)
    daily["sodium_mg"] = daily["sodium_mg"].fillna(daily["sodium_mg"].median())
    daily = daily.merge(weight, on="date", how="left")
    return daily


def consecutive_diffs(dates, values):
    """Day-to-day differences for consecutive dates only."""
    diffs = []
    for i in range(1, len(dates)):
        gap = (dates[i] - dates[i - 1]) / np.timedelta64(1, "D")
        if gap == 1:
            diffs.append(values[i] - values[i - 1])
    return np.array(diffs)


def apply_model(daily):
    """Apply glycogen + sodium corrections with 1-day lag.

    Both corrections are centered on their median levels so typical days
    get ~zero total correction. Both are lagged 1 day (morning weight
    reflects previous day's intake).
    """
    # --- Glycogen correction ---
    glycogen = simulate_glycogen(daily["carbs_g"].values)
    glycogen_lagged = np.roll(glycogen, 1)
    glycogen_lagged[0] = G_MAX * 0.67
    median_glycogen = np.median(glycogen)
    water_median = glycogen_to_water_lbs(median_glycogen)
    water_lagged = glycogen_to_water_lbs(glycogen_lagged)
    glycogen_correction = water_median - water_lagged

    # --- Sodium correction ---
    sodium = daily["sodium_mg"].values
    sodium_lagged = np.roll(sodium, 1)
    sodium_lagged[0] = np.nanmedian(sodium)
    median_sodium = np.nanmedian(sodium)
    sodium_correction = (median_sodium - sodium_lagged) * SODIUM_K

    # --- Combined ---
    total_correction = glycogen_correction + sodium_correction

    daily = daily.copy()
    daily["glycogen_g"] = np.round(glycogen_lagged, 1)
    daily["glycogen_correction_lbs"] = np.round(glycogen_correction, 2)
    daily["sodium_correction_lbs"] = np.round(sodium_correction, 2)
    daily["smoothed_weight_lbs"] = np.round(
        daily["weight_lbs"] + total_correction, 2
    )
    return daily


def eval_window(daily, start, end):
    """Evaluate variance reduction in a date window."""
    mask = ((daily["date"] >= start) & (daily["date"] <= end)
            & daily["weight_lbs"].notna())
    sub = daily[mask]
    if sub.shape[0] < 5:
        return None
    dates = sub["date"].values
    raw_d = consecutive_diffs(dates, sub["weight_lbs"].values)
    smo_d = consecutive_diffs(dates, sub["smoothed_weight_lbs"].values)
    if len(raw_d) < 3:
        return None
    rv, sv = np.var(raw_d), np.var(smo_d)
    return {
        "raw_var": rv, "smooth_var": sv,
        "reduction": (1 - sv / rv) * 100,
        "n": len(raw_d),
    }


def save_output(daily):
    """Save smoothed weight CSV."""
    cols = ["date", "weight_lbs", "glycogen_g", "glycogen_correction_lbs",
            "sodium_correction_lbs", "smoothed_weight_lbs"]
    out = daily[daily["weight_lbs"].notna()][cols].copy()
    out["date"] = out["date"].dt.strftime("%Y-%m-%d")
    out_path = ROOT / "analysis" / "smoothed_weight.csv"
    out.to_csv(out_path, index=False)
    print(f"\nWrote {len(out)} rows to {out_path}")
    return out


def plot_multiwindow(daily):
    """Multi-window diagnostic plot."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
    except ImportError:
        print("matplotlib not available, skipping plot")
        return

    windows = [
        ("Oct-Nov 2019 (weekend fasts)", "2019-10-10", "2019-12-05"),
        ("May-Jul 2022 (67d dense)", "2022-05-01", "2022-07-15"),
        ("Jan-Feb 2023 (45d dense)", "2023-01-01", "2023-02-18"),
        ("Oct-Nov 2021 (24d dense)", "2021-10-27", "2021-11-27"),
        ("Feb-Mar 2024 (36d dense)", "2024-02-05", "2024-03-18"),
        ("Aug-Sep 2025 (32d dense)", "2025-08-10", "2025-09-18"),
    ]

    fig, axes = plt.subplots(len(windows), 1, figsize=(16, 4 * len(windows)))

    for idx, (label, start, end) in enumerate(windows):
        ax = axes[idx]
        mask = (daily["date"] >= start) & (daily["date"] <= end)
        sub = daily[mask]
        w = sub["weight_lbs"].notna()

        ax.plot(sub.loc[w, "date"], sub.loc[w, "weight_lbs"],
                "o-", color="tab:blue", ms=3, lw=1, label="Observed", alpha=0.8)
        ax.plot(sub.loc[w, "date"], sub.loc[w, "smoothed_weight_lbs"],
                "s-", color="tab:red", ms=3, lw=1, label="Smoothed", alpha=0.8)

        # Carbs as background
        ax2 = ax.twinx()
        colors = ["tab:red" if c < 100 else "tab:orange" if c < 200 else "tab:green"
                  for c in sub["carbs_g"]]
        ax2.bar(sub["date"], sub["carbs_g"], color=colors, alpha=0.15, width=0.8)
        ax2.set_ylabel("Carbs (g)", fontsize=8, color="gray")
        ax2.set_ylim(0, 600)
        ax2.tick_params(labelsize=7, colors="gray")

        ev = eval_window(daily, start, end)
        stats = f"var red: {ev['reduction']:.1f}%  n={ev['n']}" if ev else "n/a"
        ax.set_ylabel("Weight (lbs)")
        ax.set_title(f"{label}  |  {stats}", fontsize=10)
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, alpha=0.2)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))

    plt.suptitle(f"Glycogen Smoothing: G_max={G_MAX}  carb_ref={CARB_REF}  "
                 f"rate_up={RATE_UP}  rate_down={RATE_DOWN}",
                 fontsize=12, y=1.01)
    plt.tight_layout()
    out = ROOT / "analysis" / "glycogen_validation.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")
    plt.close()


def main():
    print("Loading data...")
    intake, weight = load_data()
    daily = build_daily_frame(intake, weight)
    print(f"Daily frame: {len(daily)} days, {daily['weight_lbs'].notna().sum()} with weight")

    print(f"\nParameters (tuned via grid search over 6 windows):")
    print(f"  G_max    = {G_MAX}g")
    print(f"  carb_ref = {CARB_REF}g  (glycogen full when carbs >= {CARB_REF}g)")
    print(f"  rate_up  = {RATE_UP}  (refill)")
    print(f"  rate_dn  = {RATE_DOWN}  (depletion)")
    print(f"  Max water swing: {glycogen_to_water_lbs(G_MAX):.1f} lbs")

    daily = apply_model(daily)

    # Report per-window results
    windows = [
        ("Oct-Nov 2019 (fasts)", "2019-10-20", "2019-11-30"),
        ("May-Jul 2022", "2022-05-05", "2022-07-10"),
        ("Jan-Feb 2023", "2023-01-01", "2023-02-14"),
        ("Oct-Nov 2021", "2021-10-31", "2021-11-23"),
        ("Feb-Mar 2024", "2024-02-08", "2024-03-14"),
        ("Aug-Sep 2025", "2025-08-14", "2025-09-14"),
        ("Global", daily["date"].min(), daily["date"].max()),
    ]

    print(f"\n{'Window':<25} {'raw_var':>8} {'smo_var':>8} {'reduction':>9} {'n':>5}")
    print("-" * 60)
    for label, start, end in windows:
        ev = eval_window(daily, start, end)
        if ev:
            print(f"{label:<25} {ev['raw_var']:8.3f} {ev['smooth_var']:8.3f} "
                  f"{ev['reduction']:8.1f}% {ev['n']:5d}")

    # Carb→weight correlation before and after
    mask = daily["weight_lbs"].notna().values
    dates = daily["date"].values[mask]
    raw = daily["weight_lbs"].values[mask]
    smoothed = daily["smoothed_weight_lbs"].values[mask]
    carbs_prev = np.roll(daily["carbs_g"].values, 1)[mask]

    raw_d = consecutive_diffs(dates, raw)
    smo_d = consecutive_diffs(dates, smoothed)
    carb_d = []
    for i in range(1, len(dates)):
        if (dates[i] - dates[i - 1]) / np.timedelta64(1, "D") == 1:
            carb_d.append(carbs_prev[i])
    carb_d = np.array(carb_d)

    r_raw = np.corrcoef(carb_d, raw_d)[0, 1]
    r_smo = np.corrcoef(carb_d, smo_d)[0, 1]
    print(f"\nCarb→weight-change correlation:")
    print(f"  Raw:      r={r_raw:.4f}")
    print(f"  Smoothed: r={r_smo:.4f}")
    print(f"  Reduction: {abs(r_raw) - abs(r_smo):.4f}")

    save_output(daily)
    plot_multiwindow(daily)


if __name__ == "__main__":
    main()
