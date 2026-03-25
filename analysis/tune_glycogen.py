"""Parameter tuning for glycogen smoothing model.

Fine grid search around best region, multi-window diagnostic plots.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
WATER_RATIO = 3.0
GRAMS_PER_LB = 453.592


def load_data():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    weight = pd.read_csv(ROOT / "weight" / "weight.csv", parse_dates=["date"])
    weight = weight[["date", "weight_lbs"]].dropna(subset=["weight_lbs"])
    return intake, weight


def build_daily(intake, weight):
    date_range = pd.date_range(intake["date"].min(), intake["date"].max(), freq="D")
    daily = pd.DataFrame({"date": date_range})
    daily = daily.merge(intake[["date", "calories", "carbs_g"]], on="date", how="left")
    daily["carbs_g"] = daily["carbs_g"].fillna(0)
    daily["calories"] = daily["calories"].fillna(0)
    daily = daily.merge(weight, on="date", how="left")
    return daily


def simulate(carbs, g_max, carb_ref, rate_up, rate_down):
    n = len(carbs)
    glycogen = np.empty(n)
    g = g_max * 0.8
    for i in range(n):
        target = g_max * min(carbs[i] / carb_ref, 1.0)
        if target >= g:
            g += rate_up * (target - g)
        else:
            g += rate_down * (target - g)
        g = np.clip(g, 0, g_max)
        glycogen[i] = g
    return glycogen


def apply_model(daily, g_max, carb_ref, rate_up, rate_down):
    glycogen = simulate(daily["carbs_g"].values, g_max, carb_ref, rate_up, rate_down)
    glycogen_lagged = np.roll(glycogen, 1)
    glycogen_lagged[0] = g_max * 0.8
    water_full = g_max * (1 + WATER_RATIO) / GRAMS_PER_LB
    water_lagged = glycogen_lagged * (1 + WATER_RATIO) / GRAMS_PER_LB
    correction = water_full - water_lagged
    daily = daily.copy()
    daily["glycogen_g"] = glycogen_lagged
    daily["correction"] = correction
    daily["smoothed"] = daily["weight_lbs"] + correction
    return daily


def consecutive_diffs(dates, values):
    raw_d = []
    for i in range(1, len(dates)):
        gap = (dates[i] - dates[i - 1]) / np.timedelta64(1, "D")
        if gap == 1:
            raw_d.append(values[i] - values[i - 1])
    return np.array(raw_d)


def eval_window(daily, start, end):
    mask = ((daily["date"] >= start) & (daily["date"] <= end)
            & daily["weight_lbs"].notna())
    sub = daily[mask]
    if sub.shape[0] < 5:
        return None
    dates = sub["date"].values
    raw_d = consecutive_diffs(dates, sub["weight_lbs"].values)
    smo_d = consecutive_diffs(dates, sub["smoothed"].values)
    if len(raw_d) < 3:
        return None
    return {
        "raw_var": np.var(raw_d),
        "smooth_var": np.var(smo_d),
        "reduction": (1 - np.var(smo_d) / np.var(raw_d)) * 100,
        "n": len(raw_d),
        "raw_std": np.std(raw_d),
        "smooth_std": np.std(smo_d),
    }


def main():
    intake, weight = load_data()
    daily = build_daily(intake, weight)

    windows = [
        ("Oct-Nov 2019 fasts", "2019-10-20", "2019-11-30"),
        ("May-Jul 2022", "2022-05-05", "2022-07-10"),
        ("Jan-Feb 2023", "2023-01-01", "2023-02-14"),
        ("Oct-Nov 2021", "2021-10-31", "2021-11-23"),
        ("Feb-Mar 2024", "2024-02-08", "2024-03-14"),
        ("Aug-Sep 2025", "2025-08-14", "2025-09-14"),
        ("Jan-Feb 2025 (tirz)", "2025-01-02", "2025-02-09"),
        ("Oct-Nov 2024 (tirz)", "2024-10-23", "2024-11-11"),
    ]

    # Fine grid around best region
    g_maxes = [400, 450, 500, 550, 600, 650]
    carb_refs = [250, 275, 300, 325, 350]
    rate_ups = [0.40, 0.50, 0.60, 0.70]
    rate_downs = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]

    print(f"Sweeping {len(g_maxes)*len(carb_refs)*len(rate_ups)*len(rate_downs)} combos\n")

    results = []
    for gm in g_maxes:
        for cr in carb_refs:
            for ru in rate_ups:
                for rd in rate_downs:
                    d = apply_model(daily, gm, cr, ru, rd)
                    row = {"g_max": gm, "carb_ref": cr, "rate_up": ru, "rate_down": rd}
                    # Global
                    g_eval = eval_window(d, daily["date"].min(), daily["date"].max())
                    row["global_red"] = g_eval["reduction"] if g_eval else -999
                    # Per window
                    for label, s, e in windows:
                        w_eval = eval_window(d, s, e)
                        row[label] = w_eval["reduction"] if w_eval else None
                    results.append(row)

    df = pd.DataFrame(results)

    # Score: average reduction across non-tirz windows, penalize negative
    non_tirz = [w[0] for w in windows if "tirz" not in w[0]]
    df["avg_non_tirz"] = df[non_tirz].mean(axis=1)
    df["min_non_tirz"] = df[non_tirz].min(axis=1)
    # Penalize if worst non-tirz window is negative
    df["score"] = df["avg_non_tirz"] - df["min_non_tirz"].clip(upper=0).abs()

    df_sorted = df.sort_values("score", ascending=False)

    print("=== Top 20 by score (avg non-tirz reduction, penalized if any window negative) ===\n")
    cols = ["g_max", "carb_ref", "rate_up", "rate_down", "global_red", "score", "min_non_tirz"]
    header = f"{'gmax':>4} {'cr':>3} {'ru':>4} {'rd':>4} {'glob%':>6} {'score':>6} {'worst':>6}"
    for w in windows:
        header += f" {w[0][:12]:>12}"
    print(header)

    for _, r in df_sorted.head(20).iterrows():
        line = (f"{r['g_max']:4.0f} {r['carb_ref']:3.0f} {r['rate_up']:4.2f} {r['rate_down']:4.2f} "
                f"{r['global_red']:6.1f} {r['score']:6.1f} {r['min_non_tirz']:6.1f}")
        for w in windows:
            v = r[w[0]]
            line += f" {v:12.1f}" if v is not None else "          n/a"
        print(line)

    # Best params
    best = df_sorted.iloc[0]
    print(f"\n=== BEST: g_max={best['g_max']:.0f} carb_ref={best['carb_ref']:.0f} "
          f"rate_up={best['rate_up']:.2f} rate_down={best['rate_down']:.2f} ===")
    print(f"Global variance reduction: {best['global_red']:.1f}%")

    # Generate diagnostic plots
    gm = int(best["g_max"])
    cr = int(best["carb_ref"])
    ru = best["rate_up"]
    rd = best["rate_down"]

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
    except ImportError:
        print("no matplotlib")
        return

    d = apply_model(daily, gm, cr, ru, rd)

    plot_windows = [
        ("Oct-Nov 2019 (fasts)", "2019-10-10", "2019-12-05"),
        ("May-Jul 2022 (dense, no fasts)", "2022-05-01", "2022-07-15"),
        ("Jan-Feb 2023 (dense)", "2023-01-01", "2023-02-18"),
        ("Oct-Nov 2021 (dense)", "2021-10-27", "2021-11-27"),
        ("Feb-Mar 2024 (dense)", "2024-02-05", "2024-03-18"),
        ("Aug-Sep 2025 (dense)", "2025-08-10", "2025-09-18"),
    ]

    fig, axes = plt.subplots(len(plot_windows), 1, figsize=(16, 4 * len(plot_windows)))

    for idx, (label, start, end) in enumerate(plot_windows):
        ax = axes[idx]
        mask = (d["date"] >= start) & (d["date"] <= end)
        sub = d[mask]
        w = sub["weight_lbs"].notna()

        ax.plot(sub.loc[w, "date"], sub.loc[w, "weight_lbs"],
                "o-", color="tab:blue", ms=3, lw=1, label="Observed", alpha=0.8)
        ax.plot(sub.loc[w, "date"], sub.loc[w, "smoothed"],
                "s-", color="tab:red", ms=3, lw=1, label="Smoothed", alpha=0.8)

        # Show carbs as background bars (scaled to fit)
        ax2 = ax.twinx()
        colors = ["tab:red" if c < 100 else "tab:orange" if c < 200 else "tab:green"
                  for c in sub["carbs_g"]]
        ax2.bar(sub["date"], sub["carbs_g"], color=colors, alpha=0.15, width=0.8)
        ax2.set_ylabel("Carbs (g)", fontsize=8, color="gray")
        ax2.set_ylim(0, 600)
        ax2.tick_params(labelsize=7, colors="gray")

        # Stats
        ev = eval_window(d, start, end)
        stats = f"var red: {ev['reduction']:.1f}%" if ev else "n/a"

        ax.set_ylabel("Weight (lbs)")
        ax.set_title(f"{label}  |  {stats}", fontsize=10)
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, alpha=0.2)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))

    plt.suptitle(f"Glycogen Smoothing: g_max={gm}  carb_ref={cr}  "
                 f"rate_up={ru:.2f}  rate_down={rd:.2f}", fontsize=12, y=1.01)
    plt.tight_layout()
    out = ROOT / "analysis" / "glycogen_multiwindow.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nSaved {out}")
    plt.close()


if __name__ == "__main__":
    main()
