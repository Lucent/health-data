"""CS. Decompose the FM trajectory into peaks and troughs.

Each cycle has a peak (local FM maximum) and a trough (local FM minimum).
If two hidden values control the trajectory, the peaks and troughs may
move differently — for example, one could be smooth (a slow-moving
defended weight) while the other reflects an oscillating willpower budget.

This script extracts the peak series and the trough series and shows them
side by side. No model. Just looking at the shapes separately.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

ROOT = Path(__file__).resolve().parent.parent
INTAKE_MULTIPLIER = 1.09
BINGE_THRESHOLD = 1000


def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    intake["calories"] = intake["calories"] * INTAKE_MULTIPLIER
    kf = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    df = intake[["date", "calories"]].merge(
        kf[["date", "fat_mass_lbs", "tdee"]], on="date")
    df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df = df.sort_values("date").reset_index(drop=True)
    df["fm"] = df["fat_mass_lbs"]
    df["surplus"] = df["calories"] - df["tdee"]
    df["binge"] = (df["calories"] > df["tdee"] + BINGE_THRESHOLD).astype(int)
    df = df.dropna(subset=["fm", "tdee", "calories"]).reset_index(drop=True)
    pre = df[df["effective_level"] == 0].copy().reset_index(drop=True)

    # Smooth FM lightly to reduce day-noise but preserve cycle structure
    fm_smooth = pre["fm"].rolling(14, min_periods=7, center=True).mean().values

    # Find peaks (local maxima) and troughs (local minima)
    # Use prominence to filter spurious wobbles
    PROM = 1.5  # lb of prominence
    DIST = 21   # min days between peaks
    pk_idx, pk_props = find_peaks(fm_smooth, prominence=PROM, distance=DIST)
    tr_idx, tr_props = find_peaks(-fm_smooth, prominence=PROM, distance=DIST)

    peaks = pre.iloc[pk_idx][["date", "fm"]].copy().reset_index(drop=True)
    peaks["fm_smooth"] = fm_smooth[pk_idx]
    peaks["prominence"] = pk_props["prominences"]
    troughs = pre.iloc[tr_idx][["date", "fm"]].copy().reset_index(drop=True)
    troughs["fm_smooth"] = fm_smooth[tr_idx]
    troughs["prominence"] = tr_props["prominences"]

    print("=" * 100)
    print(f"PEAKS (local maxima, prominence > {PROM} lb, ≥ {DIST} d apart): n = {len(peaks)}")
    print("=" * 100)
    for _, r in peaks.iterrows():
        print(f"  {r['date'].strftime('%Y-%m-%d')}  FM = {r['fm']:>5.1f} lb  "
              f"(smoothed {r['fm_smooth']:>5.1f}, prom {r['prominence']:>4.1f})")

    print("\n" + "=" * 100)
    print(f"TROUGHS (local minima, prominence > {PROM} lb, ≥ {DIST} d apart): n = {len(troughs)}")
    print("=" * 100)
    for _, r in troughs.iterrows():
        print(f"  {r['date'].strftime('%Y-%m-%d')}  FM = {r['fm']:>5.1f} lb  "
              f"(smoothed {r['fm_smooth']:>5.1f}, prom {r['prominence']:>4.1f})")

    # Pair peaks and troughs into cycles. Use trough → next-peak as a "rise" cycle
    # and peak → next-trough as a "fall" cycle.
    print("\n" + "=" * 100)
    print("CYCLES — paired (trough, next-peak) rises:")
    print("=" * 100)
    print(f"  {'trough date':<12} {'trough FM':>9} {'peak date':<12} {'peak FM':>8} "
          f"{'rise':>5} {'days':>5} {'rate':>9}")
    rises = []
    for _, t in troughs.iterrows():
        future_peaks = peaks[peaks["date"] > t["date"]]
        if len(future_peaks) == 0:
            continue
        p = future_peaks.iloc[0]
        rise = float(p["fm_smooth"] - t["fm_smooth"])
        days = (p["date"] - t["date"]).days
        rate = rise / max(days, 1)
        rises.append({"t_date": t["date"], "t_fm": t["fm_smooth"],
                       "p_date": p["date"], "p_fm": p["fm_smooth"],
                       "rise": rise, "days": days, "rate": rate})
        print(f"  {t['date'].strftime('%Y-%m-%d'):<12} {t['fm_smooth']:>9.1f} "
              f"{p['date'].strftime('%Y-%m-%d'):<12} {p['fm_smooth']:>8.1f} "
              f"{rise:>+5.1f} {days:>5d} {rate:>+9.4f}")

    print("\nCYCLES — paired (peak, next-trough) falls:")
    print(f"  {'peak date':<12} {'peak FM':>8} {'trough date':<12} {'trough FM':>9} "
          f"{'fall':>5} {'days':>5} {'rate':>9}")
    falls = []
    for _, p in peaks.iterrows():
        future_troughs = troughs[troughs["date"] > p["date"]]
        if len(future_troughs) == 0:
            continue
        t = future_troughs.iloc[0]
        fall = float(t["fm_smooth"] - p["fm_smooth"])
        days = (t["date"] - p["date"]).days
        rate = fall / max(days, 1)
        falls.append({"p_date": p["date"], "p_fm": p["fm_smooth"],
                       "t_date": t["date"], "t_fm": t["fm_smooth"],
                       "fall": fall, "days": days, "rate": rate})
        print(f"  {p['date'].strftime('%Y-%m-%d'):<12} {p['fm_smooth']:>8.1f} "
              f"{t['date'].strftime('%Y-%m-%d'):<12} {t['fm_smooth']:>9.1f} "
              f"{fall:>+5.1f} {days:>5d} {rate:>+9.4f}")

    rises_df = pd.DataFrame(rises)
    falls_df = pd.DataFrame(falls)

    # ── Are peaks smoother than troughs (or vice versa)? ──
    # Lag-1 difference SD of each series
    if len(peaks) > 2:
        pk_diffs = peaks["fm_smooth"].diff().dropna().values
        pk_sd = float(np.std(pk_diffs))
        pk_mean = float(np.mean(pk_diffs))
    else:
        pk_sd = pk_mean = float("nan")
    if len(troughs) > 2:
        tr_diffs = troughs["fm_smooth"].diff().dropna().values
        tr_sd = float(np.std(tr_diffs))
        tr_mean = float(np.mean(tr_diffs))
    else:
        tr_sd = tr_mean = float("nan")

    print("\n" + "=" * 100)
    print("PEAK-TO-PEAK vs TROUGH-TO-TROUGH dynamics")
    print("=" * 100)
    print(f"  Peaks:   n={len(peaks)}, mean Δpeak = {pk_mean:+.2f} lb, "
          f"SD of Δpeak = {pk_sd:.2f} lb")
    print(f"  Troughs: n={len(troughs)}, mean Δtrough = {tr_mean:+.2f} lb, "
          f"SD of Δtrough = {tr_sd:.2f} lb")
    print("\n  If peaks are smoother (lower SD): peaks track a slow-moving variable,")
    print("  troughs reflect a noisy/oscillating one.")
    print("  If both are equally noisy: cycles are roughly symmetric oscillation.")

    # Cycle amplitudes over time
    if rises_df.shape[0] > 0:
        print(f"\n  Rise amplitudes: median = {rises_df['rise'].median():+.2f} lb, "
              f"IQR = {rises_df['rise'].quantile(.25):+.2f} to {rises_df['rise'].quantile(.75):+.2f}")
        print(f"  Rise durations:  median = {rises_df['days'].median():.0f} d, "
              f"IQR = {rises_df['days'].quantile(.25):.0f} to {rises_df['days'].quantile(.75):.0f}")
    if falls_df.shape[0] > 0:
        print(f"  Fall amplitudes: median = {falls_df['fall'].median():+.2f} lb, "
              f"IQR = {falls_df['fall'].quantile(.25):+.2f} to {falls_df['fall'].quantile(.75):+.2f}")
        print(f"  Fall durations:  median = {falls_df['days'].median():.0f} d, "
              f"IQR = {falls_df['days'].quantile(.25):.0f} to {falls_df['days'].quantile(.75):.0f}")

    # ── Plot ──
    fig, axes = plt.subplots(3, 1, figsize=(15, 11), sharex=True)

    ax = axes[0]
    ax.plot(pre["date"], pre["fm"], color="grey", lw=0.6, alpha=0.6, label="FM (raw)")
    ax.plot(pre["date"], fm_smooth, color="black", lw=1.0, label="FM (14d smooth)")
    ax.plot(peaks["date"], peaks["fm_smooth"], "o", color="tab:red", ms=8,
            label="peaks", zorder=5)
    ax.plot(troughs["date"], troughs["fm_smooth"], "o", color="tab:blue", ms=8,
            label="troughs", zorder=5)
    # Connect peaks with a line, troughs with a line
    if len(peaks) > 1:
        ax.plot(peaks["date"], peaks["fm_smooth"], "-", color="tab:red",
                lw=1.2, alpha=0.6)
    if len(troughs) > 1:
        ax.plot(troughs["date"], troughs["fm_smooth"], "-", color="tab:blue",
                lw=1.2, alpha=0.6)
    ax.set_ylabel("FM (lb)")
    ax.set_title("FM trajectory with extracted peaks (red) and troughs (blue)")
    ax.legend(loc="upper left")
    ax.grid(alpha=0.3)

    # Cycle amplitudes vs time
    ax = axes[1]
    if rises_df.shape[0] > 0:
        ax.plot(rises_df["t_date"], rises_df["rise"], "o-", color="tab:red",
                label="rise amplitude (trough → peak)")
    if falls_df.shape[0] > 0:
        ax.plot(falls_df["p_date"], falls_df["fall"], "o-", color="tab:blue",
                label="fall amplitude (peak → trough)")
    ax.axhline(0, color="black", lw=0.5)
    ax.set_ylabel("amplitude (lb)")
    ax.set_title("Cycle amplitudes over time")
    ax.legend()
    ax.grid(alpha=0.3)

    # Cycle durations vs time
    ax = axes[2]
    if rises_df.shape[0] > 0:
        ax.plot(rises_df["t_date"], rises_df["days"], "o-", color="tab:red",
                label="rise duration")
    if falls_df.shape[0] > 0:
        ax.plot(falls_df["p_date"], falls_df["days"], "o-", color="tab:blue",
                label="fall duration")
    ax.set_ylabel("days")
    ax.set_xlabel("Date")
    ax.set_title("Cycle durations over time")
    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout()
    plot_path = ROOT / "analysis" / "CS_cycle_peaks_troughs.png"
    fig.savefig(plot_path, dpi=140)
    print(f"\nSaved {plot_path}")

    rises_df.to_csv(ROOT / "analysis" / "CS_rises.csv", index=False)
    falls_df.to_csv(ROOT / "analysis" / "CS_falls.csv", index=False)


if __name__ == "__main__":
    main()
