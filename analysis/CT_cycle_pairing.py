"""CT. Pair adjacent falls and rises. Look for memory between cycles.

If the trajectory is composed of two coupled hidden values — call them SP (slow,
forms peaks) and W (willpower budget, forms troughs) — then we'd expect:

  - Fall depth (how far W could push FM below the recent peak) predicts the
    SIZE of the subsequent rise (W exhausts → FM bounces back).
  - Time spent at low FM (cumulative willpower-spend) predicts rise size.
  - The peak-trough symmetry per cycle: does the rise return to the previous
    peak? Or does each peak ratchet up regardless of trough depth?

No model fitting; just paired numbers and the scatter shape.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

ROOT = Path(__file__).resolve().parent.parent
INTAKE_MULTIPLIER = 1.09


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
    df = df.dropna(subset=["fm", "tdee", "calories"]).reset_index(drop=True)
    pre = df[df["effective_level"] == 0].copy().reset_index(drop=True)

    fm_smooth = pre["fm"].rolling(14, min_periods=7, center=True).mean().values

    pk_idx, _ = find_peaks(fm_smooth, prominence=1.5, distance=21)
    tr_idx, _ = find_peaks(-fm_smooth, prominence=1.5, distance=21)

    # Build alternating extrema sequence
    extrema = sorted(
        [(i, "peak", float(fm_smooth[i])) for i in pk_idx] +
        [(i, "trough", float(fm_smooth[i])) for i in tr_idx],
        key=lambda x: x[0])

    # Walk extrema, build (peak_a, trough, peak_b) triplets — one full cycle
    triplets = []
    for k in range(len(extrema) - 2):
        a, b, c = extrema[k], extrema[k + 1], extrema[k + 2]
        if a[1] == "peak" and b[1] == "trough" and c[1] == "peak":
            triplets.append({
                "peak_a_idx": a[0], "peak_a_fm": a[2], "peak_a_date": pre["date"].iloc[a[0]],
                "trough_idx": b[0], "trough_fm": b[2], "trough_date": pre["date"].iloc[b[0]],
                "peak_b_idx": c[0], "peak_b_fm": c[2], "peak_b_date": pre["date"].iloc[c[0]],
            })

    rows = []
    for t in triplets:
        fall = t["trough_fm"] - t["peak_a_fm"]            # negative
        rise = t["peak_b_fm"] - t["trough_fm"]            # positive
        net  = t["peak_b_fm"] - t["peak_a_fm"]            # net peak-to-peak
        fall_days = (t["trough_date"] - t["peak_a_date"]).days
        rise_days = (t["peak_b_date"] - t["trough_date"]).days
        # Cumulative cal during fall: trapezoidal sum of intake-tdee
        m_fall = (pre["date"] >= t["peak_a_date"]) & (pre["date"] <= t["trough_date"])
        m_rise = (pre["date"] >= t["trough_date"]) & (pre["date"] <= t["peak_b_date"])
        cum_fall = float((pre.loc[m_fall, "calories"] - pre.loc[m_fall, "tdee"]).sum())
        cum_rise = float((pre.loc[m_rise, "calories"] - pre.loc[m_rise, "tdee"]).sum())
        rows.append({
            "peak_a_date": t["peak_a_date"], "peak_a_fm": t["peak_a_fm"],
            "trough_date": t["trough_date"], "trough_fm": t["trough_fm"],
            "peak_b_date": t["peak_b_date"], "peak_b_fm": t["peak_b_fm"],
            "fall_lb": fall, "fall_days": fall_days, "cum_fall_cal": cum_fall,
            "rise_lb": rise, "rise_days": rise_days, "cum_rise_cal": cum_rise,
            "net_peak_lb": net,
        })

    cy = pd.DataFrame(rows)
    print("=" * 110)
    print(f"FULL CYCLES (peak → trough → next-peak): n = {len(cy)}")
    print("=" * 110)
    print(f"  {'peak_a':<11} {'trough':<11} {'peak_b':<11} {'fall':>6} {'rise':>6} "
          f"{'net':>6} {'cum_fall':>9} {'cum_rise':>9}")
    for _, r in cy.iterrows():
        print(f"  {r['peak_a_date'].strftime('%Y-%m-%d'):<11} "
              f"{r['trough_date'].strftime('%Y-%m-%d'):<11} "
              f"{r['peak_b_date'].strftime('%Y-%m-%d'):<11} "
              f"{r['fall_lb']:>+6.1f} {r['rise_lb']:>+6.1f} {r['net_peak_lb']:>+6.1f} "
              f"{r['cum_fall_cal']:>+9.0f} {r['cum_rise_cal']:>+9.0f}")

    print("\n" + "=" * 110)
    print("FALL → RISE PAIRING (within full cycles)")
    print("=" * 110)
    if len(cy) > 3:
        # The user's specific hypothesis: fall depth predicts subsequent rise magnitude
        r_fall_rise = float(np.corrcoef(cy["fall_lb"], cy["rise_lb"])[0, 1])
        r_fall_net  = float(np.corrcoef(cy["fall_lb"], cy["net_peak_lb"])[0, 1])
        r_cum_rise  = float(np.corrcoef(cy["cum_fall_cal"], cy["rise_lb"])[0, 1])
        # Time-spent-low (proxy for willpower spent)
        r_falldays_rise = float(np.corrcoef(cy["fall_days"], cy["rise_lb"])[0, 1])
        print(f"\n  Corr(fall depth, subsequent rise):       r = {r_fall_rise:+.3f}")
        print(f"  Corr(fall depth, net peak-to-peak):      r = {r_fall_net:+.3f}")
        print(f"  Corr(fall duration, subsequent rise):    r = {r_falldays_rise:+.3f}")
        print(f"  Corr(cum cal during fall, rise size):    r = {r_cum_rise:+.3f}")

        print("\n  Interpretation guide:")
        print("  - Strong positive r(fall depth, rise) → bigger falls produce bigger rises")
        print("    (consistent with willpower-exhaustion-and-rebound model)")
        print("  - r near zero → falls and rises are independent events")
        print("  - Negative r → bigger falls suppress subsequent rises")

    # ── Asymmetry over time: does fall/rise asymmetry change? ──
    print("\n" + "=" * 110)
    print("AMPLITUDE ASYMMETRY OVER TIME")
    print("=" * 110)
    if len(cy):
        cy["abs_fall"] = cy["fall_lb"].abs()
        cy["abs_rise"] = cy["rise_lb"]
        cy["asymmetry_lb"] = cy["abs_rise"] - cy["abs_fall"]
        cy["asymmetry_ratio"] = cy["abs_rise"] / cy["abs_fall"].replace(0, np.nan)
        print(f"\n  {'cycle peak_a':<12} {'|fall|':>6} {'|rise|':>6} {'asymm_lb':>8} {'asymm_ratio':>11}")
        for _, r in cy.iterrows():
            print(f"  {r['peak_a_date'].strftime('%Y-%m-%d'):<12} "
                  f"{r['abs_fall']:>6.1f} {r['abs_rise']:>6.1f} "
                  f"{r['asymmetry_lb']:>+8.2f} {r['asymmetry_ratio']:>11.2f}")

    # ── Plot ──
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    ax = axes[0, 0]
    if len(cy) > 0:
        ax.scatter(cy["fall_lb"].abs(), cy["rise_lb"], s=70, c=range(len(cy)),
                   cmap="viridis", edgecolor="black")
        for _, r in cy.iterrows():
            ax.annotate(r["peak_a_date"].strftime("%Y-%m"),
                        (abs(r["fall_lb"]), r["rise_lb"]), fontsize=7, alpha=0.7)
        ax.set_xlabel("|fall| (lb, peak_a → trough)")
        ax.set_ylabel("rise (lb, trough → peak_b)")
        ax.set_title(f"Fall depth vs subsequent rise  "
                     f"(r = {np.corrcoef(cy['fall_lb'].abs(), cy['rise_lb'])[0,1]:+.3f})")
        # 1:1 line
        m = max(cy["fall_lb"].abs().max(), cy["rise_lb"].max())
        ax.plot([0, m], [0, m], "k--", alpha=0.4, label="1:1")
        ax.legend()
        ax.grid(alpha=0.3)

    ax = axes[0, 1]
    if len(cy) > 0:
        ax.scatter(cy["fall_days"], cy["rise_lb"], s=70, c=range(len(cy)),
                   cmap="viridis", edgecolor="black")
        ax.set_xlabel("fall duration (days)")
        ax.set_ylabel("rise (lb)")
        ax.set_title(f"Fall duration vs rise  "
                     f"(r = {np.corrcoef(cy['fall_days'], cy['rise_lb'])[0,1]:+.3f})")
        ax.grid(alpha=0.3)

    ax = axes[1, 0]
    # Cycle amplitude over time
    if len(cy) > 0:
        ax.plot(cy["peak_a_date"], cy["abs_fall"], "o-", color="tab:blue",
                label="|fall|")
        ax.plot(cy["peak_a_date"], cy["abs_rise"], "o-", color="tab:red",
                label="|rise|")
        ax.set_ylabel("amplitude (lb)")
        ax.set_xlabel("Date")
        ax.set_title("Cycle amplitudes — falls vs rises")
        ax.legend()
        ax.grid(alpha=0.3)

    ax = axes[1, 1]
    # Net peak-to-peak movement
    if len(cy) > 0:
        ax.bar(range(len(cy)), cy["net_peak_lb"],
               color=["tab:red" if v > 0 else "tab:blue" for v in cy["net_peak_lb"]])
        ax.axhline(0, color="black", lw=0.5)
        ax.set_xticks(range(len(cy)))
        ax.set_xticklabels([d.strftime("%Y-%m") for d in cy["peak_a_date"]],
                            rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("net peak shift (lb)")
        ax.set_title("Net peak-to-peak FM shift per cycle")
        ax.grid(alpha=0.3)

    fig.tight_layout()
    plot_path = ROOT / "analysis" / "CT_cycle_pairing.png"
    fig.savefig(plot_path, dpi=140)
    print(f"\nSaved {plot_path}")

    cy.to_csv(ROOT / "analysis" / "CT_cycles.csv", index=False)


if __name__ == "__main__":
    main()
