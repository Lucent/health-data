"""CX. Test if K hysteresis drives the trough/peak cycle structure.

K's direction-dependent TDEE: falling phase has TDEE +124 cal/day, rising
phase −90 cal/day at matched FM. Spread = 214 cal/day.

Hypothesis: this directional asymmetry IS the cycle mechanism, not willpower
or SP. Predictions:

  H1. For each cycle, mean TDEE/RMR during the FALL > during the RISE.
      Sign-test should pass on 19/19 (or near). Magnitude ~5-10% of RMR.

  H2. At each trough, TDEE/RMR drops from "fall-side" to "rise-side".
      Compute TDEE/RMR mean over [trough-30d, trough] vs [trough, trough+30d].
      Pre > post.

  H3. At each peak, TDEE/RMR rises from "rise-side" to "fall-side".
      Pre < post (the next fall begins).

  H4. Cycles where fall-TDEE/RMR is much higher than rise-TDEE/RMR (large
      local K) should have bigger rises (per fall amplitude).

  H5. Plateau-mean TDEE/RMR (no clear direction) should sit between fall
      and rise values — the average across direction-mixed periods.

If H1-H4 hold cleanly, K's directional metabolism is the load-bearing
cycle mechanism. The trough is where intake matches the now-compressed
TDEE; the rebound is the K-asymmetric momentum continuing in the new
direction.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
INTAKE_MULTIPLIER = 1.09


def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    intake["calories"] = intake["calories"] * INTAKE_MULTIPLIER
    kf = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    p2 = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    cycles = pd.read_csv(ROOT / "analysis" / "CT_cycles.csv",
                          parse_dates=["peak_a_date", "trough_date", "peak_b_date"])

    df = intake[["date", "calories"]].merge(
        kf[["date", "fat_mass_lbs", "lm_lbs", "tdee"]], on="date")
    df = df.merge(p2[["date", "expected_rmr"]], on="date", how="left")
    df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df = df.sort_values("date").reset_index(drop=True)
    df["fm"] = df["fat_mass_lbs"]
    df["lm"] = df["lm_lbs"]
    df["bw"] = df["fm"] + df["lm"]
    df["ratio"] = df["tdee"] / df["expected_rmr"]
    df = df.dropna(subset=["fm", "tdee", "expected_rmr", "ratio"]).reset_index(drop=True)
    pre = df[df["effective_level"] == 0].copy().reset_index(drop=True)

    # 14-day centered smooth of TDEE/RMR ratio
    pre["ratio_14"] = pre["ratio"].rolling(14, min_periods=7, center=True).mean()

    print("=" * 100)
    print("CX. K hysteresis as cycle mechanism")
    print("=" * 100)

    # ── Per-cycle fall vs rise TDEE/RMR ──
    rows = []
    for _, c in cycles.iterrows():
        m_fall = (pre["date"] >= c["peak_a_date"]) & (pre["date"] <= c["trough_date"])
        m_rise = (pre["date"] >= c["trough_date"]) & (pre["date"] <= c["peak_b_date"])
        gf = pre[m_fall]
        gr = pre[m_rise]
        if len(gf) < 7 or len(gr) < 7:
            continue
        fall_ratio = float(gf["ratio"].mean())
        rise_ratio = float(gr["ratio"].mean())
        fall_tdee = float(gf["tdee"].mean())
        rise_tdee = float(gr["tdee"].mean())
        fall_fm_mean = float(gf["fm"].mean())
        rise_fm_mean = float(gr["fm"].mean())
        # Local K (TDEE difference, direction-of-change matched)
        k_local = fall_tdee - rise_tdee

        # Rise rate (lb/day)
        rise_amp = float(gr["fm"].iloc[-1] - gr["fm"].iloc[0])
        rise_days = max((c["peak_b_date"] - c["trough_date"]).days, 1)
        rise_rate = rise_amp / rise_days

        rows.append({
            "peak_a_date": c["peak_a_date"],
            "trough_date": c["trough_date"],
            "peak_b_date": c["peak_b_date"],
            "fall_ratio": fall_ratio,
            "rise_ratio": rise_ratio,
            "delta_ratio": fall_ratio - rise_ratio,
            "fall_tdee": fall_tdee,
            "rise_tdee": rise_tdee,
            "k_local_cal": k_local,
            "fall_fm_mean": fall_fm_mean,
            "rise_fm_mean": rise_fm_mean,
            "rise_amp_lb": rise_amp,
            "rise_rate_lb_per_d": rise_rate,
            "rise_days": rise_days,
        })

    cy = pd.DataFrame(rows)
    print(f"\nCycles analyzed: {len(cy)}")

    # ── H1: fall-ratio > rise-ratio per cycle ──
    print("\n" + "=" * 100)
    print("H1: per-cycle fall-side TDEE/RMR > rise-side TDEE/RMR?")
    print("=" * 100)
    print(f"  {'cycle':<14} {'fall ratio':>11} {'rise ratio':>11} "
          f"{'Δ ratio':>9} {'K local cal':>11}")
    for _, r in cy.iterrows():
        print(f"  {r['peak_a_date'].strftime('%Y-%m-%d'):<14} "
              f"{r['fall_ratio']:>11.4f} {r['rise_ratio']:>11.4f} "
              f"{r['delta_ratio']:>+9.4f} {r['k_local_cal']:>+11.0f}")

    n_pos = int((cy["delta_ratio"] > 0).sum())
    n_total = len(cy)
    print(f"\n  Cycles with fall > rise: {n_pos}/{n_total} ({100*n_pos/n_total:.0f}%)")
    print(f"  Median fall-rise Δ ratio: {cy['delta_ratio'].median():+.4f}")
    print(f"  Median K local (cal/day): {cy['k_local_cal'].median():+.0f}")
    print(f"  Compare to global K=+214 cal/day from FINDINGS")

    # Sign test against null of 50/50
    from scipy.stats import binomtest
    sign_test = binomtest(n_pos, n_total, 0.5, alternative='greater')
    print(f"  Sign test p-value: {sign_test.pvalue:.4f}")

    # ── H2 / H3: at trough, ratio drops; at peak, ratio rises ──
    print("\n" + "=" * 100)
    print("H2/H3: TDEE/RMR pre vs post each trough and peak")
    print("=" * 100)
    print(f"\n  TROUGHS ({len(cy)} of them):")
    print(f"  {'date':<12} {'pre30 ratio':>11} {'post30 ratio':>13} {'Δ':>8}")
    trough_deltas = []
    for _, r in cy.iterrows():
        d = r["trough_date"]
        pre_window = (pre["date"] >= d - pd.Timedelta(days=30)) & (pre["date"] < d)
        post_window = (pre["date"] > d) & (pre["date"] <= d + pd.Timedelta(days=30))
        ratio_pre = float(pre[pre_window]["ratio"].mean()) if pre_window.sum() > 0 else float("nan")
        ratio_post = float(pre[post_window]["ratio"].mean()) if post_window.sum() > 0 else float("nan")
        trough_deltas.append(ratio_pre - ratio_post)
        print(f"  {d.strftime('%Y-%m-%d'):<12} {ratio_pre:>11.4f} {ratio_post:>13.4f} "
              f"{ratio_pre - ratio_post:>+8.4f}")

    td = np.array(trough_deltas)
    td = td[~np.isnan(td)]
    n_trough_drop = int(np.sum(td > 0))
    print(f"\n  Troughs with ratio dropping (pre > post): {n_trough_drop}/{len(td)} "
          f"({100*n_trough_drop/len(td):.0f}%)")
    print(f"  Median Δ ratio at trough (pre30 - post30): {np.median(td):+.4f}")

    print(f"\n  PEAKS_b ({len(cy)} of them):")
    print(f"  {'date':<12} {'pre30 ratio':>11} {'post30 ratio':>13} {'Δ':>8}")
    peak_deltas = []
    for _, r in cy.iterrows():
        d = r["peak_b_date"]
        pre_window = (pre["date"] >= d - pd.Timedelta(days=30)) & (pre["date"] < d)
        post_window = (pre["date"] > d) & (pre["date"] <= d + pd.Timedelta(days=30))
        ratio_pre = float(pre[pre_window]["ratio"].mean()) if pre_window.sum() > 0 else float("nan")
        ratio_post = float(pre[post_window]["ratio"].mean()) if post_window.sum() > 0 else float("nan")
        peak_deltas.append(ratio_post - ratio_pre)
        print(f"  {d.strftime('%Y-%m-%d'):<12} {ratio_pre:>11.4f} {ratio_post:>13.4f} "
              f"{ratio_post - ratio_pre:>+8.4f}")

    pd_arr = np.array(peak_deltas)
    pd_arr = pd_arr[~np.isnan(pd_arr)]
    n_peak_rise = int(np.sum(pd_arr > 0))
    print(f"\n  Peaks with ratio rising (post > pre): {n_peak_rise}/{len(pd_arr)} "
          f"({100*n_peak_rise/len(pd_arr):.0f}%)")
    print(f"  Median Δ ratio at peak (post30 - pre30): {np.median(pd_arr):+.4f}")

    # ── H4: K asymmetry magnitude predicts rise rate ──
    print("\n" + "=" * 100)
    print("H4: does fall-rise TDEE Δ predict rise rate?")
    print("=" * 100)
    if len(cy) > 5:
        r_kr = float(np.corrcoef(cy["k_local_cal"], cy["rise_rate_lb_per_d"])[0, 1])
        r_kamp = float(np.corrcoef(cy["k_local_cal"], cy["rise_amp_lb"])[0, 1])
        r_dr_amp = float(np.corrcoef(cy["delta_ratio"], cy["rise_amp_lb"])[0, 1])
        print(f"  Corr(K_local_cal, rise_rate_lb_per_d): r = {r_kr:+.3f}")
        print(f"  Corr(K_local_cal, rise_amp_lb):        r = {r_kamp:+.3f}")
        print(f"  Corr(Δ ratio, rise_amp_lb):            r = {r_dr_amp:+.3f}")
        if r_kamp > 0.3:
            print("  → bigger fall-rise TDEE asymmetry predicts bigger rebound")
        elif r_kamp < -0.3:
            print("  → larger asymmetry → smaller rebound (opposite of prediction)")
        else:
            print("  → no clear relationship")

    # ── H5: plateau ratio is between fall and rise ──
    print("\n" + "=" * 100)
    print("H5: plateau-mean TDEE/RMR sits between fall-mean and rise-mean")
    print("=" * 100)
    plateau_ratios = []
    # Use the long stable plateau periods we know about
    plateau_periods = [
        ("Plateau-2 enlightened", "2015-07-01", "2017-12-31"),
        ("Plateau-3 post-keto", "2020-01-01", "2021-12-31"),
        ("Plateau-4 pre-tirz", "2023-12-01", "2024-04-30"),
    ]
    for name, start, end in plateau_periods:
        m = (pre["date"] >= start) & (pre["date"] <= end)
        g = pre[m]
        if len(g) > 30:
            r = float(g["ratio"].mean())
            plateau_ratios.append((name, r))
            print(f"  {name}: plateau ratio = {r:.4f}")

    print(f"\n  Cycle median fall ratio: {cy['fall_ratio'].median():.4f}")
    print(f"  Cycle median rise ratio: {cy['rise_ratio'].median():.4f}")

    # ── Plot ──
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    ax = axes[0, 0]
    x = np.arange(len(cy))
    width = 0.4
    ax.bar(x - width/2, cy["fall_ratio"], width, color="tab:blue", label="fall TDEE/RMR")
    ax.bar(x + width/2, cy["rise_ratio"], width, color="tab:red", label="rise TDEE/RMR")
    ax.set_xticks(x)
    ax.set_xticklabels([d.strftime("%Y-%m") for d in cy["peak_a_date"]],
                        rotation=45, ha="right", fontsize=7)
    ax.axhline(cy["fall_ratio"].median(), color="tab:blue", ls="--", alpha=0.5)
    ax.axhline(cy["rise_ratio"].median(), color="tab:red", ls="--", alpha=0.5)
    ax.set_ylabel("TDEE/RMR")
    ax.set_title(f"H1: per-cycle fall vs rise TDEE/RMR  ({n_pos}/{n_total} fall>rise)")
    ax.legend()
    ax.grid(alpha=0.3, axis="y")

    ax = axes[0, 1]
    ax.scatter(cy["k_local_cal"], cy["rise_amp_lb"], s=70,
                c=range(len(cy)), cmap="viridis", edgecolor="black")
    for _, r in cy.iterrows():
        ax.annotate(r["peak_a_date"].strftime("%Y-%m"),
                    (r["k_local_cal"], r["rise_amp_lb"]), fontsize=7, alpha=0.6)
    ax.set_xlabel("K local (fall TDEE − rise TDEE, cal/day)")
    ax.set_ylabel("rise amplitude (lb)")
    if len(cy) > 5:
        r_kamp = float(np.corrcoef(cy["k_local_cal"], cy["rise_amp_lb"])[0, 1])
        ax.set_title(f"H4: K asymmetry vs rise amplitude (r = {r_kamp:+.3f})")
    ax.axvline(0, color="black", lw=0.5)
    ax.grid(alpha=0.3)

    ax = axes[1, 0]
    ax.plot(pre["date"], pre["ratio_14"], color="black", lw=0.8, label="TDEE/RMR (14d)")
    for _, r in cy.iterrows():
        ax.axvline(r["trough_date"], color="tab:blue", alpha=0.2)
        ax.axvline(r["peak_b_date"], color="tab:red", alpha=0.2)
    ax.set_xlabel("Date")
    ax.set_ylabel("TDEE/RMR")
    ax.set_title("TDEE/RMR over time, troughs (blue) and peaks (red)")
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    # Plot trough_pre - trough_post and peak_post - peak_pre as histograms
    ax.hist(td, bins=12, alpha=0.6, color="tab:blue",
            label=f"trough (pre30 − post30), median={np.median(td):+.4f}")
    ax.hist(pd_arr, bins=12, alpha=0.6, color="tab:red",
            label=f"peak (post30 − pre30), median={np.median(pd_arr):+.4f}")
    ax.axvline(0, color="black", lw=0.5)
    ax.set_xlabel("Δ ratio (pre to post extremum)")
    ax.set_ylabel("count")
    ax.set_title("H2/H3: TDEE/RMR transitions at extrema")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    plot_path = ROOT / "analysis" / "CX_k_drives_cycles.png"
    fig.savefig(plot_path, dpi=140)
    print(f"\nSaved {plot_path}")

    cy.to_csv(ROOT / "analysis" / "CX_cycle_k.csv", index=False)


if __name__ == "__main__":
    main()
