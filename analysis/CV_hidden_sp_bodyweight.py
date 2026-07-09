"""CV. Re-run the hidden-SP mass-balance using total body weight (FM + LM).

CU's per-cycle SP_i clustered around 200-300 lb in FM-units, with recent
cycles (2022-2023) landing at 285-294. Subject's pre-data peak weight was
~285-300 lb. Hypothesis: SP is defended on TOTAL BODY WEIGHT, not FM.

Same mass balance as CU, with body weight replacing FM in the drain:

    α · ∫_cycle max(0, SP − BW) dt  =  β · ∫_cycle max(0, intake − TDEE) dt

If SP_BW pins down to a sharp value (say 285-300) with LOW SD across
cycles, that's strong evidence for a defended total-weight set point at
the subject's lifetime peak.

Convincing thresholds (preregistered):
  - sharp interior minimum in CV(SP_BW)
  - best CV < 0.5 with SP in plausible range (250-350 lb)
  - per-cycle SP_i clusters within SD < 15 lb
  - cluster center near 285-300 lb (the subject's documented peak)
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent.parent
INTAKE_MULTIPLIER = 1.09


def cycle_drain_bw(bw_arr, sp):
    return float(np.sum(np.maximum(0.0, sp - bw_arr)))


def cycle_recharge(surplus_arr):
    return float(np.sum(np.maximum(0.0, surplus_arr)))


def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    intake["calories"] = intake["calories"] * INTAKE_MULTIPLIER
    kf = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    cycles = pd.read_csv(ROOT / "analysis" / "CT_cycles.csv",
                          parse_dates=["peak_a_date", "trough_date", "peak_b_date"])

    df = intake[["date", "calories"]].merge(
        kf[["date", "fat_mass_lbs", "lm_lbs", "tdee"]], on="date")
    df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df = df.sort_values("date").reset_index(drop=True)
    df["fm"] = df["fat_mass_lbs"]
    df["lm"] = df["lm_lbs"]
    df["bw"] = df["fm"] + df["lm"]
    df["surplus"] = df["calories"] - df["tdee"]
    df = df.dropna(subset=["fm", "lm", "tdee", "calories"]).reset_index(drop=True)
    pre = df[df["effective_level"] == 0].copy().reset_index(drop=True)

    print("=" * 100)
    print("CV. Hidden SP from cycle mass balance, computed on TOTAL BODY WEIGHT")
    print("=" * 100)
    print(f"\nPre-tirz days: {len(pre)}")
    print(f"BW range: {pre['bw'].min():.1f} to {pre['bw'].max():.1f} lb")
    print(f"FM range: {pre['fm'].min():.1f} to {pre['fm'].max():.1f} lb")
    print(f"LM range: {pre['lm'].min():.1f} to {pre['lm'].max():.1f} lb")

    cycle_data = []
    for _, c in cycles.iterrows():
        m = (pre["date"] >= c["peak_a_date"]) & (pre["date"] <= c["peak_b_date"])
        g = pre[m]
        if len(g) < 30:
            continue
        cycle_data.append({
            "peak_a_date": c["peak_a_date"], "peak_b_date": c["peak_b_date"],
            "trough_date": c["trough_date"],
            "bw": g["bw"].values, "fm": g["fm"].values, "lm": g["lm"].values,
            "surplus": g["surplus"].values,
            "n": len(g),
            "bw_max": float(g["bw"].max()),
            "fm_max": float(g["fm"].max()),
            "recharge": cycle_recharge(g["surplus"].values),
        })
    print(f"\nUsing {len(cycle_data)} cycles")

    # ── Sweep SP_BW ──
    SP_grid = np.arange(150, 401, 1.0)
    rows = []
    for sp in SP_grid:
        ratios = []
        for c in cycle_data:
            d = cycle_drain_bw(c["bw"], sp)
            r = c["recharge"]
            if d > 0 and r > 0:
                ratios.append(d / r)
        if len(ratios) < 5:
            continue
        ra = np.array(ratios)
        rows.append({
            "SP_BW": sp,
            "mean_R": float(ra.mean()),
            "cv_R": float(ra.std() / ra.mean()) if ra.mean() > 0 else float("inf"),
            "log_std_R": float(np.std(np.log(ra))),
        })
    sweep = pd.DataFrame(rows)

    best_idx = sweep["cv_R"].idxmin()
    best_sp = float(sweep.loc[best_idx, "SP_BW"])
    best_cv = float(sweep.loc[best_idx, "cv_R"])
    best_logstd_idx = sweep["log_std_R"].idxmin()
    best_sp_log = float(sweep.loc[best_logstd_idx, "SP_BW"])

    print(f"\nBest SP_BW by CV(R):     {best_sp:.0f} lb  (CV = {best_cv:.3f})")
    print(f"Best SP_BW by SD(log R): {best_sp_log:.0f} lb  "
          f"(SD log R = {sweep.loc[best_logstd_idx, 'log_std_R']:.3f})")
    print(f"\nCV(R) at SP=200:  {sweep[sweep['SP_BW']==200]['cv_R'].values[0]:.3f}")
    print(f"CV(R) at SP=250:  {sweep[sweep['SP_BW']==250]['cv_R'].values[0]:.3f}")
    print(f"CV(R) at SP=285:  {sweep[sweep['SP_BW']==285]['cv_R'].values[0]:.3f}")
    print(f"CV(R) at SP=300:  {sweep[sweep['SP_BW']==300]['cv_R'].values[0]:.3f}")
    print(f"CV(R) at SP=350:  {sweep[sweep['SP_BW']==350]['cv_R'].values[0]:.3f}")
    print(f"CV(R) at SP=400:  {sweep[sweep['SP_BW']==400]['cv_R'].values[0]:.3f}")

    print("\nPreregistered checks:")
    sharp_min = (SP_grid[0] + 5 < best_sp < SP_grid[-1] - 5)
    print(f"  [{('PASS' if sharp_min else 'FAIL')}] sharp interior minimum: best={best_sp:.0f}")
    print(f"  [{('PASS' if best_cv < 0.5 else 'FAIL')}] best CV < 0.5: {best_cv:.3f}")
    in_range = 250 <= best_sp <= 350
    print(f"  [{('PASS' if in_range else 'FAIL')}] best SP in 250-350 plausible band: "
          f"{best_sp:.0f}")

    # ── Per-cycle free SP_i_BW ──
    kappa = float(sweep.loc[best_idx, "mean_R"])
    print(f"\n  Using κ = mean_R(SP={best_sp:.0f}) = {kappa:.4f}\n")

    print(f"  {'cycle (peak_a)':<14} {'SP_i_BW':>8} {'BW max':>7} {'drain':>9} {'recharge':>9}")
    sp_per_cycle = []
    for c in cycle_data:
        target_drain = kappa * c["recharge"]

        def f(sp):
            return cycle_drain_bw(c["bw"], sp) - target_drain

        sp_lo = c["bw_max"] - 2.0
        sp_hi = 1000.0
        try:
            if f(sp_lo) < 0 < f(sp_hi):
                sp_i = brentq(f, sp_lo, sp_hi)
            else:
                sp_i = float("nan")
        except Exception:
            sp_i = float("nan")
        sp_per_cycle.append({
            "peak_a_date": pd.Timestamp(c["peak_a_date"]),
            "trough_date": pd.Timestamp(c["trough_date"]),
            "peak_b_date": pd.Timestamp(c["peak_b_date"]),
            "sp_i_bw": sp_i,
            "bw_max": c["bw_max"],
            "fm_max": c["fm_max"],
        })
        print(f"  {pd.Timestamp(c['peak_a_date']).strftime('%Y-%m-%d'):<14} "
              f"{sp_i:>8.1f} {c['bw_max']:>7.1f} {cycle_drain_bw(c['bw'], sp_i):>9.0f} "
              f"{c['recharge']:>9.0f}")

    sp_df = pd.DataFrame(sp_per_cycle).dropna(subset=["sp_i_bw"])
    print(f"\n  SP_i_BW: mean = {sp_df['sp_i_bw'].mean():.1f} lb, "
          f"median = {sp_df['sp_i_bw'].median():.1f}, "
          f"SD = {sp_df['sp_i_bw'].std():.2f}, "
          f"min = {sp_df['sp_i_bw'].min():.1f}, max = {sp_df['sp_i_bw'].max():.1f}")
    print(f"  IQR: {sp_df['sp_i_bw'].quantile(0.25):.1f} to "
          f"{sp_df['sp_i_bw'].quantile(0.75):.1f}")
    if len(sp_df) > 5:
        r_drift = float(np.corrcoef(sp_df["bw_max"], sp_df["sp_i_bw"])[0, 1])
        print(f"  Corr(SP_i_BW, BW_max in cycle) = {r_drift:+.3f}")

    sp_clustered = bool(sp_df["sp_i_bw"].std() < 15)
    print(f"  [{('PASS' if sp_clustered else 'FAIL')}] per-cycle SP_i SD < 15 lb: "
          f"{sp_df['sp_i_bw'].std():.2f}")

    near_peak = bool(280 <= sp_df["sp_i_bw"].median() <= 305)
    print(f"  [{('PASS' if near_peak else 'FAIL')}] median SP_i_BW near subject peak "
          f"(285-300): {sp_df['sp_i_bw'].median():.1f}")

    # ── Subset: post-2014 cycles only (after the willpower-era anomaly) ──
    print("\n" + "=" * 100)
    print("POST-2014 CYCLES ONLY (excluding willpower-era anomaly)")
    print("=" * 100)
    post_cycles = [c for c in cycle_data
                    if pd.Timestamp(c["peak_a_date"]) >= pd.Timestamp("2015-01-01")]
    rows2 = []
    for sp in SP_grid:
        ratios = []
        for c in post_cycles:
            d = cycle_drain_bw(c["bw"], sp)
            r = c["recharge"]
            if d > 0 and r > 0:
                ratios.append(d / r)
        if len(ratios) < 5:
            continue
        ra = np.array(ratios)
        rows2.append({
            "SP_BW": sp,
            "mean_R": float(ra.mean()),
            "cv_R": float(ra.std() / ra.mean()) if ra.mean() > 0 else float("inf"),
        })
    sweep2 = pd.DataFrame(rows2)
    if len(sweep2):
        best2 = sweep2.loc[sweep2["cv_R"].idxmin()]
        print(f"\n  Best SP_BW (post-2014 only): {best2['SP_BW']:.0f} lb, CV = {best2['cv_R']:.3f}")

        # Per-cycle SP_i for post-2014
        kappa2 = float(best2["mean_R"])
        sp_post = []
        for c in post_cycles:
            target = kappa2 * c["recharge"]

            def f2(sp):
                return cycle_drain_bw(c["bw"], sp) - target

            sp_lo = c["bw_max"] - 2.0
            try:
                if f2(sp_lo) < 0 < f2(1000.0):
                    sp_post.append(brentq(f2, sp_lo, 1000.0))
            except Exception:
                pass
        if sp_post:
            sp_post = np.array(sp_post)
            print(f"  Per-cycle SP_i (post-2014, n={len(sp_post)}): "
                  f"mean = {sp_post.mean():.1f}, median = {np.median(sp_post):.1f}, "
                  f"SD = {sp_post.std():.2f}, IQR = {np.percentile(sp_post, 25):.1f} to "
                  f"{np.percentile(sp_post, 75):.1f}")

    # ── Plot ──
    fig, axes = plt.subplots(3, 1, figsize=(14, 11))

    ax = axes[0]
    ax.plot(sweep["SP_BW"], sweep["cv_R"], color="tab:blue", lw=1.5, label="all cycles")
    if len(sweep2):
        ax.plot(sweep2["SP_BW"], sweep2["cv_R"], color="tab:green", lw=1.5,
                label="post-2014 only")
    ax.axvspan(285, 300, color="tab:red", alpha=0.15, label="subject lifetime peak 285-300")
    ax.axvline(best_sp, color="tab:red", ls="--", label=f"min CV at SP={best_sp:.0f}")
    ax.axhline(0.5, color="grey", ls=":", label="prereg CV=0.5")
    ax.set_xlabel("Candidate SP (total body weight, lb)")
    ax.set_ylabel("CV of R = drain/recharge across cycles")
    ax.set_title(f"Mass-balance scan on body weight. Best CV = {best_cv:.3f} at SP = {best_sp:.0f} lb")
    ax.legend()
    ax.grid(alpha=0.3)

    ax = axes[1]
    if len(sp_df):
        ax.plot(sp_df["peak_a_date"], sp_df["sp_i_bw"], "o-", color="tab:purple",
                label="SP_i_BW (per cycle)")
        ax.axhline(best_sp, color="tab:red", ls="--", label=f"global SP={best_sp:.0f}")
        ax.axhspan(285, 300, color="tab:red", alpha=0.15, label="subject peak 285-300")
        ax.plot(sp_df["peak_a_date"], sp_df["bw_max"], "s-", color="grey",
                alpha=0.5, label="BW max in cycle")
        ax.set_ylabel("body weight (lb)")
        ax.set_xlabel("Cycle date")
        ax.set_title(f"Per-cycle SP_i_BW: SD = {sp_df['sp_i_bw'].std():.2f} lb, "
                     f"corr with BW_max = {np.corrcoef(sp_df['bw_max'], sp_df['sp_i_bw'])[0,1]:+.2f}")
        ax.legend(loc="upper left")
        ax.grid(alpha=0.3)

    # Histogram of SP_i_BW values
    ax = axes[2]
    if len(sp_df):
        ax.hist(sp_df["sp_i_bw"], bins=15, color="tab:purple", alpha=0.7,
                edgecolor="black")
        ax.axvspan(285, 300, color="tab:red", alpha=0.15, label="subject peak")
        ax.axvline(sp_df["sp_i_bw"].median(), color="black", lw=2,
                   label=f"median = {sp_df['sp_i_bw'].median():.0f}")
        ax.set_xlabel("SP_i_BW (lb)")
        ax.set_ylabel("count")
        ax.set_title("Distribution of per-cycle SP_i (body weight)")
        ax.legend()
        ax.grid(alpha=0.3)

    fig.tight_layout()
    plot_path = ROOT / "analysis" / "CV_hidden_sp_bodyweight.png"
    fig.savefig(plot_path, dpi=140)
    print(f"\nSaved {plot_path}")

    sweep.to_csv(ROOT / "analysis" / "CV_sp_sweep_bw.csv", index=False)
    sp_df.to_csv(ROOT / "analysis" / "CV_sp_per_cycle_bw.csv", index=False)


if __name__ == "__main__":
    main()
