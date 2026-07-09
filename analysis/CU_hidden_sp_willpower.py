"""CU. Reverse-engineer hidden SP from cycle-level willpower mass-balance.

Model:

    dW/dt = -α · max(0, SP − FM)   (drain when below SP)
            + β · max(0, surplus)   (recharge when overeating)

Mass balance over a closed cycle (W returns to starting state):

    α · ∫_cycle max(0, SP − FM) dt  =  β · ∫_cycle max(0, intake − TDEE) dt

=>  drain_integral(SP) / recharge_integral  =  β/α  ≡ a single constant κ

If the model is correct with a fixed SP, then the ratio
R_i(SP) = drain_integral_i(SP) / recharge_integral_i  should be CONSTANT
across cycles (equal to β/α).

We sweep SP and find the value that minimizes CV(R) across cycles. That
is the SP at which the cycles are MOST mass-consistent.

Pre-registered convincing thresholds:
  - sharp interior minimum in CV(SP)
  - best CV < 0.5
  - W(t) integrated under the fitted (SP, κ) is bounded oscillation
  - per-cycle free SP_i clusters near the global SP (low SD, < 10 lb)
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent.parent
INTAKE_MULTIPLIER = 1.09


def cycle_drain(fm_arr, sp):
    """∫ max(0, SP − FM) dt summed over the cycle (lb·days)."""
    return float(np.sum(np.maximum(0.0, sp - fm_arr)))


def cycle_recharge(surplus_arr):
    """∫ max(0, surplus) dt summed over the cycle (cal)."""
    return float(np.sum(np.maximum(0.0, surplus_arr)))


def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    intake["calories"] = intake["calories"] * INTAKE_MULTIPLIER
    kf = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    cycles = pd.read_csv(ROOT / "analysis" / "CT_cycles.csv",
                          parse_dates=["peak_a_date", "trough_date", "peak_b_date"])

    df = intake[["date", "calories"]].merge(
        kf[["date", "fat_mass_lbs", "tdee"]], on="date")
    df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df = df.sort_values("date").reset_index(drop=True)
    df["fm"] = df["fat_mass_lbs"]
    df["surplus"] = df["calories"] - df["tdee"]
    df = df.dropna(subset=["fm", "tdee", "calories"]).reset_index(drop=True)
    pre = df[df["effective_level"] == 0].copy().reset_index(drop=True)

    # Build per-cycle daily slices once for fast SP sweeps
    cycle_data = []
    for _, c in cycles.iterrows():
        m = (pre["date"] >= c["peak_a_date"]) & (pre["date"] <= c["peak_b_date"])
        g = pre[m].copy()
        if len(g) < 30:
            continue
        cycle_data.append({
            "peak_a_date": c["peak_a_date"], "peak_b_date": c["peak_b_date"],
            "trough_date": c["trough_date"],
            "fm": g["fm"].values,
            "surplus": g["surplus"].values,
            "date": g["date"].values,
            "n": len(g),
            "fm_max": float(g["fm"].max()),
            "recharge": cycle_recharge(g["surplus"].values),
        })

    print("=" * 100)
    print(f"CU. Hidden SP from cycle-level mass balance.  n cycles = {len(cycle_data)}")
    print("=" * 100)

    # ── Sweep SP ──
    SP_grid = np.arange(50, 251, 1.0)
    rows = []
    for sp in SP_grid:
        ratios = []
        for c in cycle_data:
            d = cycle_drain(c["fm"], sp)
            r = c["recharge"]
            if d > 0 and r > 0:
                ratios.append(d / r)
        if len(ratios) < 5:
            continue
        ra = np.array(ratios)
        rows.append({
            "SP": sp,
            "n": len(ra),
            "mean_R": float(ra.mean()),
            "std_R": float(ra.std()),
            "cv_R": float(ra.std() / ra.mean()) if ra.mean() > 0 else float("inf"),
            "log_mean_R": float(np.mean(np.log(ra))),
            "log_std_R": float(np.std(np.log(ra))),
        })
    sweep = pd.DataFrame(rows)

    # ── Best SP by CV ──
    best_cv_idx = sweep["cv_R"].idxmin()
    best_sp = sweep.loc[best_cv_idx, "SP"]
    best_cv = sweep.loc[best_cv_idx, "cv_R"]
    best_logstd_idx = sweep["log_std_R"].idxmin()
    best_sp_log = sweep.loc[best_logstd_idx, "SP"]
    best_logstd = sweep.loc[best_logstd_idx, "log_std_R"]

    print(f"\nBest SP by CV(R):     {best_sp:.0f} lb  (CV = {best_cv:.3f})")
    print(f"Best SP by SD(log R): {best_sp_log:.0f} lb  (SD log R = {best_logstd:.3f})")

    print(f"\nCV(R) at SP=80:   {sweep[sweep['SP']==80]['cv_R'].values[0]:.3f}")
    print(f"CV(R) at SP=100:  {sweep[sweep['SP']==100]['cv_R'].values[0]:.3f}")
    print(f"CV(R) at SP=150:  {sweep[sweep['SP']==150]['cv_R'].values[0]:.3f}")
    print(f"CV(R) at SP=200:  {sweep[sweep['SP']==200]['cv_R'].values[0]:.3f}")
    print(f"CV(R) at SP=250:  {sweep[sweep['SP']==250]['cv_R'].values[0]:.3f}")

    print("\nPreregistered checks:")
    print(f"  Sharp interior minimum?  Best SP = {best_sp:.0f}, "
          f"between {SP_grid[0]:.0f} and {SP_grid[-1]:.0f} bounds: "
          f"{'YES' if (SP_grid[0] + 5 < best_sp < SP_grid[-1] - 5) else 'NO'}")
    print(f"  Best CV < 0.5?            {best_cv:.3f} → {'YES' if best_cv < 0.5 else 'NO'}")

    # ── Per-cycle free SP_i: solve cycle-by-cycle for SP_i such that R_i = κ_global ──
    # κ_global = mean R at best_sp
    print("\n" + "=" * 100)
    print("PER-CYCLE free SP_i — does it cluster?")
    print("=" * 100)
    # Use mean_R at best_sp as the shared κ
    kappa = float(sweep.loc[best_cv_idx, "mean_R"])
    print(f"  Using κ = mean_R(SP={best_sp:.0f}) = {kappa:.4f}\n")

    print(f"  {'cycle (peak_a)':<14} {'SP_i':>6} {'drain':>8} {'recharge':>9}")
    sp_per_cycle = []
    for c in cycle_data:
        # Find SP_i such that cycle_drain(fm, SP_i) / recharge = κ
        target_drain = kappa * c["recharge"]

        def f(sp):
            return cycle_drain(c["fm"], sp) - target_drain

        # cycle_drain is monotone increasing in SP
        sp_lo, sp_hi = c["fm_max"] - 5, 500.0
        # f at sp_lo: drain ≈ small; f at sp_hi: drain very large
        if f(sp_lo) < 0 < f(sp_hi):
            sp_i = brentq(f, sp_lo, sp_hi)
        else:
            sp_i = float("nan")
        sp_per_cycle.append({
            "peak_a_date": pd.Timestamp(c["peak_a_date"]),
            "trough_date": pd.Timestamp(c["trough_date"]),
            "peak_b_date": pd.Timestamp(c["peak_b_date"]),
            "sp_i": sp_i,
            "drain_i": cycle_drain(c["fm"], sp_i) if not np.isnan(sp_i) else np.nan,
            "recharge": c["recharge"],
            "fm_max_in_cycle": c["fm_max"],
        })
        print(f"  {pd.Timestamp(c['peak_a_date']).strftime('%Y-%m-%d'):<14} "
              f"{sp_i:>6.1f} {cycle_drain(c['fm'], sp_i):>8.0f} {c['recharge']:>9.0f}")

    sp_df = pd.DataFrame(sp_per_cycle)
    sp_df = sp_df.dropna(subset=["sp_i"])
    print(f"\n  SP_i: mean = {sp_df['sp_i'].mean():.1f} lb, "
          f"median = {sp_df['sp_i'].median():.1f}, "
          f"SD = {sp_df['sp_i'].std():.2f}, "
          f"min = {sp_df['sp_i'].min():.1f}, max = {sp_df['sp_i'].max():.1f}")
    print(f"  IQR: {sp_df['sp_i'].quantile(0.25):.1f} to {sp_df['sp_i'].quantile(0.75):.1f}")

    # Does SP_i drift with FM_max in cycle?
    if len(sp_df) > 5:
        r_drift = float(np.corrcoef(sp_df["fm_max_in_cycle"], sp_df["sp_i"])[0, 1])
        print(f"  Corr(SP_i, FM_max_in_cycle) = {r_drift:+.3f}  "
              f"({'tracks FM' if abs(r_drift) > 0.5 else 'independent of FM'})")

    # ── Compute W(t) under fitted (SP, κ) ──
    print("\n" + "=" * 100)
    print("W(t) trajectory under fitted (SP, κ)")
    print("=" * 100)
    # dW/dt = -α·max(0, SP−FM) + β·max(0, surplus)
    # Set α = 1, β = α·κ ⇒ dW/dt in lb-day units
    # W(t) cumulative integral. Initial W = 0; track shape.
    SP_used = best_sp
    pre["drain_term"] = np.maximum(0, SP_used - pre["fm"])
    # Convert: drain in lb·day units; recharge in cal, divide by 3500 to lb-equiv
    # Use β/α = κ from drain/recharge ratio directly
    # If drain has units lb·day and recharge has units cal, then
    # for them to balance we need β/α with units (lb·day)/cal
    # so that α·drain (lb·day) = β·recharge·(lb·day/cal)
    # Means dW/dt in units of lb·day/day = lb after time-integrating
    pre["dW"] = -pre["drain_term"] + (1.0 / kappa) * np.maximum(0, pre["surplus"])
    pre["W"] = pre["dW"].cumsum()

    # Detrend: any constant drift means α and β aren't perfectly balanced; that's fine
    pre["W_centered"] = pre["W"] - pre["W"].rolling(365, min_periods=30).mean()

    # Stats on W
    print(f"  W absolute range: {pre['W'].min():.1f} to {pre['W'].max():.1f}")
    print(f"  W centered (annual-detrended) range: {pre['W_centered'].min():.1f} to "
          f"{pre['W_centered'].max():.1f}")
    print(f"  W at troughs vs peaks (centered):")
    for _, c in cycles.iterrows():
        try:
            t_idx = pre[pre["date"] == c["trough_date"]].index[0]
            p_idx = pre[pre["date"] == c["peak_b_date"]].index[0]
            print(f"    cycle {c['peak_a_date'].strftime('%Y-%m-%d')}: "
                  f"W_trough = {pre.loc[t_idx, 'W_centered']:+.1f}, "
                  f"W_peak = {pre.loc[p_idx, 'W_centered']:+.1f}")
        except (IndexError, KeyError):
            continue

    # Check: SD of W_centered at troughs vs at peaks
    trough_dates = pd.to_datetime(cycles["trough_date"]).tolist()
    peak_dates = pd.to_datetime(cycles["peak_b_date"]).tolist()
    pre_idx = pre.set_index("date")
    w_at_troughs = []
    w_at_peaks = []
    for d in trough_dates:
        if d in pre_idx.index:
            v = pre_idx.loc[d, "W_centered"]
            if not np.isnan(v):
                w_at_troughs.append(float(v))
    for d in peak_dates:
        if d in pre_idx.index:
            v = pre_idx.loc[d, "W_centered"]
            if not np.isnan(v):
                w_at_peaks.append(float(v))

    if w_at_troughs and w_at_peaks:
        wt = np.array(w_at_troughs)
        wp = np.array(w_at_peaks)
        print(f"\n  W at troughs: mean = {wt.mean():+.2f}, SD = {wt.std():.2f}, n = {len(wt)}")
        print(f"  W at peaks:   mean = {wp.mean():+.2f}, SD = {wp.std():.2f}, n = {len(wp)}")
        print(f"  Mean(W_peak) − Mean(W_trough) = {wp.mean() - wt.mean():+.2f}")
        print(f"  → battery range estimate: {wp.mean() - wt.mean():+.2f}")

    # ── Plots ──
    fig, axes = plt.subplots(3, 1, figsize=(14, 11))

    ax = axes[0]
    ax.plot(sweep["SP"], sweep["cv_R"], color="tab:blue", lw=1.5)
    ax.axvline(best_sp, color="tab:red", ls="--", label=f"min CV at SP={best_sp:.0f}")
    ax.axhline(0.5, color="grey", ls=":", label="prereg threshold CV=0.5")
    ax.set_xlabel("Candidate SP (lb)")
    ax.set_ylabel("CV of R = drain/recharge across cycles")
    ax.set_title(f"Mass-balance scan. Best CV = {best_cv:.3f} at SP = {best_sp:.0f} lb")
    ax.legend()
    ax.grid(alpha=0.3)

    ax = axes[1]
    # Plot SP_i over time
    if len(sp_df):
        ax.plot(sp_df["peak_a_date"], sp_df["sp_i"], "o-", color="tab:purple", label="SP_i (per cycle)")
        ax.axhline(best_sp, color="tab:red", ls="--", label=f"global SP={best_sp:.0f}")
        ax.plot(sp_df["peak_a_date"], sp_df["fm_max_in_cycle"], "s-", color="grey",
                alpha=0.5, label="FM_max in cycle")
        ax.set_ylabel("FM (lb)")
        ax.set_xlabel("Cycle date")
        ax.set_title(f"Per-cycle SP_i: SD = {sp_df['sp_i'].std():.2f} lb, "
                     f"corr with FM_max = {np.corrcoef(sp_df['fm_max_in_cycle'], sp_df['sp_i'])[0,1]:+.2f}")
        ax.legend()
        ax.grid(alpha=0.3)

    ax = axes[2]
    ax.plot(pre["date"], pre["W_centered"], color="tab:orange", lw=1.0)
    for _, c in cycles.iterrows():
        ax.axvline(pd.Timestamp(c["trough_date"]), color="tab:blue", alpha=0.2)
        ax.axvline(pd.Timestamp(c["peak_b_date"]), color="tab:red", alpha=0.2)
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlabel("Date")
    ax.set_ylabel("W (lb-equivalent, annual-detrended)")
    ax.set_title("W(t) integrated under fitted (SP, κ) — battery oscillation?")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    plot_path = ROOT / "analysis" / "CU_hidden_sp_willpower.png"
    fig.savefig(plot_path, dpi=140)
    print(f"\nSaved {plot_path}")

    sweep.to_csv(ROOT / "analysis" / "CU_sp_sweep.csv", index=False)
    sp_df.to_csv(ROOT / "analysis" / "CU_sp_per_cycle.csv", index=False)

    # ── Final verdict ──
    print("\n" + "=" * 100)
    print("VERDICT vs PREREGISTRATION")
    print("=" * 100)
    sharp_min = (SP_grid[0] + 5 < best_sp < SP_grid[-1] - 5)
    cv_low = best_cv < 0.5
    sp_clustered = sp_df["sp_i"].std() < 10 if len(sp_df) > 5 else False
    print(f"  [{'PASS' if sharp_min else 'FAIL'}] Sharp interior minimum: "
          f"best SP = {best_sp:.0f} (bounds {SP_grid[0]:.0f}-{SP_grid[-1]:.0f})")
    print(f"  [{'PASS' if cv_low else 'FAIL'}] Best CV < 0.5: {best_cv:.3f}")
    print(f"  [{'PASS' if sp_clustered else 'FAIL'}] Per-cycle SP_i SD < 10 lb: "
          f"{sp_df['sp_i'].std():.2f} lb")


if __name__ == "__main__":
    main()
