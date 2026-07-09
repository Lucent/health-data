"""CW. Find a willpower variable that exhausts at the right time.

For each candidate drain rate function f(t), compute cumulative drain along
each FALL (peak_a → trough). If the model is correct, cumulative drain at
trough should be CONSTANT across all 19 falls (= total willpower budget W_max).

The candidate drain rate functions:

  Model A — deficit only:
      drain rate = max(0, TDEE − intake)   [cal/day, only on deficit days]

  Model B — pure FM exposure:
      B1: drain rate = 1                    [time-only; CV = CV of n_days]
      B2: drain rate = FM                   [FM × time]
      B3: drain rate = FM × deficit         [combined]
      B4: drain rate = deficit / FM         [deficit harder at low FM]

  Model C — gap-driven, SP free per cycle:
      drain rate = max(0, SP_i − BW)        [BW units]
      For each fall: SP_i × n − Σ BW = W_max
      → SP_i = (W_max + Σ BW) / n
      Sweep W_max, find value where SD of SP_i is minimum.

  Model D — gap × deficit:
      drain rate = max(0, SP_i − BW) × max(0, TDEE − intake)
      Two-parameter sweep (W_max, α) — but here just sweep W_max with
      gap-weighted deficit integral.

The model with lowest CV of total-drain-at-trough has the most consistent
willpower budget. Model C/D additionally produce per-cycle SP_i values
which we examine for clustering.
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
    df["deficit_pos"] = np.maximum(0.0, -df["surplus"])
    df = df.dropna(subset=["fm", "lm", "tdee", "calories"]).reset_index(drop=True)
    pre = df[df["effective_level"] == 0].copy().reset_index(drop=True)

    # Extract FALL portion of each cycle (peak_a → trough)
    falls = []
    for _, c in cycles.iterrows():
        m = (pre["date"] >= c["peak_a_date"]) & (pre["date"] <= c["trough_date"])
        g = pre[m].copy().reset_index(drop=True)
        if len(g) < 14:
            continue
        falls.append({
            "peak_a_date": c["peak_a_date"], "trough_date": c["trough_date"],
            "fm": g["fm"].values, "bw": g["bw"].values,
            "deficit_pos": g["deficit_pos"].values,
            "surplus": g["surplus"].values,
            "n_days": len(g),
            "fall_amplitude_lb": float(g["bw"].iloc[0] - g["bw"].iloc[-1]),
        })
    print("=" * 100)
    print(f"CW. Willpower exhaustion test, {len(falls)} falls")
    print("=" * 100)

    def report_cv(label, drains, units=""):
        d = np.array(drains)
        mean, sd, cv = d.mean(), d.std(), d.std() / abs(d.mean()) if d.mean() != 0 else float("inf")
        med = np.median(d)
        iqr = np.percentile(d, 75) - np.percentile(d, 25)
        print(f"  {label:<48} mean = {mean:>10.0f} {units} "
              f"SD = {sd:>9.0f} CV = {cv:.3f}  IQR/median = {iqr/abs(med):.3f}")
        return mean, sd, cv

    cv_summary = {}

    # ── Model A — deficit only ──
    print("\nMODEL A — drain ∝ deficit, summed over deficit days:")
    drains_A = [float(f["deficit_pos"].sum()) for f in falls]
    _, _, cv_summary["A: deficit"] = report_cv("Σ max(0, TDEE − intake)", drains_A, "cal")

    # ── Model B variants ──
    print("\nMODEL B — drain ∝ functions of FM (no SP, no deficit gating):")
    drains_B1 = [float(f["n_days"]) for f in falls]
    _, _, cv_summary["B1: time"] = report_cv("Σ 1 (= fall duration)", drains_B1, "days")

    drains_B2 = [float(f["fm"].sum()) for f in falls]
    _, _, cv_summary["B2: FM × time"] = report_cv("Σ FM (= mean FM × n_days)", drains_B2, "lb·d")

    drains_B3 = [float((f["fm"] * f["deficit_pos"]).sum()) for f in falls]
    _, _, cv_summary["B3: FM × deficit"] = report_cv("Σ FM × deficit", drains_B3, "lb·cal")

    drains_B4 = [float((f["deficit_pos"] / np.maximum(f["fm"], 5)).sum()) for f in falls]
    _, _, cv_summary["B4: deficit / FM"] = report_cv("Σ deficit / FM", drains_B4, "cal/lb")

    drains_B5 = [float((f["bw"] * f["deficit_pos"]).sum()) for f in falls]
    _, _, cv_summary["B5: BW × deficit"] = report_cv("Σ BW × deficit", drains_B5, "lb·cal")

    # ── Model C — gap-driven (SP free per cycle) ──
    print("\nMODEL C — drain ∝ (SP_i − BW) with SP_i free per cycle:")
    print("  Sweeping W_max (in lb·days). For each, SP_i = (W_max + Σ BW) / n_days.")
    print("  Find W_max minimizing SD of SP_i across cycles.")
    W_grid = np.logspace(2.5, 6, 200)  # 316 to 1M lb·days
    sp_records = []
    for W in W_grid:
        sps = np.array([(W + f["bw"].sum()) / f["n_days"] for f in falls])
        sp_records.append({
            "W_max": W,
            "sp_mean": sps.mean(),
            "sp_sd": sps.std(),
            "sp_cv": sps.std() / sps.mean(),
        })
    sp_sweep = pd.DataFrame(sp_records)
    best_idx_sd = sp_sweep["sp_sd"].idxmin()
    best_W = sp_sweep.loc[best_idx_sd, "W_max"]
    best_SD = sp_sweep.loc[best_idx_sd, "sp_sd"]
    best_mean = sp_sweep.loc[best_idx_sd, "sp_mean"]
    print(f"\n  Best W_max (min SD of SP_i): {best_W:.0f} lb·days")
    print(f"  At best W_max:")
    print(f"    SP_i mean = {best_mean:.1f} lb (BW)")
    print(f"    SP_i SD   = {best_SD:.2f} lb")
    print(f"    SP_i CV   = {best_SD/best_mean:.3f}")

    sps_at_best = np.array([(best_W + f["bw"].sum()) / f["n_days"] for f in falls])
    print(f"\n  Per-fall SP_i at best W_max:")
    for f, sp_i in zip(falls, sps_at_best):
        print(f"    {f['peak_a_date'].strftime('%Y-%m-%d')} → "
              f"{f['trough_date'].strftime('%Y-%m-%d')}  "
              f"({f['n_days']:>3}d, ΔBW {f['fall_amplitude_lb']:>+5.1f}): "
              f"SP_i = {sp_i:>6.1f} lb")

    # Try W_max anchored to give SP_i ≈ 290 (subject's lifetime peak hypothesis)
    target_sp = 290.0
    # mean SP = W_mean(1/n) + mean(mean_bw_i). Solve for W:
    one_over_n = np.array([1.0 / f["n_days"] for f in falls])
    mean_bw_i = np.array([f["bw"].mean() for f in falls])
    W_for_290 = (target_sp - mean_bw_i.mean()) / one_over_n.mean()
    if W_for_290 > 0:
        sps_at_290 = np.array([(W_for_290 + f["bw"].sum()) / f["n_days"] for f in falls])
        print(f"\n  At W_max = {W_for_290:.0f} (chosen so mean SP_i = {target_sp:.0f}):")
        print(f"    SP_i mean = {sps_at_290.mean():.1f}, SD = {sps_at_290.std():.2f}")
        print(f"    Per-fall SP_i: {sorted([round(x, 1) for x in sps_at_290])}")

    # ── Model D — gap × deficit ──
    print("\nMODEL D — drain ∝ (SP_i − BW) × deficit_pos with SP_i free:")
    # For each fall: ∫(SP_i - BW) · deficit_pos dt = W_max
    # SP_i · Σ(deficit_pos) - Σ(BW · deficit_pos) = W_max
    # SP_i = (W_max + Σ(BW · deficit_pos)) / Σ(deficit_pos)
    sp_records_D = []
    W_grid_D = np.logspace(3, 8, 200)  # cal-lb units, larger
    for W in W_grid_D:
        sps = []
        for f in falls:
            denom = float(f["deficit_pos"].sum())
            num = float((f["bw"] * f["deficit_pos"]).sum())
            if denom < 1e-3:
                continue
            sps.append((W + num) / denom)
        if len(sps) < 5:
            continue
        sps = np.array(sps)
        sp_records_D.append({
            "W_max": W,
            "sp_mean": sps.mean(),
            "sp_sd": sps.std(),
        })
    sp_sweep_D = pd.DataFrame(sp_records_D)
    if len(sp_sweep_D):
        best_idx_D = sp_sweep_D["sp_sd"].idxmin()
        bestD_W = sp_sweep_D.loc[best_idx_D, "W_max"]
        bestD_SD = sp_sweep_D.loc[best_idx_D, "sp_sd"]
        bestD_mean = sp_sweep_D.loc[best_idx_D, "sp_mean"]
        print(f"  Best W_max: {bestD_W:.0f} (lb·cal)")
        print(f"  SP_i mean = {bestD_mean:.1f}, SD = {bestD_SD:.2f}, CV = {bestD_SD/bestD_mean:.3f}")
        sps_at_best_D = []
        for f in falls:
            denom = float(f["deficit_pos"].sum())
            num = float((f["bw"] * f["deficit_pos"]).sum())
            if denom < 1e-3:
                sps_at_best_D.append(np.nan)
                continue
            sps_at_best_D.append((bestD_W + num) / denom)
        print(f"\n  Per-fall SP_i at best W_max (gap × deficit):")
        for f, sp_i in zip(falls, sps_at_best_D):
            print(f"    {f['peak_a_date'].strftime('%Y-%m-%d')}: SP_i = {sp_i:>7.1f} lb")

    # ── Summary table ──
    print("\n" + "=" * 100)
    print("CV SUMMARY (lower = more consistent willpower budget across cycles)")
    print("=" * 100)
    cv_summary["C: SP free, no deficit gating"] = best_SD / best_mean
    if len(sp_sweep_D):
        cv_summary["D: SP free, deficit-gated"] = bestD_SD / bestD_mean
    for k, v in sorted(cv_summary.items(), key=lambda x: x[1]):
        marker = "  ★" if v == min(cv_summary.values()) else "   "
        print(f"  {marker} {k:<40} CV = {v:.3f}")

    # ── Plots ──
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel 1: cumulative drain at trough across falls, normalized for each model
    ax = axes[0, 0]
    models_to_plot = [
        ("A: deficit", np.array(drains_A) / np.mean(drains_A)),
        ("B1: time", np.array(drains_B1) / np.mean(drains_B1)),
        ("B2: FM × time", np.array(drains_B2) / np.mean(drains_B2)),
        ("B3: FM × deficit", np.array(drains_B3) / np.mean(drains_B3)),
    ]
    for label, vals in models_to_plot:
        ax.plot(range(len(vals)), vals, "o-", label=label, alpha=0.7)
    ax.axhline(1.0, color="black", lw=0.5)
    ax.set_xlabel("fall index")
    ax.set_ylabel("drain at trough (normalized to mean)")
    ax.set_title("Cumulative drain at trough — does each model give a flat line?")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # Panel 2: SD of SP_i vs W_max for Model C
    ax = axes[0, 1]
    ax.plot(sp_sweep["W_max"], sp_sweep["sp_sd"], color="tab:purple")
    ax.axvline(best_W, color="tab:red", ls="--", label=f"min SD at W={best_W:.0f}")
    ax.set_xscale("log")
    ax.set_xlabel("W_max (lb·days)")
    ax.set_ylabel("SD of SP_i across cycles (lb)")
    ax.set_title(f"Model C: SD of per-cycle SP_i vs W_max  (min SD = {best_SD:.1f} lb)")
    ax.legend()
    ax.grid(alpha=0.3)

    # Panel 3: SP_i values per cycle at best W_max (Model C)
    ax = axes[1, 0]
    fall_dates = [f["peak_a_date"] for f in falls]
    ax.plot(fall_dates, sps_at_best, "o-", color="tab:purple", label=f"SP_i at W_max={best_W:.0f}")
    ax.axhline(sps_at_best.mean(), color="tab:red", ls="--",
               label=f"mean = {sps_at_best.mean():.0f}")
    ax.axhspan(285, 300, color="orange", alpha=0.2, label="hypothesized 285-300")
    ax.set_xlabel("Fall start date")
    ax.set_ylabel("SP_i (BW, lb)")
    ax.set_title("Model C: per-fall SP_i at best W_max")
    ax.legend()
    ax.grid(alpha=0.3)

    # Panel 4: CV bar chart
    ax = axes[1, 1]
    items = sorted(cv_summary.items(), key=lambda x: x[1])
    labels = [k for k, _ in items]
    vals = [v for _, v in items]
    colors = ["tab:green" if v == min(vals) else "tab:blue" for v in vals]
    ax.barh(range(len(labels)), vals, color=colors)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0.3, color="grey", ls=":", label="rough good-fit threshold")
    ax.set_xlabel("CV of cumulative drain at trough (or SP_i SD/mean for C/D)")
    ax.set_title("Model comparison — lower CV = constant willpower budget")
    ax.legend()
    ax.grid(alpha=0.3, axis="x")

    fig.tight_layout()
    plot_path = ROOT / "analysis" / "CW_willpower_exhaustion.png"
    fig.savefig(plot_path, dpi=140)
    print(f"\nSaved {plot_path}")


if __name__ == "__main__":
    main()
