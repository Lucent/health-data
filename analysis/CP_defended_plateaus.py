"""CP. Are FM plateaus defended?

Test the defended-plateau hypothesis directly:

1. Detect FM trajectory change points (plateaus and transitions) without
   manual labels — use piecewise-constant DP fit, K chosen by BIC elbow.
2. For each plateau (stable segment, low slope, ≥6 months):
     a. Compute defended-band statistics: plateau mean, IQR, SD around mean.
     b. Find perturbation runs INSIDE the plateau: 7-day windows where
        mean surplus magnitude > 500 cal/day. Track FM excursion and
        FM 60 days after perturbation ends. Compute snap-back fraction.
     c. Compute lag-30 autocorrelation of within-segment FM residual.
        Defended → low ACF (fast snap-back). Passive → high ACF (drift).
3. For each transition: characterize what made it different from
   within-plateau perturbations the body successfully absorbed.

Defended-plateau predicts:
  - Within-plateau perturbations snap back: median snap-back > 0.5
  - FM lag-30 ACF inside plateaus is low (< 0.5)
  - Transition-triggering perturbations are LARGER or LONGER than
    within-plateau perturbations the body absorbed
  - Plateau mean shifts ARE the set point movements; segment-by-segment.

Passive-resting predicts:
  - Within-plateau perturbations don't snap back; they leave a permanent
    offset proportional to cumulative cal surplus / 3500
  - Lag-30 ACF inside plateaus is high (≥ 0.8)
  - No threshold separates plateau-internal from transition-triggering
    perturbations — it's just that some happen to be larger
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
INTAKE_MULTIPLIER = 1.09


def piecewise_constant_dp(signal, max_k=12, min_size=60):
    """DP for optimal K-segmentation minimizing within-segment SS.
    Returns (cost_by_k, reconstruct_fn) where reconstruct_fn(K) gives breakpoints."""
    n = len(signal)
    cs = np.concatenate([[0.0], np.cumsum(signal.astype(float))])
    cs2 = np.concatenate([[0.0], np.cumsum(signal.astype(float) ** 2)])

    def seg_cost(i, j):
        if j - i <= 0:
            return 0.0
        m = (cs[j] - cs[i]) / (j - i)
        return float(cs2[j] - cs2[i] - m * (cs[j] - cs[i]))

    INF = 1e18
    f = np.full((max_k + 1, n + 1), INF)
    f[0, 0] = 0.0
    bp = np.full((max_k + 1, n + 1), -1, dtype=int)

    for k in range(1, max_k + 1):
        lo_i = k * min_size
        for i in range(lo_i, n + 1):
            best = INF
            best_j = -1
            for j in range((k - 1) * min_size, i - min_size + 1):
                if f[k - 1, j] >= INF:
                    continue
                c = f[k - 1, j] + seg_cost(j, i)
                if c < best:
                    best = c
                    best_j = j
            f[k, i] = best
            bp[k, i] = best_j

    def reconstruct(k):
        if f[k, n] >= INF:
            return None
        breaks = [n]
        cur = n
        kk = k
        while kk > 0:
            cur = bp[kk, cur]
            breaks.append(cur)
            kk -= 1
        breaks.reverse()
        return breaks

    return f[:, n], reconstruct


def main():
    print("=" * 80)
    print("CP. Defended-plateau test")
    print("=" * 80)

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
    pre = df[df["effective_level"] == 0].reset_index(drop=True)

    # Drop rows with NaN FM (typically pre-first-scan)
    pre = pre.dropna(subset=["fm", "tdee", "calories"]).reset_index(drop=True)
    print(f"\nValid pre-tirz days: {len(pre)}")

    # ── Subsample to weekly resolution for DP speed ──
    pre["week"] = (pre.index // 7)
    weekly = pre.groupby("week").agg(
        date=("date", "first"),
        fm=("fm", "mean"),
        surplus=("surplus", "mean"),
        calories=("calories", "mean"),
        tdee=("tdee", "mean"),
    ).reset_index(drop=True)
    weekly = weekly.dropna(subset=["fm"]).reset_index(drop=True)
    n_weeks = len(weekly)
    print(f"Weekly resolution: {n_weeks} weeks ({n_weeks * 7 / 365:.1f} years)")

    # ── DP segmentation ──
    print("\nRunning piecewise-constant DP segmentation...")
    fm_w = weekly["fm"].values
    cost_by_k, reconstruct = piecewise_constant_dp(
        fm_w, max_k=10, min_size=12)  # min 12 weeks ≈ 90 days

    # BIC-style elbow: pick K where marginal improvement < threshold
    print("\nDP cost by K:")
    for k in range(1, 11):
        marginal = (cost_by_k[k - 1] - cost_by_k[k]) if k > 0 else 0
        print(f"  K={k:>2}: cost={cost_by_k[k]:>10.1f}  Δ={marginal:>+10.1f}")

    K = 7
    print(f"\nUsing K={K} segments")
    breaks = reconstruct(K)

    # Build segment list
    segments = []
    for i in range(len(breaks) - 1):
        s, e = breaks[i], breaks[i + 1]
        seg_w = weekly.iloc[s:e]
        seg_d = pre[(pre["date"] >= seg_w["date"].iloc[0]) &
                    (pre["date"] <= seg_w["date"].iloc[-1])]
        if len(seg_d) < 60:
            continue
        fm_mean = float(seg_d["fm"].mean())
        fm_sd = float(seg_d["fm"].std())
        fm_iqr = float(seg_d["fm"].quantile(0.75) - seg_d["fm"].quantile(0.25))
        slope_lb_per_day = (seg_d["fm"].iloc[-1] - seg_d["fm"].iloc[0]) / len(seg_d)
        # Lag-30 autocorrelation of FM residual (FM minus segment-local trend)
        x = np.arange(len(seg_d))
        coefs = np.polyfit(x, seg_d["fm"].values, 1)
        resid = seg_d["fm"].values - (coefs[0] * x + coefs[1])
        if len(resid) > 60:
            r_acf = float(np.corrcoef(resid[:-30], resid[30:])[0, 1])
        else:
            r_acf = float("nan")
        segments.append({
            "idx": i,
            "start": seg_d["date"].iloc[0],
            "end": seg_d["date"].iloc[-1],
            "n_days": len(seg_d),
            "fm_mean": fm_mean,
            "fm_sd": fm_sd,
            "fm_iqr": fm_iqr,
            "slope": slope_lb_per_day,
            "intake": float(seg_d["calories"].mean()),
            "tdee": float(seg_d["tdee"].mean()),
            "lag30_acf": r_acf,
            "is_plateau": abs(slope_lb_per_day) < 0.02 and len(seg_d) >= 180,
        })

    seg_df = pd.DataFrame(segments)
    print("\n" + "=" * 110)
    print("SEGMENTS")
    print("=" * 110)
    print(f"  {'idx':>3} {'start':<11} {'end':<11} {'days':>5} "
          f"{'FM mean':>8} {'FM SD':>6} {'IQR':>5} {'slope':>9} "
          f"{'intake':>7} {'TDEE':>6} {'lag30 ACF':>9} {'plateau?':>9}")
    for _, r in seg_df.iterrows():
        print(f"  {int(r['idx']):>3} {r['start'].strftime('%Y-%m-%d'):<11} "
              f"{r['end'].strftime('%Y-%m-%d'):<11} {int(r['n_days']):>5} "
              f"{r['fm_mean']:>8.1f} {r['fm_sd']:>6.2f} {r['fm_iqr']:>5.2f} "
              f"{r['slope']:>+9.4f} {r['intake']:>7.0f} {r['tdee']:>6.0f} "
              f"{r['lag30_acf']:>+9.3f} {'YES' if r['is_plateau'] else 'no':>9}")

    # ── Within-plateau perturbation analysis ──
    print("\n" + "=" * 110)
    print("WITHIN-PLATEAU PERTURBATIONS — does FM snap back?")
    print("=" * 110)
    pre["surplus_7d"] = pre["surplus"].rolling(7, min_periods=7).mean()

    all_perturbations = []
    for _, s in seg_df.iterrows():
        if not s["is_plateau"]:
            continue
        seg = pre[(pre["date"] >= s["start"]) & (pre["date"] <= s["end"])].copy()
        plateau_mean = s["fm_mean"]
        # Find perturbation runs: contiguous days where |surplus_7d| > 500
        sp = seg["surplus_7d"].values
        in_run = False
        run_start = None
        runs = []
        for i in range(len(sp)):
            if not np.isnan(sp[i]) and abs(sp[i]) > 500:
                if not in_run:
                    in_run = True
                    run_start = i
            else:
                if in_run:
                    if i - run_start >= 5:  # at least 5 days
                        runs.append((run_start, i))
                    in_run = False
        if in_run and len(sp) - run_start >= 5:
            runs.append((run_start, len(sp)))

        for r_start, r_end in runs:
            t_end = seg["date"].iloc[r_end - 1]
            fm_at_end = seg["fm"].iloc[r_end - 1]
            t_60 = t_end + pd.Timedelta(days=60)
            future_mask = (pre["date"] >= t_60 - pd.Timedelta(days=3)) & \
                          (pre["date"] <= t_60 + pd.Timedelta(days=3))
            if future_mask.sum() == 0:
                continue
            fm_60d = float(pre.loc[future_mask, "fm"].mean())
            mean_surplus = float(np.mean(sp[r_start:r_end]))
            cum_surplus = float(np.sum(seg["surplus"].iloc[r_start:r_end]))
            excursion = fm_at_end - plateau_mean
            residual_60d = fm_60d - plateau_mean
            if abs(excursion) < 0.5:
                continue  # not a real perturbation
            snap_back = 1.0 - residual_60d / excursion
            all_perturbations.append({
                "segment": int(s["idx"]),
                "start": seg["date"].iloc[r_start],
                "end": t_end,
                "duration_d": r_end - r_start,
                "mean_surplus_cal_d": mean_surplus,
                "cum_surplus_cal": cum_surplus,
                "fm_at_end": fm_at_end,
                "fm_60d_after": fm_60d,
                "plateau_mean": plateau_mean,
                "excursion_lb": excursion,
                "residual_60d_lb": residual_60d,
                "snap_back": snap_back,
            })

    perts = pd.DataFrame(all_perturbations)
    if not perts.empty:
        print(f"\nFound {len(perts)} within-plateau perturbations (≥5 days, |surplus| > 500 cal/d)")
        print(f"\n  {'seg':>3} {'start':<11} {'days':>4} "
              f"{'mean_surp':>9} {'cum_surp':>9} {'excursion':>9} "
              f"{'resid_60d':>9} {'snap_back':>9}")
        for _, p in perts.iterrows():
            print(f"  {int(p['segment']):>3} {p['start'].strftime('%Y-%m-%d'):<11} "
                  f"{int(p['duration_d']):>4} "
                  f"{p['mean_surplus_cal_d']:>+9.0f} {p['cum_surplus_cal']:>+9.0f} "
                  f"{p['excursion_lb']:>+9.2f} {p['residual_60d_lb']:>+9.2f} "
                  f"{p['snap_back']:>+9.2f}")

        print(f"\nMedian snap-back across all within-plateau perturbations: "
              f"{perts['snap_back'].median():+.3f}")
        print(f"Fraction with snap_back > 0.5 (substantial return): "
              f"{(perts['snap_back'] > 0.5).mean():.2%}")
        print(f"Fraction with snap_back > 0.8 (near-full return): "
              f"{(perts['snap_back'] > 0.8).mean():.2%}")

        # Stratify by direction
        up = perts[perts["excursion_lb"] > 0.5]
        down = perts[perts["excursion_lb"] < -0.5]
        if len(up) > 0:
            print(f"\nUpward perturbations (FM rose):     n={len(up):>2}, "
                  f"median snap-back = {up['snap_back'].median():+.3f}")
        if len(down) > 0:
            print(f"Downward perturbations (FM fell):    n={len(down):>2}, "
                  f"median snap-back = {down['snap_back'].median():+.3f}")

    # ── Plot ──
    fig, axes = plt.subplots(2, 1, figsize=(14, 9), sharex=True)

    ax = axes[0]
    ax.plot(pre["date"], pre["fm"], color="black", lw=0.8, alpha=0.8)
    plateau_color = "#aed6f1"
    transition_color = "#f5b7b1"
    for _, s in seg_df.iterrows():
        c = plateau_color if s["is_plateau"] else transition_color
        ax.axvspan(s["start"], s["end"], alpha=0.4, color=c)
        # Plateau-mean horizontal line
        ax.hlines(s["fm_mean"], s["start"], s["end"], color="darkblue",
                  lw=1.5, alpha=0.7)
        # Annotate
        ax.text(s["start"] + (s["end"] - s["start"]) / 2,
                s["fm_mean"] + 2,
                f"P{int(s['idx'])}\n{s['fm_mean']:.0f}\nACF{s['lag30_acf']:.2f}",
                ha="center", fontsize=7, alpha=0.8)
    ax.set_ylabel("FM (lb)")
    ax.set_title("FM with auto-detected segments. Plateau = blue, transition = red")
    ax.grid(alpha=0.3)

    ax = axes[1]
    if not perts.empty:
        # Plot snap-back vs excursion
        sc = ax.scatter(perts["excursion_lb"], perts["snap_back"],
                        c=perts["duration_d"], cmap="viridis", s=50,
                        edgecolor="black")
        plt.colorbar(sc, ax=ax, label="duration (days)")
        ax.axhline(1.0, color="green", ls="--", alpha=0.5, label="full snap-back")
        ax.axhline(0.0, color="red", ls="--", alpha=0.5, label="no snap-back")
        ax.axvline(0, color="black", lw=0.5)
        ax.set_xlabel("FM excursion at perturbation end (lb)")
        ax.set_ylabel("Snap-back fraction (1 = full return to plateau mean)")
        ax.set_title(f"Within-plateau perturbation snap-back  "
                     f"(median = {perts['snap_back'].median():+.2f})")
        ax.legend()
        ax.grid(alpha=0.3)

    fig.tight_layout()
    plot_path = ROOT / "analysis" / "CP_defended_plateaus.png"
    fig.savefig(plot_path, dpi=140)
    print(f"\nSaved {plot_path}")

    # Save data
    seg_df.to_csv(ROOT / "analysis" / "CP_segments.csv", index=False)
    if not perts.empty:
        perts.to_csv(ROOT / "analysis" / "CP_perturbations.csv", index=False)


if __name__ == "__main__":
    main()
