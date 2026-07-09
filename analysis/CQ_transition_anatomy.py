"""CQ. Anatomy of the disjoint plateau-to-plateau jumps.

CP found discrete plateaus separated by transition periods. This script looks
at the transitions themselves and asks two questions:

  (1) Absorption vs breakthrough. For each plateau, what was the LARGEST
      cumulative surplus (over 7/14/30/60/90-day windows) the plateau
      absorbed without breaking? For each transition, what was the cumulative
      surplus that drove it? Is there a clean threshold?

  (2) What does the jump look like inside? Is FM rising smoothly through the
      transition, or in discrete steps? Where does the cumulative surplus
      come from — a few big binge episodes, or sustained low-grade surplus?
      Are there sub-plateaus within transitions?

No model. Just numbers and shape.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
INTAKE_MULTIPLIER = 1.09
BINGE_THRESHOLD = 1000


def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    intake["calories"] = intake["calories"] * INTAKE_MULTIPLIER
    kf = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    epochs = pd.read_csv(ROOT / "intake" / "diet_epochs.csv", parse_dates=["start", "end"])
    seg = pd.read_csv(ROOT / "analysis" / "CP_segments.csv", parse_dates=["start", "end"])

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

    # ── identify transitions: between consecutive plateaus, or non-plateau segs ──
    # Use CP segmentation. Bridge consecutive plateaus by the segment between them.
    # Treat non-plateau segments themselves as transitions.

    # Compose list of transition events: ordered pairs of (plateau_before, plateau_after)
    # and capture the "transition zone" between them.
    transitions = []
    plateaus = seg[seg["is_plateau"]].sort_values("start").reset_index(drop=True)
    for i in range(len(plateaus) - 1):
        before = plateaus.iloc[i]
        after = plateaus.iloc[i + 1]
        transitions.append({
            "name": f"P{int(before['idx'])} ({before['fm_mean']:.0f}) → P{int(after['idx'])} ({after['fm_mean']:.0f})",
            "before_idx": int(before["idx"]),
            "after_idx": int(after["idx"]),
            "before_fm": before["fm_mean"],
            "after_fm": after["fm_mean"],
            "delta_fm": after["fm_mean"] - before["fm_mean"],
            "trans_start": before["end"],
            "trans_end": after["start"],
        })

    # Also include leading transition into P2 (the floor) — the willpower descent
    # from P0 (~80 lb mean) all the way down. Captured by including segments 0 and 1
    # as the transition into P2.
    if len(plateaus) > 0:
        first_p = plateaus.iloc[0]
        # Use the start of the data
        first_date = pre["date"].iloc[0]
        transitions.insert(0, {
            "name": f"start → P{int(first_p['idx'])} ({first_p['fm_mean']:.0f})",
            "before_idx": -1,
            "after_idx": int(first_p["idx"]),
            "before_fm": float(pre.iloc[0]["fm"]),
            "after_fm": first_p["fm_mean"],
            "delta_fm": first_p["fm_mean"] - float(pre.iloc[0]["fm"]),
            "trans_start": first_date,
            "trans_end": first_p["start"],
        })

    # ── PART 1: Per-transition characterization ──
    print("=" * 110)
    print("PART 1 — TRANSITION ANATOMY (per-transition stats)")
    print("=" * 110)
    print(f"  {'transition':<32} {'duration':>8} {'ΔFM':>6} {'rate':>9} "
          f"{'mean surp':>10} {'cum surp':>10} {'binge%':>7} {'TDEE':>6}")
    for t in transitions:
        m = (pre["date"] >= t["trans_start"]) & (pre["date"] <= t["trans_end"])
        g = pre[m]
        if len(g) < 7:
            continue
        dur_days = (t["trans_end"] - t["trans_start"]).days
        rate = t["delta_fm"] / max(dur_days, 1)
        mean_surplus = float(g["surplus"].mean())
        cum_surplus = float(g["surplus"].sum())
        binge_pct = float(g["binge"].mean() * 100)
        tdee_mean = float(g["tdee"].mean())
        t["dur_days"] = dur_days
        t["rate"] = rate
        t["mean_surplus"] = mean_surplus
        t["cum_surplus"] = cum_surplus
        t["binge_pct"] = binge_pct
        t["tdee_mean"] = tdee_mean
        print(f"  {t['name']:<32} {dur_days:>5d} d {t['delta_fm']:>+6.1f} {rate:>+9.4f} "
              f"{mean_surplus:>+10.0f} {cum_surplus:>+10.0f} {binge_pct:>6.1f}% {tdee_mean:>6.0f}")

    # ── PART 2: Absorption vs breakthrough ──
    print("\n" + "=" * 110)
    print("PART 2 — ABSORPTION vs BREAKTHROUGH")
    print("=" * 110)
    print("\nFor each plateau, what was the LARGEST cumulative surplus over rolling")
    print("7/14/30/60/90-day windows that the plateau absorbed without breaking?")
    print()

    plateau_max = []
    for _, p in plateaus.iterrows():
        m = (pre["date"] >= p["start"]) & (pre["date"] <= p["end"])
        g = pre[m]
        record = {"plateau": f"P{int(p['idx'])}", "fm_mean": p["fm_mean"], "n_days": len(g)}
        for w in (7, 14, 30, 60, 90):
            roll = g["surplus"].rolling(w, min_periods=w).sum()
            record[f"max_cum_surp_{w}d"] = float(roll.max()) if roll.notna().any() else float("nan")
            record[f"min_cum_surp_{w}d"] = float(roll.min()) if roll.notna().any() else float("nan")
        plateau_max.append(record)
    pmax = pd.DataFrame(plateau_max)
    print(f"  {'plateau':<8} {'FM':>5} {'n_days':>6}  "
          f"{'+7d':>7} {'+14d':>7} {'+30d':>8} {'+60d':>9} {'+90d':>9} | "
          f"{'-7d':>7} {'-14d':>7} {'-30d':>8} {'-60d':>9} {'-90d':>9}")
    for _, r in pmax.iterrows():
        print(f"  {r['plateau']:<8} {r['fm_mean']:>5.0f} {int(r['n_days']):>6}  "
              f"{r['max_cum_surp_7d']:>+7.0f} {r['max_cum_surp_14d']:>+7.0f} "
              f"{r['max_cum_surp_30d']:>+8.0f} {r['max_cum_surp_60d']:>+9.0f} {r['max_cum_surp_90d']:>+9.0f} | "
              f"{r['min_cum_surp_7d']:>+7.0f} {r['min_cum_surp_14d']:>+7.0f} "
              f"{r['min_cum_surp_30d']:>+8.0f} {r['min_cum_surp_60d']:>+9.0f} {r['min_cum_surp_90d']:>+9.0f}")

    print("\nFor each transition, what was the cumulative surplus during it?")
    print(f"\n  {'transition':<32} {'days':>5}  "
          f"{'cum surp':>10} {'cum / days':>10} {'best 30d':>10} {'best 90d':>10}")
    for t in transitions:
        if "dur_days" not in t:
            continue
        m = (pre["date"] >= t["trans_start"]) & (pre["date"] <= t["trans_end"])
        g = pre[m]
        if len(g) < 30:
            continue
        cum = float(g["surplus"].sum())
        roll30 = g["surplus"].rolling(30, min_periods=30).sum()
        roll90 = g["surplus"].rolling(90, min_periods=90).sum()
        best30 = float(roll30.max()) if roll30.notna().any() else float("nan")
        best90 = float(roll90.max()) if roll90.notna().any() else float("nan")
        print(f"  {t['name']:<32} {t['dur_days']:>5} {cum:>+10.0f} "
              f"{cum / max(t['dur_days'], 1):>+10.0f} {best30:>+10.0f} {best90:>+10.0f}")

    # ── PART 3: Inside the transitions ──
    print("\n" + "=" * 110)
    print("PART 3 — INSIDE EACH TRANSITION (smooth ramp, sub-plateaus, or bursts?)")
    print("=" * 110)
    for t in transitions:
        if "dur_days" not in t:
            continue
        m = (pre["date"] >= t["trans_start"]) & (pre["date"] <= t["trans_end"])
        g = pre[m].copy().reset_index(drop=True)
        if len(g) < 30:
            continue
        print(f"\n  {t['name']}  [{t['trans_start'].strftime('%Y-%m-%d')} → "
              f"{t['trans_end'].strftime('%Y-%m-%d')}, {t['dur_days']}d, ΔFM {t['delta_fm']:+.1f} lb]")

        # Identify the periods of fastest FM movement: 30-day rolling FM slope
        g["fm_slope_30d"] = (g["fm"] - g["fm"].shift(30)) / 30.0
        slope_max = float(g["fm_slope_30d"].abs().max())
        # Identify when the steepest move happened
        if g["fm_slope_30d"].notna().any():
            idx_max = int(g["fm_slope_30d"].abs().idxmax())
            steep_date = g["date"].iloc[idx_max]
            steep_fm = g["fm"].iloc[idx_max]
            print(f"    Steepest 30-d slope: {g['fm_slope_30d'].iloc[idx_max]:+.4f} lb/day "
                  f"around {steep_date.strftime('%Y-%m-%d')} (FM {steep_fm:.1f})")

        # Bin transition into 4 quarters and report mean intake/surplus per quarter
        n = len(g)
        chunks = [g.iloc[i * n // 4: (i + 1) * n // 4] for i in range(4)]
        print(f"    {'quarter':<8} {'days':>5} {'FM start':>9} {'FM end':>7} "
              f"{'intake':>7} {'TDEE':>6} {'mean surp':>10} {'binge%':>7}")
        for qi, q in enumerate(chunks):
            if len(q) < 5:
                continue
            print(f"    Q{qi + 1:<7} {len(q):>5} {q['fm'].iloc[0]:>9.1f} "
                  f"{q['fm'].iloc[-1]:>7.1f} {q['calories'].mean():>7.0f} "
                  f"{q['tdee'].mean():>6.0f} {(q['calories'] - q['tdee']).mean():>+10.0f} "
                  f"{q['binge'].mean() * 100:>6.1f}%")

        # Identify diet epochs that overlap this transition
        overlap = epochs[(epochs["start"] <= t["trans_end"]) & (epochs["end"] >= t["trans_start"])]
        if len(overlap):
            print(f"    Overlapping diet epochs:")
            for _, e in overlap.iterrows():
                print(f"      {e['start'].date()} → {e['end'].date()}: {e['label']} ({e['detail']})")

    # ── Plot per-transition zoom ──
    n_trans = len([t for t in transitions if "dur_days" in t and t["dur_days"] > 30])
    cols = 2
    rows = (n_trans + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(15, 3.2 * rows), squeeze=False)
    axflat = axes.flatten()
    plot_idx = 0
    for t in transitions:
        if "dur_days" not in t or t["dur_days"] <= 30:
            continue
        ax = axflat[plot_idx]
        plot_idx += 1
        # Pad the window for context
        pad = pd.Timedelta(days=60)
        m = (pre["date"] >= t["trans_start"] - pad) & (pre["date"] <= t["trans_end"] + pad)
        g = pre[m].copy()
        ax2 = ax.twinx()
        ax.plot(g["date"], g["fm"], color="black", lw=1.0, label="FM")
        cum = (g["surplus"]).cumsum() / 3500.0  # in lb units
        ax2.plot(g["date"], cum - cum.iloc[0], color="tab:purple", lw=0.8,
                 alpha=0.6, label="cum surplus / 3500 (lb)")
        ax.axvspan(t["trans_start"], t["trans_end"], color="#f5b7b1", alpha=0.3)
        ax.axhline(t["before_fm"], color="tab:blue", ls="--", lw=1, alpha=0.6)
        ax.axhline(t["after_fm"], color="tab:red", ls="--", lw=1, alpha=0.6)
        # Mark binge days
        binges = g[g["binge"] == 1]
        ax.scatter(binges["date"], binges["fm"], s=8, color="tab:orange",
                   alpha=0.5, label="binge")
        ax.set_title(f"{t['name']}  ({t['dur_days']}d, ΔFM {t['delta_fm']:+.1f})", fontsize=10)
        ax.set_ylabel("FM (lb)")
        ax2.set_ylabel("cum surp (lb)", color="tab:purple")
        ax.grid(alpha=0.3)

    for j in range(plot_idx, len(axflat)):
        axflat[j].set_visible(False)

    fig.tight_layout()
    plot_path = ROOT / "analysis" / "CQ_transition_anatomy.png"
    fig.savefig(plot_path, dpi=140)
    print(f"\nSaved {plot_path}")


if __name__ == "__main__":
    main()
