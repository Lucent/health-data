"""CO. What happens between maintenance plateaus?

CN found maintenance FM rises 30 → 55 → 72 lb over 15 years in monotone steps.
The body idles at one FM for months/years, then ratchets up to a new idle.
Look at the actual daily numbers across each transition vs the plateaus that
bracket it. Specifically:

  - When FM is idling, what are intake / TDEE / binge / surplus running at?
  - During the transition, what changes? Is it intake going up? TDEE going down?
    Both? Concentrated in days/weeks/months?
  - After the transition, what settles into the new idle? Did intake drop back?
    Did TDEE recover? Or did the equilibrium just move?

Output: per-period table of the actual numbers, plot of FM with plateaus and
transitions colored, and a focused look at each transition's intake/TDEE.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

ROOT = Path(__file__).resolve().parent.parent
INTAKE_MULTIPLIER = 1.09
BINGE_THRESHOLD = 1000


def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    intake["calories"] = intake["calories"] * INTAKE_MULTIPLIER
    kf = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    p2 = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])

    df = intake[["date", "calories"]].merge(
        kf[["date", "fat_mass_lbs", "tdee", "lm_lbs"]], on="date")
    df = df.merge(p2[["date", "expected_rmr"]], on="date", how="left")
    df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df = df.sort_values("date").reset_index(drop=True)
    df["fm"] = df["fat_mass_lbs"]
    df["binge"] = (df["calories"] > df["tdee"] + BINGE_THRESHOLD).astype(int)
    df["surplus"] = df["calories"] - df["tdee"]
    df["tdee_rmr_ratio"] = df["tdee"] / df["expected_rmr"]

    pre = df[df["effective_level"] == 0].copy()

    # Visually-anchored period definitions, from inspecting FM trajectory
    # and CN's maintenance distribution
    periods = [
        ("Plateau-1: post-aggressive-loss",   "2013-07-01", "2014-06-30"),  # FM ~22
        ("Transition-1: post-floor regain",   "2014-07-01", "2015-06-30"),  # 22 → 36
        ("Plateau-2: enlightened era",        "2015-07-01", "2017-12-31"),  # FM ~50
        ("Transition-2: keto era",            "2018-01-01", "2019-12-31"),  # 51 → 68
        ("Plateau-3: post-keto stable",       "2020-01-01", "2021-12-31"),  # FM ~71
        ("Transition-3: pre-tirz drift",      "2022-01-01", "2023-12-31"),  # 74 → 86
        ("Plateau-4: pre-tirz peak",          "2023-12-01", "2024-04-30"),  # FM ~85
    ]

    print("=" * 100)
    print(f"  {'Period':<35} {'FM start':>9} {'FM end':>8} {'ΔFM':>6} "
          f"{'intake':>7} {'TDEE':>6} {'surplus':>8} {'TDEE/RMR':>9} "
          f"{'binge%':>7} {'days':>5}")
    print("=" * 100)

    rows = []
    for label, start, end in periods:
        mask = (pre["date"] >= start) & (pre["date"] <= end)
        g = pre[mask]
        if len(g) < 30:
            continue
        fm_start = float(g["fm"].iloc[0])
        fm_end = float(g["fm"].iloc[-1])
        intake_mean = float(g["calories"].mean())
        tdee_mean = float(g["tdee"].mean())
        surplus_mean = intake_mean - tdee_mean
        ratio_mean = float(g["tdee_rmr_ratio"].mean())
        binge_pct = float(g["binge"].mean() * 100)
        rows.append({
            "label": label, "start": start, "end": end,
            "fm_start": fm_start, "fm_end": fm_end,
            "intake": intake_mean, "tdee": tdee_mean,
            "surplus": surplus_mean, "ratio": ratio_mean,
            "binge_pct": binge_pct, "n": len(g),
        })
        print(f"  {label:<35} {fm_start:>9.1f} {fm_end:>8.1f} {fm_end-fm_start:>+6.1f} "
              f"{intake_mean:>7.0f} {tdee_mean:>6.0f} {surplus_mean:>+8.0f} {ratio_mean:>9.3f} "
              f"{binge_pct:>6.1f}% {len(g):>5d}")

    summary = pd.DataFrame(rows)

    # ── Now compare adjacent plateaus' steady-state intake and TDEE ──
    print("\n" + "=" * 100)
    print("ADJACENT PLATEAU COMPARISON — what shifted, intake or TDEE?")
    print("=" * 100)
    plateaus = summary[summary["label"].str.startswith("Plateau")].reset_index(drop=True)
    for i in range(len(plateaus) - 1):
        a, b = plateaus.iloc[i], plateaus.iloc[i + 1]
        d_intake = b["intake"] - a["intake"]
        d_tdee   = b["tdee"]   - a["tdee"]
        d_fm     = b["fm_start"] - a["fm_end"]
        d_ratio  = b["ratio"]  - a["ratio"]
        print(f"\n  {a['label']:<35} → {b['label']}")
        print(f"    FM:        {a['fm_end']:>5.1f} lb → {b['fm_start']:>5.1f} lb   ΔFM = {d_fm:+.1f} lb")
        print(f"    Intake:    {a['intake']:>5.0f}    → {b['intake']:>5.0f}      Δ   = {d_intake:+.0f} cal/d")
        print(f"    TDEE:      {a['tdee']:>5.0f}    → {b['tdee']:>5.0f}      Δ   = {d_tdee:+.0f} cal/d")
        print(f"    TDEE/RMR:  {a['ratio']:>5.3f}    → {b['ratio']:>5.3f}      Δ   = {d_ratio:+.3f}")
        print(f"    binge %:   {a['binge_pct']:>5.1f}    → {b['binge_pct']:>5.1f}")

    # ── During transitions, where does the surplus accumulate? ──
    print("\n" + "=" * 100)
    print("WITHIN-TRANSITION DAILY DETAIL — intake_t vs TDEE_t")
    print("=" * 100)
    transitions = summary[summary["label"].str.startswith("Transition")].reset_index(drop=True)
    for _, t in transitions.iterrows():
        mask = (pre["date"] >= t["start"]) & (pre["date"] <= t["end"])
        g = pre[mask].copy()
        # Bin into thirds
        n = len(g)
        t1 = g.iloc[:n // 3]
        t2 = g.iloc[n // 3: 2 * n // 3]
        t3 = g.iloc[2 * n // 3:]
        print(f"\n  {t['label']}")
        print(f"    {'segment':<12} {'days':>5} {'intake':>8} {'TDEE':>7} {'surplus':>9} {'binge%':>7} {'TDEE/RMR':>9} {'mean_FM':>8}")
        for name, sub in [("first 1/3", t1), ("middle", t2), ("last 1/3", t3)]:
            print(f"    {name:<12} {len(sub):>5d} {sub['calories'].mean():>8.0f} "
                  f"{sub['tdee'].mean():>7.0f} {(sub['calories'] - sub['tdee']).mean():>+9.0f} "
                  f"{sub['binge'].mean() * 100:>6.1f}% {sub['tdee_rmr_ratio'].mean():>9.3f} "
                  f"{sub['fm'].mean():>8.1f}")

    # ── Plot: FM trajectory with periods colored, plus intake / TDEE / surplus rolling ──
    fig, axes = plt.subplots(3, 1, figsize=(14, 11), sharex=True)
    pre30 = pre.copy()
    pre30["intake_30d"] = pre30["calories"].rolling(30, min_periods=15).mean()
    pre30["tdee_30d"] = pre30["tdee"].rolling(30, min_periods=15).mean()
    pre30["surplus_30d"] = pre30["surplus"].rolling(30, min_periods=15).mean()

    ax = axes[0]
    ax.plot(pre30["date"], pre30["fm"], color="black", lw=1.0)
    colors = ["#cce5ff", "#ffcccc", "#cce5ff", "#ffcccc", "#cce5ff", "#ffcccc", "#cce5ff"]
    for (label, start, end), c in zip(periods, colors):
        ax.axvspan(pd.Timestamp(start), pd.Timestamp(end), color=c, alpha=0.5)
    ax.set_ylabel("FM (lb)")
    ax.set_title("FM with plateaus (blue) and transitions (red)")
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(pre30["date"], pre30["intake_30d"], label="intake (30d)", color="tab:green")
    ax.plot(pre30["date"], pre30["tdee_30d"], label="TDEE (30d)", color="tab:red")
    for (label, start, end), c in zip(periods, colors):
        ax.axvspan(pd.Timestamp(start), pd.Timestamp(end), color=c, alpha=0.3)
    ax.set_ylabel("cal/day")
    ax.legend()
    ax.grid(alpha=0.3)

    ax = axes[2]
    ax.plot(pre30["date"], pre30["surplus_30d"], color="tab:purple")
    ax.axhline(0, color="black", lw=0.5)
    for (label, start, end), c in zip(periods, colors):
        ax.axvspan(pd.Timestamp(start), pd.Timestamp(end), color=c, alpha=0.3)
    ax.set_ylabel("surplus (30d, cal/day)")
    ax.set_xlabel("Date")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    plot_path = ROOT / "analysis" / "CO_plateau_transitions.png"
    fig.savefig(plot_path, dpi=140)
    print(f"\nSaved {plot_path}")


if __name__ == "__main__":
    main()
