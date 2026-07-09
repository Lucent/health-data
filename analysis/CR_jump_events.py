"""CR. Find the discrete FM jumps directly. No segmentation prior.

For each day, compute a 30-day forward FM slope. Identify contiguous periods
where the slope exceeds a threshold (rising or falling). These are "jump
events". Characterise each: duration, ΔFM, rate, intake/TDEE/binge during
the jump, what happened just before, what diet epoch overlaps.

The discovery question: are FM jumps discrete events with idle periods
between them, or a continuous gradient? If discrete, what triggers them?
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

    # 30-day backward + forward symmetric slope at each day
    pre["fm_slope_30d"] = (pre["fm"].shift(-15) - pre["fm"].shift(15)) / 30.0

    # Identify jump events: contiguous days with |slope| above threshold
    THRESHOLD = 0.04   # lb/day sustained
    MIN_DURATION = 21  # days
    pre["above_thresh"] = pre["fm_slope_30d"].abs() > THRESHOLD

    runs = []
    in_run = False
    start_i = None
    for i, val in enumerate(pre["above_thresh"]):
        if val and not in_run:
            in_run = True
            start_i = i
        elif (not val or i == len(pre) - 1) and in_run:
            end_i = i if not val else i + 1
            if end_i - start_i >= MIN_DURATION:
                runs.append((start_i, end_i))
            in_run = False

    print("=" * 110)
    print(f"DISCRETE FM JUMP EVENTS  (|30-d slope| > {THRESHOLD} lb/d, sustained ≥ {MIN_DURATION} d)")
    print("=" * 110)
    print(f"  {'#':>3} {'start':<11} {'end':<11} {'days':>5} "
          f"{'FM start':>9} {'FM end':>7} {'ΔFM':>6} {'rate':>9} "
          f"{'intake':>7} {'TDEE':>6} {'mean surp':>10} {'binge%':>7} | overlap")

    jumps = []
    for ji, (s, e) in enumerate(runs):
        g = pre.iloc[s:e]
        if len(g) < MIN_DURATION:
            continue
        fm_start = float(g["fm"].iloc[0])
        fm_end = float(g["fm"].iloc[-1])
        delta = fm_end - fm_start
        rate = delta / max(len(g), 1)
        intake_mean = float(g["calories"].mean())
        tdee_mean = float(g["tdee"].mean())
        surplus_mean = float(g["surplus"].mean())
        binge_pct = float(g["binge"].mean() * 100)
        # Overlapping diet epochs
        overlap = epochs[(epochs["start"] <= g["date"].iloc[-1]) &
                          (epochs["end"] >= g["date"].iloc[0])]
        overlap_str = ", ".join(overlap["label"].tolist()) if len(overlap) else "—"
        jumps.append({
            "idx": ji,
            "start": g["date"].iloc[0],
            "end": g["date"].iloc[-1],
            "days": len(g),
            "fm_start": fm_start,
            "fm_end": fm_end,
            "delta": delta,
            "rate": rate,
            "intake": intake_mean,
            "tdee": tdee_mean,
            "surplus": surplus_mean,
            "binge_pct": binge_pct,
            "overlap": overlap_str,
        })
        print(f"  {ji:>3} {g['date'].iloc[0].strftime('%Y-%m-%d'):<11} "
              f"{g['date'].iloc[-1].strftime('%Y-%m-%d'):<11} {len(g):>5} "
              f"{fm_start:>9.1f} {fm_end:>7.1f} {delta:>+6.1f} {rate:>+9.4f} "
              f"{intake_mean:>7.0f} {tdee_mean:>6.0f} {surplus_mean:>+10.0f} "
              f"{binge_pct:>6.1f}% | {overlap_str[:50]}")

    jdf = pd.DataFrame(jumps)
    out_csv = ROOT / "analysis" / "CR_jump_events.csv"
    jdf.to_csv(out_csv, index=False)
    print(f"\nSaved {out_csv}")

    if not jumps:
        return

    # ── Aggregate stats ──
    print("\n" + "=" * 110)
    print("DISTRIBUTION")
    print("=" * 110)
    rises = jdf[jdf["delta"] > 0]
    falls = jdf[jdf["delta"] < 0]
    print(f"  Rising jumps:  n={len(rises)}, median ΔFM = {rises['delta'].median():+.1f} lb, "
          f"median duration = {rises['days'].median():.0f} d, "
          f"median rate = {rises['rate'].median():+.4f} lb/day")
    if len(falls):
        print(f"  Falling jumps: n={len(falls)}, median ΔFM = {falls['delta'].median():+.1f} lb, "
              f"median duration = {falls['days'].median():.0f} d, "
              f"median rate = {falls['rate'].median():+.4f} lb/day")

    # ── Plot ──
    fig, axes = plt.subplots(3, 1, figsize=(15, 9), sharex=True)
    pre["intake_30d"] = pre["calories"].rolling(30, min_periods=15).mean()
    pre["tdee_30d"] = pre["tdee"].rolling(30, min_periods=15).mean()
    pre["surplus_30d"] = pre["surplus"].rolling(30, min_periods=15).mean()

    ax = axes[0]
    ax.plot(pre["date"], pre["fm"], color="black", lw=0.8)
    for j in jumps:
        c = "#ff6b6b" if j["delta"] > 0 else "#4ecdc4"
        ax.axvspan(j["start"], j["end"], alpha=0.3, color=c)
        # annotate ΔFM
        mid = j["start"] + (j["end"] - j["start"]) / 2
        ax.annotate(f"{j['delta']:+.0f}", (mid, (j["fm_start"] + j["fm_end"]) / 2),
                    ha="center", fontsize=8, fontweight="bold")
    ax.set_ylabel("FM (lb)")
    ax.set_title(f"Discrete FM jumps (|30d slope| > {THRESHOLD} lb/d, ≥ {MIN_DURATION} d). "
                 f"Rises = red, falls = teal")
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(pre["date"], pre["intake_30d"], color="tab:green", label="intake (30d)")
    ax.plot(pre["date"], pre["tdee_30d"], color="tab:red", label="TDEE (30d)")
    for j in jumps:
        c = "#ff6b6b" if j["delta"] > 0 else "#4ecdc4"
        ax.axvspan(j["start"], j["end"], alpha=0.2, color=c)
    ax.set_ylabel("cal/day")
    ax.legend()
    ax.grid(alpha=0.3)

    ax = axes[2]
    ax.plot(pre["date"], pre["surplus_30d"], color="tab:purple")
    ax.axhline(0, color="black", lw=0.5)
    for j in jumps:
        c = "#ff6b6b" if j["delta"] > 0 else "#4ecdc4"
        ax.axvspan(j["start"], j["end"], alpha=0.2, color=c)
    ax.set_ylabel("surplus (30d, cal/day)")
    ax.set_xlabel("Date")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    plot_path = ROOT / "analysis" / "CR_jump_events.png"
    fig.savefig(plot_path, dpi=140)
    print(f"\nSaved {plot_path}")


if __name__ == "__main__":
    main()
