"""CM. Is AG's set-point-distance signal really FM velocity in disguise?

Mathematical motivation. At HL=60d, EMA(FM) is approximately a 43-day lag of FM
(since EMA's "centre of mass" sits HL · 1/ln(2) ≈ 43d back). So:

    sp_distance(t) = EMA(FM, 60d)(t) − FM(t)
                   ≈ FM(t − 43d) − FM(t)
                   ≈ −43 · dFM/dt          (for slow trends)

If AG's r=−0.62 is essentially a velocity signal, then:
  - r(binge_rate, sp_distance) ≈ r(binge_rate, −velocity) [direction]
  - Corr(sp_distance, −velocity_~43d) ≈ 1
  - Joint regression: only ONE of them retains a coefficient
  - FM-banded analysis: within FM bands, only velocity matters

If there's something beyond velocity (a real history-dependent SP):
  - Joint regression: BOTH retain coefficients
  - sp_distance partial r > velocity partial r at long lags
  - FM-banded relationship persists for sp_distance after partialing velocity

This script does the discriminative test.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
INTAKE_MULTIPLIER = 1.09
BINGE_THRESHOLD = 1000


def partial_corr(y, x_target, x_controls):
    """Partial correlation of y with x_target controlling for x_controls (cols)."""
    X_c = np.column_stack([np.ones(len(y))] + list(x_controls))
    # residualize y
    by, *_ = np.linalg.lstsq(X_c, y, rcond=None)
    y_resid = y - X_c @ by
    # residualize x_target
    bx, *_ = np.linalg.lstsq(X_c, x_target, rcond=None)
    x_resid = x_target - X_c @ bx
    return float(np.corrcoef(y_resid, x_resid)[0, 1])


def main():
    print("=" * 80)
    print("CM. Velocity vs set-point distance — which actually predicts binges?")
    print("=" * 80)

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
    df["binge"] = (df["calories"] > df["tdee"] + BINGE_THRESHOLD).astype(int)
    df["binge_rate_90d"] = df["binge"].rolling(90, min_periods=90).mean()

    # Pre-tirz only
    pre = df[df["effective_level"] == 0].copy().reset_index(drop=True)

    # Velocity at multiple windows (lb/day)
    for w in (14, 30, 43, 60, 90):
        pre[f"vel_{w}d"] = (pre["fm"] - pre["fm"].shift(w)) / w

    # SP distance at multiple HL
    fm_s = pre.set_index("date")["fm"]
    for hl in (30, 60, 100, 180):
        ema = fm_s.ewm(halflife=hl, min_periods=30).mean()
        pre[f"sp_dist_{hl}"] = pre["date"].map(ema.to_dict()) - pre["fm"]

    # Phase tag
    def phase_of(d):
        for _, row in epochs.iterrows():
            if row["start"] <= d <= row["end"]:
                return row["label"]
        return "other"
    pre["phase"] = pre["date"].apply(phase_of)

    valid = pre.dropna(subset=["binge_rate_90d", "vel_30d", "sp_dist_60",
                               "vel_60d", "vel_43d", "sp_dist_30"]).copy()
    print(f"\nValid pre-tirz days: {len(valid)}")

    # ── 1. Mathematical equivalence check ──
    print("\n--- Mathematical equivalence: sp_distance vs -velocity ---")
    print("  if sp_dist_HL ≈ -velocity_(0.69·HL), correlation should be ≈ +1")
    for hl in (30, 60, 100, 180):
        for w in (14, 30, 43, 60, 90):
            r = float(np.corrcoef(valid[f"sp_dist_{hl}"], -valid[f"vel_{w}d"])[0, 1])
            print(f"  Corr(sp_dist_HL{hl}, -vel_{w}d) = {r:+.4f}")

    # ── 2. Marginal correlations with binge_rate_90d ──
    print("\n--- Marginal correlations with 90d binge rate ---")
    targets = []
    for w in (14, 30, 43, 60, 90):
        c = f"vel_{w}d"
        r = float(np.corrcoef(valid[c], valid["binge_rate_90d"])[0, 1])
        print(f"  binge_rate vs vel_{w}d        r = {r:+.4f}")
        targets.append((c, r))
    for hl in (30, 60, 100, 180):
        c = f"sp_dist_{hl}"
        r = float(np.corrcoef(valid[c], valid["binge_rate_90d"])[0, 1])
        print(f"  binge_rate vs sp_dist_{hl:>3}d  r = {r:+.4f}")
        targets.append((c, r))
    r_fm = float(np.corrcoef(valid["fm"], valid["binge_rate_90d"])[0, 1])
    print(f"  binge_rate vs FM (absolute)  r = {r_fm:+.4f}")

    # ── 3. Joint partial correlations ──
    print("\n--- Joint test: who survives partialing the other? ---")
    # SP-distance vs velocity, controlling for absolute FM
    y = valid["binge_rate_90d"].values
    fm_v = valid["fm"].values
    for hl, w in [(60, 43), (60, 30), (100, 60), (30, 14)]:
        sp = valid[f"sp_dist_{hl}"].values
        vel = valid[f"vel_{w}d"].values
        r_sp_alone   = partial_corr(y, sp, [fm_v])
        r_vel_alone  = partial_corr(y, vel, [fm_v])
        r_sp_partial = partial_corr(y, sp, [fm_v, vel])
        r_vel_partial = partial_corr(y, vel, [fm_v, sp])
        print(f"  HL={hl}d  vel={w}d:")
        print(f"    sp_dist  alone (|FM)         r = {r_sp_alone:+.4f}")
        print(f"    vel      alone (|FM)         r = {r_vel_alone:+.4f}")
        print(f"    sp_dist  partial (|FM, vel)  r = {r_sp_partial:+.4f}")
        print(f"    vel      partial (|FM, sp)   r = {r_vel_partial:+.4f}")

    # ── 4. Phase-stratified ──
    print("\n--- Phase-stratified marginal r (sp_dist_60, vel_30d) ---")
    print(f"  {'phase':<26} {'n':>5} {'r(sp_60)':>10} {'r(vel_30)':>10}")
    phase_results = []
    for phase, g in valid.groupby("phase"):
        if len(g) < 100:
            continue
        r_sp = float(np.corrcoef(g["sp_dist_60"], g["binge_rate_90d"])[0, 1])
        r_v  = float(np.corrcoef(g["vel_30d"],   g["binge_rate_90d"])[0, 1])
        print(f"  {phase:<26} {len(g):>5} {r_sp:>+10.3f} {r_v:>+10.3f}")
        phase_results.append({"phase": phase, "n": len(g), "r_sp": r_sp, "r_vel": r_v})

    # ── 5. FM-banded ──
    print("\n--- FM-banded marginal r (controls for absolute FM by stratification) ---")
    print(f"  {'FM band (lb)':<14} {'n':>5} {'r(sp_60)':>10} {'r(vel_30)':>10} {'r(vel_43)':>10}")
    bands = [(15, 30), (30, 50), (50, 70), (70, 100)]
    for lo, hi in bands:
        g = valid[(valid["fm"] >= lo) & (valid["fm"] < hi)]
        if len(g) < 100:
            continue
        r_sp = float(np.corrcoef(g["sp_dist_60"], g["binge_rate_90d"])[0, 1])
        r_v30  = float(np.corrcoef(g["vel_30d"],  g["binge_rate_90d"])[0, 1])
        r_v43  = float(np.corrcoef(g["vel_43d"],  g["binge_rate_90d"])[0, 1])
        print(f"  {lo:>3}-{hi:<3}        {len(g):>5} {r_sp:>+10.3f} {r_v30:>+10.3f} {r_v43:>+10.3f}")

    # ── 6. Plot ──
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    valid["year_phase"] = valid["date"].dt.year

    ax = axes[0, 0]
    ax.scatter(valid["sp_dist_60"], valid["binge_rate_90d"], s=2,
               c=valid["year_phase"], cmap="viridis", alpha=0.4)
    ax.set_xlabel("sp_distance (HL=60d)  [SP − FM, lb]")
    ax.set_ylabel("90d binge rate")
    r = float(np.corrcoef(valid["sp_dist_60"], valid["binge_rate_90d"])[0, 1])
    ax.set_title(f"Binge rate vs SP-distance  (r = {r:+.3f})")
    ax.axvline(0, color="black", lw=0.5)

    ax = axes[0, 1]
    ax.scatter(valid["vel_30d"], valid["binge_rate_90d"], s=2,
               c=valid["year_phase"], cmap="viridis", alpha=0.4)
    ax.set_xlabel("FM velocity (Δ30d, lb/day)")
    ax.set_ylabel("90d binge rate")
    r = float(np.corrcoef(valid["vel_30d"], valid["binge_rate_90d"])[0, 1])
    ax.set_title(f"Binge rate vs FM velocity (Δ30d) (r = {r:+.3f})")
    ax.axvline(0, color="black", lw=0.5)

    ax = axes[1, 0]
    ax.scatter(valid["sp_dist_60"], -valid["vel_43d"], s=2, alpha=0.4)
    ax.set_xlabel("sp_distance (HL=60d)")
    ax.set_ylabel("−velocity (Δ43d)")
    r = float(np.corrcoef(valid["sp_dist_60"], -valid["vel_43d"])[0, 1])
    ax.set_title(f"Equivalence test: sp_dist vs −velocity  (r = {r:+.3f})")
    lim = max(valid["sp_dist_60"].abs().max(), 0.05) * 0.5
    ax.plot([-lim, lim], [-lim/43, lim/43], "k--", alpha=0.5, label="theoretical: −1/43 slope")
    ax.legend()

    ax = axes[1, 1]
    if phase_results:
        pr = pd.DataFrame(phase_results)
        x = np.arange(len(pr))
        ax.bar(x - 0.2, pr["r_sp"], 0.4, label="r(sp_dist_60)", color="tab:blue")
        ax.bar(x + 0.2, pr["r_vel"], 0.4, label="r(vel_30d)", color="tab:orange")
        ax.set_xticks(x)
        ax.set_xticklabels([p[:12] for p in pr["phase"]], rotation=45, ha="right", fontsize=8)
        ax.axhline(0, color="black", lw=0.5)
        ax.set_ylabel("r with 90d binge rate")
        ax.set_title("Per-phase marginal r")
        ax.legend()

    fig.tight_layout()
    plot_path = ROOT / "analysis" / "CM_velocity_vs_setpoint.png"
    fig.savefig(plot_path, dpi=140)
    print(f"\nSaved {plot_path}")


if __name__ == "__main__":
    main()
