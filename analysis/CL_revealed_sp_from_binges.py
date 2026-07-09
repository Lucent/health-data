"""CL. Revealed set point from binge-rate inversion.

# Premise

AG fits binge_rate = f(SP - FM) where SP = EMA(FM, HL=60d). The autocorrelation
worry: high-binge periods raise FM, FM smooths to EMA-SP, EMA-SP creates apparent
SP-distance. AG could be measuring binge momentum, not a latent set point.

# Inversion

If f is monotone, given observed binge_rate(t) and FM(t):

    d(t) = SP(t) - FM(t) = f^{-1}(binge_rate(t))
    SP(t) = FM(t) + d(t)

This computes SP without assuming any FM smoothing. SP gets to be whatever number,
when subtracted from FM, yields the right binge propensity for that day.

# What success vs failure looks like (preregistered)

If binge-derived SP is a real latent state independent of FM smoothing:
  1. Corr(SP_revealed, EMA-SP) is MODERATE: 0.50 to 0.90. Not perfect.
  2. SP and FM decouple in identifiable periods (e.g., post-restriction, drug start).
  3. SP_revealed - FM responds to exogenous events (binge episode, fast) ahead of FM.

If we're just rediscovering EMA-SP via a different smoothing:
  - Corr(SP_revealed, EMA-SP) > 0.95.
  - SP_revealed is a noisier version of EMA-SP with no novel structure.

If "binge momentum" explains AG (the user's worry):
  - SP_revealed - FM tracks lagged binge density and nothing else.
  - Residualizing SP_revealed on past binge counts collapses the signal.

# Construction note

f is fit using AG's binned (SP_dist_EMA, binge_rate) curve as a calibration.
This is mildly circular (f's shape comes from EMA-SP) but the *inversion* uses
only daily binge rate and FM — no FM smoothing in SP_revealed itself. The shape
of f is borrowed; the daily values of SP_revealed are independent.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parent.parent

INTAKE_MULTIPLIER = 1.09
BINGE_THRESHOLD = 1000  # cal above TDEE


def logistic(d, base, amp, k, d0):
    return base + amp / (1.0 + np.exp(-k * (d - d0)))


def invert(b, base, amp, k, d0):
    """f^{-1}: binge_rate → SP_dist (positive when below SP)."""
    z = (b - base) / amp
    z = np.clip(z, 1e-3, 1.0 - 1e-3)
    return d0 + (1.0 / k) * np.log(z / (1.0 - z))


def main():
    print("=" * 80)
    print("CL. Revealed set point from binge-rate inversion (preregistered)")
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
    df["binge"] = (df["calories"] > df["tdee"] + BINGE_THRESHOLD).astype(int)

    pre = df[df["effective_level"] == 0].copy().reset_index(drop=True)
    print(f"\nPre-tirz days: {len(pre)}, total binges: {pre['binge'].sum()}, "
          f"rate: {pre['binge'].mean()*100:.2f}%")

    # ── Calibration: fit f on BINNED (SP_dist, binge_rate) data with bounded k>0 ──
    fm_series = pre.set_index("date")["fm"]
    ema_sp_60 = fm_series.ewm(halflife=60, min_periods=30).mean()
    pre["ema_sp_60"] = pre["date"].map(ema_sp_60.to_dict())
    pre["d_ema"] = pre["ema_sp_60"] - pre["fm"]  # positive below SP

    # Bin per-day binge rate by SP-distance, then fit. Matches AG.
    fit_data = pre.dropna(subset=["d_ema", "binge"])
    bins = np.arange(-15, 16, 1.0)
    fit_data["d_bin"] = pd.cut(fit_data["d_ema"], bins=bins, include_lowest=True)
    bin_summary = (fit_data.groupby("d_bin", observed=True)
                   .agg(d_mid=("d_ema", "mean"), br=("binge", "mean"), n=("binge", "size"))
                   .dropna()
                   .reset_index(drop=True))
    bin_summary = bin_summary[bin_summary["n"] >= 30]
    print("\nBinned calibration (≥30 days per bin):")
    for _, r in bin_summary.iterrows():
        print(f"  d ≈ {r['d_mid']:+5.1f}  binge_rate = {r['br']*100:5.2f}%  n={int(r['n'])}")

    popt, _ = curve_fit(
        logistic, bin_summary["d_mid"], bin_summary["br"],
        p0=[0.03, 0.05, 0.2, 5.0],
        bounds=([0.005, 0.005, 0.05, -10.0], [0.10, 0.20, 2.0, 30.0]),
        maxfev=10000)
    base, amp, k, d0 = popt
    print(f"\nCalibrated f(d): B = {base:.4f} + {amp:.4f} / (1 + exp(-{k:.3f}·(d - {d0:.2f})))")
    print(f"  baseline {base*100:.2f}%  saturation {(base+amp)*100:.2f}%  "
          f"midpoint at d={d0:.2f} lb below SP, slope k={k:.3f}/lb")

    # ── Invert at two smoothing windows: 30d (fast), 90d (slow) ──
    for window in (30, 90):
        col = f"br_{window}c"
        pre[col] = pre["binge"].rolling(window, min_periods=window // 2, center=True).mean()
        d_col = f"d_revealed_{window}"
        sp_col = f"sp_revealed_{window}"
        pre[d_col] = pre[col].apply(
            lambda b: invert(b, base, amp, k, d0) if pd.notna(b) else np.nan)
        pre[sp_col] = pre["fm"] + pre[d_col]

    # ── Diagnostic correlations ──
    print("\n--- Correlations (where both defined) ---")
    for window in (30, 90):
        sp = pre[f"sp_revealed_{window}"]
        v = sp.notna() & pre["ema_sp_60"].notna()
        r_ema = float(np.corrcoef(sp[v], pre.loc[v, "ema_sp_60"])[0, 1])
        r_fm = float(np.corrcoef(sp[v], pre.loc[v, "fm"])[0, 1])
        print(f"  window={window}d  Corr(SP_revealed, EMA-SP-60d) = {r_ema:+.3f}, "
              f"Corr(SP_revealed, FM) = {r_fm:+.3f}")

    # AR control test: does SP_revealed gap survive lagged binge controls?
    pre["binge_lag7_30"] = pre["binge"].rolling(24, min_periods=12).sum().shift(7)
    pre["binge_lag30_90"] = pre["binge"].rolling(60, min_periods=30).sum().shift(30)

    for window in (30, 90):
        gap = pre[f"d_revealed_{window}"]
        # Predict next-7d binge frequency from gap, with AR controls
        future = pre["binge"].rolling(7, min_periods=4).sum().shift(-7)
        v = gap.notna() & future.notna() & pre["binge_lag7_30"].notna() & pre["binge_lag30_90"].notna()
        if v.sum() < 100:
            continue
        sub = pre.loc[v, [f"d_revealed_{window}", "binge_lag7_30", "binge_lag30_90"]].copy()
        sub["future"] = future[v]
        # Bivariate r and partial r
        r_marginal = float(np.corrcoef(sub[f"d_revealed_{window}"], sub["future"])[0, 1])
        # Partial r controlling for AR terms
        from numpy.linalg import lstsq
        X_ar = np.column_stack([np.ones(len(sub)), sub["binge_lag7_30"], sub["binge_lag30_90"]])
        gap_resid = sub[f"d_revealed_{window}"].values - X_ar @ lstsq(X_ar, sub[f"d_revealed_{window}"].values, rcond=None)[0]
        fut_resid = sub["future"].values - X_ar @ lstsq(X_ar, sub["future"].values, rcond=None)[0]
        r_partial = float(np.corrcoef(gap_resid, fut_resid)[0, 1])
        print(f"  window={window}d  Predicting NEXT-7d binge count from SP gap:")
        print(f"     marginal r = {r_marginal:+.3f}")
        print(f"     partial  r = {r_partial:+.3f}  (controlling for prior 7-30d, 30-90d binges)")

    # ── Plot ──
    fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    ax = axes[0]
    ax.plot(pre["date"], pre["fm"], "-", color="grey", alpha=0.4, label="FM (Kalman)")
    ax.plot(pre["date"], pre["ema_sp_60"], "-", color="tab:blue", lw=2, label="EMA SP (HL=60d)")
    ax.plot(pre["date"], pre["sp_revealed_90"], "-", color="tab:red", lw=2,
            label="Revealed SP (90d binge window)")
    ax.plot(pre["date"], pre["sp_revealed_30"], "-", color="tab:orange", lw=1, alpha=0.7,
            label="Revealed SP (30d binge window)")
    ax.set_ylabel("FM (lb)")
    ax.set_title("CL — Revealed SP from binge-rate inversion vs EMA SP")
    ax.legend(loc="upper left")
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(pre["date"], pre["d_revealed_90"], "-", color="tab:red", lw=2,
            label="Revealed SP - FM (90d window)")
    ax.plot(pre["date"], pre["d_revealed_30"], "-", color="tab:orange", lw=1, alpha=0.7,
            label="Revealed SP - FM (30d window)")
    ax.plot(pre["date"], pre["d_ema"], "-", color="tab:blue", lw=2, alpha=0.6,
            label="EMA SP - FM (HL=60d)")
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlabel("Date")
    ax.set_ylabel("SP gap (lb, positive = below SP)")
    ax.legend(loc="upper left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    plot_path = ROOT / "analysis" / "CL_revealed_sp.png"
    fig.savefig(plot_path, dpi=140)
    print(f"\nSaved {plot_path}")

    out = pre[["date", "fm", "ema_sp_60", "d_ema",
               "br_30c", "br_90c",
               "sp_revealed_30", "sp_revealed_90",
               "d_revealed_30", "d_revealed_90"]]
    out_csv = ROOT / "analysis" / "CL_revealed_sp.csv"
    out.to_csv(out_csv, index=False)
    print(f"Saved {out_csv}")


if __name__ == "__main__":
    main()
