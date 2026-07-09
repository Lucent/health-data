"""CK. Reveal a hidden set point from a storage-efficiency regime change.

# Preregistration (committed before running the model)

## Construction

For each day t in the pre-tirz era:

    baseline_TDEE_t  = expected_RMR_t × 1.2          (sedentary reference, NOT Kalman TDEE)
    naive_surplus_t  = intake_t − baseline_TDEE_t    (canonical 1.09× already in intake)
    cal_per_lb_FM_t  = 3500 + 1100 × (1 − ff)/ff     (Forbes correction)
    predicted_ΔFM_t  = naive_surplus_t / cal_per_lb_FM_t
    actual_ΔFM_t     = Kalman FM trajectory diff
    clip_signal_t    = predicted_ΔFM_t − actual_ΔFM_t   (>0 = calories disappeared)

Aggregated to 7-day rolling sums (to span the metabolic ~15d HL).

For each 90-day pre-tirz window (sliding, step 30 days), fit:

    clip_signal = α + γ · surplus + β · surplus · max(0, FM − SP*)

with (α, γ, β) linear given SP*; SP* found by golden-section search across the
window's observed FM range. Note: SP* is identified PURELY from a metabolic
regime change. No EMA, no smoothing of FM. The within-window threshold cannot
inherit the BR autocorrelation because it doesn't reference FM history.

## What we expect if set-point-mediated clipping is real

1. β > 0 in most windows (target: > 70% of windows have β > 0).
2. Median SP* lies in 25-90 lb FM (subject's lifetime range), not at the bounds
   of each window's FM span.
3. SP*(t) shows phase structure aligned with FM history:
   - low (~15-25 lb) during the 2013-14 willpower bottom;
   - rising through 2014-2018 regain to ~50-80 lb;
   - tracking FM in 2018-2024 maintenance/regain at the FM band.
4. Median ΔR² (full vs no-SP null) > 0.05.
5. SP*(t) correlates with EMA(FM, 50d) at r ∈ [0.4, 0.9] — i.e. positive but
   not a perfect smoothing. A r > 0.95 means we just rediscovered EMA from a
   different angle (still useful but no decircularisation). r < 0.4 means the
   two are tracking different things and revealed-SP is the new ground truth.
6. β magnitude: at a +5 lb above-SP gap and surplus = 500 cal/day (typical
   binge-period day), the clipping should be on the order of tens to hundreds
   of cal/day, consistent with AI's 612 cal residual at the deepest above-SP
   bin. In window units (7d cumulative), expect β · 5 · 3500 ≈ 100-500 cal of
   clipping per 7-day window per lb above SP per cal of surplus — translating
   to β on the order of 1e-5 to 1e-4 (lb / cal²).

## What we expect if set-point-mediated clipping is absent

- β CIs span zero; sign flips between windows; fraction with β > 0 near 0.5.
- SP* often hits the window's FM range bounds (no real interior threshold).
- Median ΔR² ≈ 0.
- Then the SP-as-threshold framing is retired; the AM continuous-pressure
  reading takes precedence; CK becomes a registered negative result.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parent.parent

LBS_TO_KG = 0.45359237
FORBES_C_KG = 10.4

def cal_per_lb_FM(fm_lbs):
    fm_kg = max(float(fm_lbs) * LBS_TO_KG, 5.0)
    ff = 1.0 - (FORBES_C_KG ** 2) / (FORBES_C_KG + fm_kg) ** 2
    if ff <= 1e-3:
        return 4500.0
    return 3500.0 + 1100.0 * (1.0 - ff) / ff


def fit_window(surplus, fm, clip):
    """Fit clip = α + γ·surplus + β·surplus·max(0, FM - SP*) on a window.
    Returns dict with sp_star, coefs, R²s.
    """
    fm_lo, fm_hi = float(np.min(fm)), float(np.max(fm))
    if fm_hi - fm_lo < 1.0:
        return None  # too little FM variation to identify SP*

    def ss_at(sp_star):
        gate = np.maximum(0.0, fm - sp_star)
        X = np.column_stack([np.ones_like(surplus), surplus, surplus * gate])
        try:
            coefs, *_ = np.linalg.lstsq(X, clip, rcond=None)
        except np.linalg.LinAlgError:
            return 1e18
        pred = X @ coefs
        return float(np.sum((clip - pred) ** 2))

    res = minimize_scalar(ss_at, bounds=(fm_lo + 0.5, fm_hi - 0.5), method="bounded",
                          options={"xatol": 0.1})
    sp_star = float(res.x)

    gate = np.maximum(0.0, fm - sp_star)
    X = np.column_stack([np.ones_like(surplus), surplus, surplus * gate])
    coefs, *_ = np.linalg.lstsq(X, clip, rcond=None)
    pred = X @ coefs
    ss_res = float(np.sum((clip - pred) ** 2))
    ss_tot = float(np.sum((clip - clip.mean()) ** 2))
    r2_full = 1 - ss_res / ss_tot if ss_tot > 1e-10 else 0.0

    X0 = np.column_stack([np.ones_like(surplus), surplus])
    coefs0, *_ = np.linalg.lstsq(X0, clip, rcond=None)
    pred0 = X0 @ coefs0
    ss_res0 = float(np.sum((clip - pred0) ** 2))
    r2_null = 1 - ss_res0 / ss_tot if ss_tot > 1e-10 else 0.0

    return {
        "sp_star": sp_star,
        "alpha": float(coefs[0]),
        "gamma": float(coefs[1]),
        "beta": float(coefs[2]),
        "r2_full": r2_full,
        "r2_null": r2_null,
        "delta_r2": r2_full - r2_null,
        "n": len(surplus),
        "fm_lo": fm_lo,
        "fm_hi": fm_hi,
    }


def main():
    print("=" * 80)
    print("CK. Revealed set point from clip-signal regime fit (preregistered)")
    print("=" * 80)

    kf = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    INTAKE_MULTIPLIER = 1.09
    intake = intake.copy()
    intake["calories"] = intake["calories"] * INTAKE_MULTIPLIER

    df = kf[["date", "fat_mass_lbs"]].rename(columns={"fat_mass_lbs": "fm"})
    df = df.merge(intake[["date", "calories"]], on="date", how="left")
    df = df.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df = df.sort_values("date").reset_index(drop=True)
    df = df.dropna(subset=["fm", "calories", "expected_rmr"])

    # Keep only pre-tirz days for the SP fit
    pre = df[df["effective_level"] == 0].copy().reset_index(drop=True)

    pre["baseline_tdee"] = pre["expected_rmr"] * 1.2
    pre["naive_surplus"] = pre["calories"] - pre["baseline_tdee"]
    pre["cal_per_lb"] = pre["fm"].apply(cal_per_lb_FM)
    pre["predicted_dfm_daily"] = pre["naive_surplus"] / pre["cal_per_lb"]

    WINDOW_AGG = 7
    pre["surplus_7d"] = pre["naive_surplus"].rolling(WINDOW_AGG, min_periods=WINDOW_AGG).sum()
    pre["fm_change_7d"] = pre["fm"].diff(WINDOW_AGG)
    pre["predicted_dfm_7d"] = pre["surplus_7d"] / pre["cal_per_lb"]
    pre["clip_signal_7d"] = pre["predicted_dfm_7d"] - pre["fm_change_7d"]

    valid = pre.dropna(subset=["surplus_7d", "fm_change_7d"]).reset_index(drop=True)
    print(f"\nPre-tirz, post-roll valid days: {len(valid)}")
    print(f"Surplus range (7d cum): {valid['surplus_7d'].min():.0f} to "
          f"{valid['surplus_7d'].max():.0f} cal")
    print(f"FM range: {valid['fm'].min():.1f} to {valid['fm'].max():.1f} lb")

    # Sliding 90-day windows, step 30 days
    WINDOW_FIT = 90
    STEP = 30
    rows = []
    for s in range(0, len(valid) - WINDOW_FIT, STEP):
        chunk = valid.iloc[s:s + WINDOW_FIT]
        if len(chunk) < 60:
            continue
        fit = fit_window(
            chunk["surplus_7d"].values,
            chunk["fm"].values,
            chunk["clip_signal_7d"].values,
        )
        if fit is None:
            continue
        fit["window_start"] = chunk["date"].iloc[0]
        fit["window_mid"] = chunk["date"].iloc[len(chunk) // 2]
        fit["window_end"] = chunk["date"].iloc[-1]
        fit["fm_med"] = float(chunk["fm"].median())
        rows.append(fit)

    res = pd.DataFrame(rows)
    out_csv = ROOT / "analysis" / "CK_revealed_sp.csv"
    res.to_csv(out_csv, index=False)
    print(f"\nFit {len(res)} sliding 90-day windows. Saved {out_csv}")

    # ─── headline preregistered checks ───
    print("\n--- Preregistered checks ---")
    n = len(res)
    n_beta_pos = int((res["beta"] > 0).sum())
    print(f"  Windows with β > 0:           {n_beta_pos}/{n}  ({100*n_beta_pos/n:.0f}%)  "
          f"[predicted > 70%]")
    print(f"  Median SP*:                   {res['sp_star'].median():.1f} lb FM  "
          f"[predicted 25-90 lb]")
    at_bound = ((res["sp_star"] - res["fm_lo"]).abs() < 1.0).sum() + \
               ((res["fm_hi"] - res["sp_star"]).abs() < 1.0).sum()
    print(f"  Windows with SP* at bound:    {at_bound}/{n}  ({100*at_bound/n:.0f}%)  "
          f"[predicted few]")
    print(f"  Median ΔR² (full vs null):    {res['delta_r2'].median():+.4f}  "
          f"[predicted > 0.05]")
    print(f"  Median β:                     {res['beta'].median():+.2e}  "
          f"[predicted ~1e-5 to 1e-4]")

    # SP* vs EMA SP correlation (on fit window mids)
    fm_series = pre.set_index("date")["fm"]
    ema_sp = fm_series.ewm(halflife=50, min_periods=30).mean()
    res["ema_sp_at_mid"] = pd.to_datetime(res["window_mid"]).map(ema_sp.to_dict())
    valid_corr = res.dropna(subset=["ema_sp_at_mid"])
    if len(valid_corr) > 5:
        r_spema = float(np.corrcoef(valid_corr["sp_star"], valid_corr["ema_sp_at_mid"])[0, 1])
        print(f"  Corr(SP*, EMA-SP-50d):        r = {r_spema:+.3f}  "
              f"[predicted +0.4 to +0.9]")

    # Plot — single panel: FM trajectory, EMA SP, revealed SP*
    fig, ax = plt.subplots(figsize=(13, 5))
    ax.plot(fm_series.index, fm_series.values, "-", color="grey", alpha=0.4,
            label="FM (Kalman, pre-tirz)")
    ax.plot(ema_sp.index, ema_sp.values, "-", color="tab:blue", lw=2,
            label="EMA SP (HL=50d)")
    sp_x = pd.to_datetime(res["window_mid"])
    ax.plot(sp_x, res["sp_star"], "o-", color="tab:red", lw=2,
            label="Revealed SP* (90d windows)")
    # Bound markers where SP* hit the window range
    bound_mask = (((res["sp_star"] - res["fm_lo"]).abs() < 1.0) |
                  ((res["fm_hi"] - res["sp_star"]).abs() < 1.0))
    ax.scatter(sp_x[bound_mask], res["sp_star"][bound_mask], marker="x",
               color="black", s=70, zorder=5, label="SP* at window bound")

    ax.set_xlabel("Date")
    ax.set_ylabel("FM (lb)")
    ax.set_title("CK — Revealed SP* from clip-signal regime fit  vs  EMA SP")
    ax.legend(loc="upper left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    plot_path = ROOT / "analysis" / "CK_revealed_sp.png"
    fig.savefig(plot_path, dpi=140)
    print(f"\nSaved {plot_path}")

    # Diagnostic detail per window
    print("\n--- Per-window detail (first/middle/last) ---")
    print(f"  {'window_mid':<12} {'SP*':>6} {'β':>10} {'ΔR²':>8} {'β>0?':>6}")
    show = list(range(0, len(res), max(1, len(res) // 10)))
    for i in show:
        r = res.iloc[i]
        print(f"  {pd.to_datetime(r['window_mid']).strftime('%Y-%m-%d'):<12} "
              f"{r['sp_star']:>6.1f} {r['beta']:>+10.2e} {r['delta_r2']:>+8.3f} "
              f"{'yes' if r['beta'] > 0 else 'no':>6}")


if __name__ == "__main__":
    main()
