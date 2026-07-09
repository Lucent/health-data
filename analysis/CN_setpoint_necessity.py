"""CN. Two diagnostics: is a defended set point necessary at all?

Test 1 — Intake autoregression vs SP residual.
  Predict 7-day mean intake at week t from weekly intake history (lags 1,2,4,8 weeks)
  plus current FM and sp_distance. If sp_distance retains predictive power after AR
  controls, the set-point framework adds something. If it collapses, intake is just
  habit + autocorrelation and SP is unnecessary.

Test 2 — Long-term attractor test.
  Find every sliding 90-day pre-tirz window where the subject was approximately at
  maintenance (|mean daily surplus| < 50 cal/day, where surplus = intake − Kalman TDEE).
  Plot the FM values during these windows. If they cluster around a single value,
  there's a defended weight. If they span the subject's lifetime FM range, there is
  no attractor — FM sits wherever it happens to.

The first answers "does SP carry information beyond intake habit?" The second answers
"does the body return to a single value when behaviorally neutral?" Both are needed.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.linalg import lstsq

ROOT = Path(__file__).resolve().parent.parent
INTAKE_MULTIPLIER = 1.09


def partial_corr(y, x_target, x_controls):
    Xc = np.column_stack([np.ones(len(y))] + list(x_controls))
    by, *_ = lstsq(Xc, y, rcond=None)
    y_r = y - Xc @ by
    bx, *_ = lstsq(Xc, x_target, rcond=None)
    x_r = x_target - Xc @ bx
    return float(np.corrcoef(y_r, x_r)[0, 1])


def main():
    print("=" * 80)
    print("CN. Is a defended set point necessary?")
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

    pre = df[df["effective_level"] == 0].copy().reset_index(drop=True)
    pre["surplus"] = pre["calories"] - pre["tdee"]

    # Weekly aggregates
    pre["intake_7d"] = pre["calories"].rolling(7, min_periods=7).mean()
    pre["fm_7d"] = pre["fm"].rolling(7, min_periods=7).mean()
    pre["surplus_7d"] = pre["surplus"].rolling(7, min_periods=7).mean()

    # SP distance at HL=60d
    fm_s = pre.set_index("date")["fm"]
    ema60 = fm_s.ewm(halflife=60, min_periods=30).mean()
    pre["sp_dist_60"] = pre["date"].map(ema60.to_dict()) - pre["fm"]

    # Lagged intake (in days; 7, 14, 28, 56)
    for lag in (7, 14, 28, 56):
        pre[f"intake_lag_{lag}"] = pre["intake_7d"].shift(lag)

    valid = pre.dropna(subset=[
        "intake_7d", "intake_lag_7", "intake_lag_14",
        "intake_lag_28", "intake_lag_56",
        "fm", "sp_dist_60"]).copy()
    print(f"\nValid pre-tirz days: {len(valid)}")

    # ================================================================
    # Test 1 — Intake AR vs SP residual
    # ================================================================
    print("\n" + "=" * 80)
    print("TEST 1 — Predicting 7-day mean intake from AR + FM + sp_distance")
    print("=" * 80)

    y = valid["intake_7d"].values
    ar1 = valid["intake_lag_7"].values
    ar2 = valid["intake_lag_14"].values
    ar3 = valid["intake_lag_28"].values
    ar4 = valid["intake_lag_56"].values
    fm_v = valid["fm"].values
    sp = valid["sp_dist_60"].values

    # Models, in order of what they include
    def r2(y_obs, X):
        b, *_ = lstsq(X, y_obs, rcond=None)
        pred = X @ b
        ss_res = np.sum((y_obs - pred) ** 2)
        ss_tot = np.sum((y_obs - y_obs.mean()) ** 2)
        return 1 - ss_res / ss_tot

    intercept = np.ones(len(y))
    X_const  = np.column_stack([intercept])
    X_ar1    = np.column_stack([intercept, ar1])
    X_ar4    = np.column_stack([intercept, ar1, ar2, ar3, ar4])
    X_ar_fm  = np.column_stack([intercept, ar1, ar2, ar3, ar4, fm_v])
    X_ar_full = np.column_stack([intercept, ar1, ar2, ar3, ar4, fm_v, sp])

    print(f"  Constant only           R² = {r2(y, X_const):.4f}")
    print(f"  + intake_{{t-7d}}         R² = {r2(y, X_ar1):.4f}  (AR1)")
    print(f"  + intake AR(7,14,28,56) R² = {r2(y, X_ar4):.4f}  (full AR)")
    print(f"  + FM                    R² = {r2(y, X_ar_fm):.4f}")
    print(f"  + sp_distance_60        R² = {r2(y, X_ar_full):.4f}")

    delta_sp = r2(y, X_ar_full) - r2(y, X_ar_fm)
    print(f"\n  ΔR² from adding sp_distance: {delta_sp:+.4f}")

    r_sp_partial = partial_corr(y, sp, [ar1, ar2, ar3, ar4, fm_v])
    r_sp_marginal = float(np.corrcoef(y, sp)[0, 1])
    print(f"  r(intake, sp_dist)  marginal:                    {r_sp_marginal:+.4f}")
    print(f"  r(intake, sp_dist)  partial (|AR + FM):          {r_sp_partial:+.4f}")

    if abs(r_sp_partial) < 0.10:
        verdict_1 = "SP carries no information beyond intake AR. Framework not needed."
    elif abs(r_sp_partial) < 0.20:
        verdict_1 = "SP carries marginal residual information after AR."
    else:
        verdict_1 = "SP retains substantial information after AR. Genuine signal."
    print(f"  → {verdict_1}")

    # ================================================================
    # Test 2 — Long-term attractor test
    # ================================================================
    print("\n" + "=" * 80)
    print("TEST 2 — Where does FM sit during maintenance?")
    print("=" * 80)
    WINDOW = 90
    THRESHOLD = 50
    pre["surplus_90d_mean"] = pre["surplus"].rolling(WINDOW, min_periods=WINDOW).mean()
    pre["fm_90d_mean"] = pre["fm"].rolling(WINDOW, min_periods=WINDOW).mean()

    maint = pre[(pre["surplus_90d_mean"].abs() < THRESHOLD) &
                pre["surplus_90d_mean"].notna()].copy()

    # Deduplicate: take one row per ~30-day stride to reduce overlap
    if len(maint):
        maint = maint.iloc[::30].copy()

    print(f"  Maintenance windows (|mean surplus| < {THRESHOLD} cal/day): {len(maint)}")
    print(f"  FM during maintenance: min={maint['fm_90d_mean'].min():.1f}, "
          f"max={maint['fm_90d_mean'].max():.1f}, "
          f"median={maint['fm_90d_mean'].median():.1f}")

    print("\n  Per-decade maintenance FM:")
    for label, lo, hi in [("2011-2014", "2011-01-01", "2014-12-31"),
                          ("2015-2019", "2015-01-01", "2019-12-31"),
                          ("2020-2024", "2020-01-01", "2024-12-31")]:
        m = maint[(maint["date"] >= lo) & (maint["date"] <= hi)]
        if len(m):
            print(f"    {label}: n={len(m):>2}  FM range "
                  f"{m['fm_90d_mean'].min():.1f}-{m['fm_90d_mean'].max():.1f} lb, "
                  f"median {m['fm_90d_mean'].median():.1f}")

    fm_iqr = maint["fm_90d_mean"].quantile(0.75) - maint["fm_90d_mean"].quantile(0.25)
    fm_full_range = pre["fm"].max() - pre["fm"].min()
    spread_frac = fm_iqr / fm_full_range
    print(f"\n  Maintenance-FM IQR: {fm_iqr:.1f} lb")
    print(f"  Lifetime FM range:  {fm_full_range:.1f} lb")
    print(f"  Maintenance IQR / lifetime range: {spread_frac:.2%}")

    if spread_frac < 0.15:
        verdict_2 = "Maintenance FM clusters tightly. Defended weight likely exists."
    elif spread_frac < 0.35:
        verdict_2 = "Maintenance FM somewhat clustered. Soft attractor possible."
    else:
        verdict_2 = "Maintenance FM is spread across the lifetime range. NO single attractor."
    print(f"  → {verdict_2}")

    # ── Plots ──
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    ax = axes[0, 0]
    obs = valid["intake_7d"].values
    pred_ar = X_ar4 @ lstsq(X_ar4, y, rcond=None)[0]
    ax.scatter(pred_ar, obs, s=2, alpha=0.4)
    lo, hi = obs.min(), obs.max()
    ax.plot([lo, hi], [lo, hi], "k--", alpha=0.5)
    ax.set_xlabel("Predicted intake (AR only)")
    ax.set_ylabel("Observed 7-day intake")
    ax.set_title(f"AR(7,14,28,56) prediction  R² = {r2(y, X_ar4):.3f}")

    ax = axes[0, 1]
    pred_full = X_ar_full @ lstsq(X_ar_full, y, rcond=None)[0]
    ax.scatter(pred_full, obs, s=2, alpha=0.4)
    ax.plot([lo, hi], [lo, hi], "k--", alpha=0.5)
    ax.set_xlabel("Predicted intake (AR + FM + sp_dist)")
    ax.set_ylabel("Observed")
    ax.set_title(f"AR + FM + sp_dist  R² = {r2(y, X_ar_full):.3f}  (Δ {delta_sp:+.4f})")

    ax = axes[1, 0]
    ax.plot(pre["date"], pre["fm"], color="grey", alpha=0.4, label="FM")
    ax.scatter(maint["date"], maint["fm_90d_mean"], s=40, color="tab:red",
               zorder=5, label="Maintenance windows")
    ax.set_xlabel("Date")
    ax.set_ylabel("FM (lb)")
    ax.set_title("FM trajectory with maintenance points highlighted")
    ax.legend()

    ax = axes[1, 1]
    ax.hist(maint["fm_90d_mean"], bins=20, color="tab:red", edgecolor="black")
    ax.set_xlabel("Mean FM during maintenance window (lb)")
    ax.set_ylabel("count")
    ax.set_title(f"Distribution of FM at maintenance  IQR={fm_iqr:.1f} / range={fm_full_range:.1f} lb")

    fig.tight_layout()
    plot_path = ROOT / "analysis" / "CN_setpoint_necessity.png"
    fig.savefig(plot_path, dpi=140)
    print(f"\nSaved {plot_path}")

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"  Test 1: {verdict_1}")
    print(f"  Test 2: {verdict_2}")


if __name__ == "__main__":
    main()
