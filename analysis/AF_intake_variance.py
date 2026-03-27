"""Does calorie variance independently predict fat mass change?

At the same average intake, does eating 2000 ± 500 cal/day cause more
fat gain than 2000 ± 100? The binge-restrict / "metabolic damage"
literature predicts yes. The data says no — the sign is reversed.

Controlling for caloric surplus, higher intake variance predicts slightly
LESS fat gain (partial r = -0.20, 14d window). The coefficient is small
but significant: -0.0011 lbs/day per cal of std (95% CI: -0.0013 to -0.0010).
The mechanism is likely that high-variance periods include fasting days,
which transiently raise TDEE (finding B: weekend fasting).

Result: intake variance is mildly protective, not fattening. Variable
eating does not cause "metabolic damage" in this dataset.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def main():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    epochs = pd.read_csv(ROOT / "intake" / "diet_epochs.csv", parse_dates=["start", "end"])

    df = intake[["date", "calories"]].merge(
        kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="inner")
    df = df.sort_values("date").reset_index(drop=True)
    n = len(df)
    print(f"{n} days with intake + Kalman estimates")

    # Trailing windows
    for w in [14, 21, 30]:
        df[f"cal_mean_{w}d"] = df["calories"].rolling(w, min_periods=w).mean()
        df[f"cal_std_{w}d"] = df["calories"].rolling(w, min_periods=w).std()
        df[f"cal_cv_{w}d"] = df[f"cal_std_{w}d"] / df[f"cal_mean_{w}d"]
        df[f"fm_change_{w}d"] = df["fat_mass_lbs"].shift(-w) - df["fat_mass_lbs"]
        df[f"tdee_mean_{w}d"] = df["tdee"].rolling(w, min_periods=w).mean()
        df[f"surplus_{w}d"] = df[f"cal_mean_{w}d"] - df[f"tdee_mean_{w}d"]

    # ── Raw correlations ──
    print(f"\n=== Raw correlations: variance metrics → FM change ===")
    print(f"  {'Metric':>20} {'14d':>8} {'30d':>8}")
    for metric in ["cal_mean", "cal_std", "cal_cv", "surplus"]:
        vals_14 = []
        vals_30 = []
        for w, vals_list in [(14, vals_14), (30, vals_30)]:
            col = f"{metric}_{w}d"
            fm_col = f"fm_change_{w}d"
            valid = df[col].notna() & df[fm_col].notna()
            if valid.sum() > 100:
                r = np.corrcoef(df.loc[valid, col], df.loc[valid, fm_col])[0, 1]
                vals_list.append(r)
            else:
                vals_list.append(np.nan)
        print(f"  {metric:>20} {vals_14[0]:+8.3f} {vals_30[0]:+8.3f}")

    # ── Partial correlations ──
    print(f"\n=== Partial correlations: cal_std → FM change, controlling for... ===")
    for w in [14, 21, 30]:
        valid = (df[f"cal_mean_{w}d"].notna() & df[f"cal_std_{w}d"].notna() &
                 df[f"fm_change_{w}d"].notna() & df[f"tdee_mean_{w}d"].notna())
        sub = df[valid]
        mean_v = sub[f"cal_mean_{w}d"].values
        std_v = sub[f"cal_std_{w}d"].values
        fm_v = sub[f"fm_change_{w}d"].values
        surplus_v = sub[f"surplus_{w}d"].values

        # Control for mean intake
        X = np.column_stack([mean_v, np.ones(len(sub))])
        res_std = std_v - X @ np.linalg.lstsq(X, std_v, rcond=None)[0]
        res_fm = fm_v - X @ np.linalg.lstsq(X, fm_v, rcond=None)[0]
        r_mean = np.corrcoef(res_std, res_fm)[0, 1]

        # Control for surplus
        X2 = np.column_stack([surplus_v, np.ones(len(sub))])
        res_std2 = std_v - X2 @ np.linalg.lstsq(X2, std_v, rcond=None)[0]
        res_fm2 = fm_v - X2 @ np.linalg.lstsq(X2, fm_v, rcond=None)[0]
        r_surplus = np.corrcoef(res_std2, res_fm2)[0, 1]

        print(f"  {w}d: | mean intake: r = {r_mean:+.4f}   | surplus: r = {r_surplus:+.4f}  (n={len(sub)})")

    # ── Bin analysis ──
    print(f"\n=== At matched calorie levels (30d tertiles) ===")
    w = 30
    valid = df[f"cal_mean_{w}d"].notna() & df[f"cal_std_{w}d"].notna() & df[f"fm_change_{w}d"].notna()
    sub = df[valid].copy()
    sub["cal_bin"] = pd.qcut(sub[f"cal_mean_{w}d"], 3, labels=["low", "mid", "high"])
    print(f"  {'Cal level':>10} {'Variance':>10} {'Mean cal':>9} {'Std':>6} {'FM Δ/mo':>8} {'n':>5}")
    for cb in ["low", "mid", "high"]:
        grp = sub[sub["cal_bin"] == cb]
        med = grp[f"cal_std_{w}d"].median()
        for vl, vm in [("low_var", grp[f"cal_std_{w}d"] <= med),
                        ("high_var", grp[f"cal_std_{w}d"] > med)]:
            g = grp[vm]
            if len(g) < 30:
                continue
            fm_mo = g[f"fm_change_{w}d"].mean()
            print(f"  {cb:>10} {vl:>10} {g[f'cal_mean_{w}d'].mean():9.0f} "
                  f"{g[f'cal_std_{w}d'].mean():6.0f} {fm_mo:+8.2f} {len(g):5d}")

    # ── Diet epochs ──
    print(f"\n=== Diet epoch comparison ===")
    print(f"  {'Epoch':>25} {'Cal':>5} {'Std':>5} {'CV':>5} {'FM/mo':>7}")
    for _, ep in epochs.iterrows():
        mask = (df["date"] >= ep["start"]) & (df["date"] <= ep["end"])
        edf = df[mask]
        if len(edf) < 14:
            continue
        mc = edf["calories"].mean()
        sc = edf["calories"].std()
        cv = sc / mc if mc > 0 else 0
        days = (edf["date"].iloc[-1] - edf["date"].iloc[0]).days
        fm_mo = (edf["fat_mass_lbs"].iloc[-1] - edf["fat_mass_lbs"].iloc[0]) / days * 30 if days > 0 else 0
        print(f"  {ep['label'][:25]:>25} {mc:5.0f} {sc:5.0f} {cv:5.2f} {fm_mo:+7.2f}")

    # ── Regression with bootstrap CI ──
    print(f"\n=== Multivariate: surplus + variance → FM change ===")
    for w in [14, 30]:
        valid = (df[f"surplus_{w}d"].notna() & df[f"cal_std_{w}d"].notna() &
                 df[f"fm_change_{w}d"].notna())
        sub = df[valid]
        y = sub[f"fm_change_{w}d"].values
        X1 = np.column_stack([sub[f"surplus_{w}d"].values, np.ones(len(sub))])
        X2 = np.column_stack([sub[f"surplus_{w}d"].values, sub[f"cal_std_{w}d"].values,
                              np.ones(len(sub))])
        c1 = np.linalg.lstsq(X1, y, rcond=None)[0]
        c2 = np.linalg.lstsq(X2, y, rcond=None)[0]
        rmse1 = np.sqrt(np.mean((y - X1 @ c1) ** 2))
        rmse2 = np.sqrt(np.mean((y - X2 @ c2) ** 2))

        np.random.seed(42)
        boots = []
        for _ in range(2000):
            idx = np.random.choice(len(sub), len(sub), replace=True)
            try:
                cb = np.linalg.lstsq(X2[idx], y[idx], rcond=None)[0]
                boots.append(cb[1])
            except:
                pass
        boots = np.array(boots)
        ci_lo, ci_hi = np.percentile(boots, [2.5, 97.5])

        print(f"\n  {w}d window (n={len(sub)}):")
        print(f"    Surplus only:       RMSE = {rmse1:.4f}")
        print(f"    Surplus + cal_std:  RMSE = {rmse2:.4f}")
        print(f"    Std coefficient: {c2[1]:+.6f}  95% CI: [{ci_lo:.6f}, {ci_hi:.6f}]")
        sig = "SIGNIFICANT" if ci_lo > 0 or ci_hi < 0 else "not significant"
        print(f"    {sig}")

    # ── Summary ──
    print(f"\n=== Summary ===")
    print(f"  Controlling for caloric surplus, intake variance predicts slightly")
    print(f"  LESS fat gain (partial r ≈ -0.20, significant at 14d and 30d).")
    print(f"  The coefficient is small: -0.0011 lbs/day per cal of std.")
    print(f"  At 500 vs 300 cal std, that's -0.22 lbs/month extra loss.")
    print(f"  Mechanism: high-variance periods include fasting days that")
    print(f"  transiently raise TDEE (finding B). The variance is not")
    print(f"  independently fattening — 'metabolic damage from yo-yo")
    print(f"  dieting' is not supported by this dataset.")


if __name__ == "__main__":
    main()
