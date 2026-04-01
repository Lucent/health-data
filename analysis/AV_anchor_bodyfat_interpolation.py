#!/usr/bin/env python3
"""AV. Daily FM/FFM interpolation from weight + composition anchors only.

Build a daily body-composition path without using intake, TDEE, or RMR:
  - On weigh-in days, corrected weight constrains FM + FFM.
  - On composition-scan days, FM and FFM are anchored directly.
  - Between anchors, both series are regularized to be smooth, with FFM
    constrained to vary more slowly than FM.

This is intended as a less endogenous alternative to the intake-informed
Kalman pipeline when testing set-point style appetite hypotheses.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import lsqr

ROOT = Path(__file__).resolve().parent.parent

# Observation weights
W_WEIGHT = 4.0
W_COMP_FM = 18.0
W_COMP_FFM = 18.0

# Smoothness penalties
L2_FM_SECOND = 0.20
L2_FFM_SECOND = 1.25
L2_FM_FIRST = 0.015
L2_FFM_FIRST = 0.080

# Weak composition prior to keep trajectories plausible between long gaps.
L2_PRIOR_FM = 0.05
L2_PRIOR_FFM = 0.10


def load_data():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    weight = pd.read_csv(ROOT / "analysis" / "P1_smoothed_weight.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "composition" / "composition.csv", parse_dates=["date"])
    return intake, weight, comp


def interp_prior(daily_dates, comp_dates, values):
    x = np.array([(d - daily_dates[0]) / np.timedelta64(1, "D") for d in daily_dates], dtype=float)
    xc = np.array([(d - daily_dates[0]) / np.timedelta64(1, "D") for d in comp_dates], dtype=float)
    return np.interp(x, xc, values)


def add_row_triplets(rows, cols, data, r, cvals):
    for c, v in cvals:
        rows.append(r)
        cols.append(c)
        data.append(v)


def main():
    intake, weight, comp = load_data()
    daily = pd.DataFrame({"date": pd.date_range(intake["date"].min(), intake["date"].max(), freq="D")})
    daily = daily.merge(weight[["date", "smoothed_weight_lbs"]], on="date", how="left")

    comp = comp[["date", "fat_mass_lbs", "lean_mass_lbs", "weight_lbs"]].dropna(
        subset=["fat_mass_lbs", "lean_mass_lbs"]
    ).copy()
    comp = comp.sort_values("date").drop_duplicates("date", keep="last").reset_index(drop=True)

    n = len(daily)
    fm_offset = 0
    ffm_offset = n

    date_to_idx = {d: i for i, d in enumerate(daily["date"].values)}
    weight_idx = []
    weight_vals = []
    for _, row in daily.dropna(subset=["smoothed_weight_lbs"]).iterrows():
        idx = date_to_idx[row["date"].to_datetime64()]
        weight_idx.append(idx)
        weight_vals.append(row["smoothed_weight_lbs"])

    comp_idx = []
    comp_fm = []
    comp_ffm = []
    for _, row in comp.iterrows():
        idx = date_to_idx.get(row["date"].to_datetime64())
        if idx is not None:
            comp_idx.append(idx)
            comp_fm.append(row["fat_mass_lbs"])
            comp_ffm.append(row["lean_mass_lbs"])

    comp_dates = np.array([daily["date"].values[i] for i in comp_idx])
    fm_prior = interp_prior(daily["date"].values, comp_dates, np.array(comp_fm))
    ffm_prior = interp_prior(daily["date"].values, comp_dates, np.array(comp_ffm))

    rows = []
    cols = []
    data = []
    b = []
    r = 0

    # Weight observations: FM + FFM = corrected scale weight.
    for idx, wt in zip(weight_idx, weight_vals):
        s = np.sqrt(W_WEIGHT)
        add_row_triplets(rows, cols, data, r, [(fm_offset + idx, s), (ffm_offset + idx, s)])
        b.append(s * wt)
        r += 1

    # Composition anchors.
    for idx, fm in zip(comp_idx, comp_fm):
        s = np.sqrt(W_COMP_FM)
        add_row_triplets(rows, cols, data, r, [(fm_offset + idx, s)])
        b.append(s * fm)
        r += 1
    for idx, ffm in zip(comp_idx, comp_ffm):
        s = np.sqrt(W_COMP_FFM)
        add_row_triplets(rows, cols, data, r, [(ffm_offset + idx, s)])
        b.append(s * ffm)
        r += 1

    # Weak prior to composition-only linear interpolation.
    for i in range(n):
        s_fm = np.sqrt(L2_PRIOR_FM)
        s_ffm = np.sqrt(L2_PRIOR_FFM)
        add_row_triplets(rows, cols, data, r, [(fm_offset + i, s_fm)])
        b.append(s_fm * fm_prior[i])
        r += 1
        add_row_triplets(rows, cols, data, r, [(ffm_offset + i, s_ffm)])
        b.append(s_ffm * ffm_prior[i])
        r += 1

    # First-difference penalties.
    for i in range(1, n):
        s_fm = np.sqrt(L2_FM_FIRST)
        s_ffm = np.sqrt(L2_FFM_FIRST)
        add_row_triplets(rows, cols, data, r, [(fm_offset + i - 1, -s_fm), (fm_offset + i, s_fm)])
        b.append(0.0)
        r += 1
        add_row_triplets(rows, cols, data, r, [(ffm_offset + i - 1, -s_ffm), (ffm_offset + i, s_ffm)])
        b.append(0.0)
        r += 1

    # Second-difference smoothness penalties.
    for i in range(1, n - 1):
        s_fm = np.sqrt(L2_FM_SECOND)
        s_ffm = np.sqrt(L2_FFM_SECOND)
        add_row_triplets(
            rows, cols, data, r,
            [(fm_offset + i - 1, s_fm), (fm_offset + i, -2 * s_fm), (fm_offset + i + 1, s_fm)]
        )
        b.append(0.0)
        r += 1
        add_row_triplets(
            rows, cols, data, r,
            [(ffm_offset + i - 1, s_ffm), (ffm_offset + i, -2 * s_ffm), (ffm_offset + i + 1, s_ffm)]
        )
        b.append(0.0)
        r += 1

    A = sparse.csr_matrix((data, (rows, cols)), shape=(r, 2 * n))
    b = np.asarray(b)

    sol = lsqr(A, b, atol=1e-8, btol=1e-8, iter_lim=20000)
    x = sol[0]
    fm = x[fm_offset:fm_offset + n]
    ffm = x[ffm_offset:ffm_offset + n]

    # Clamp only for reporting; solution itself is unconstrained.
    fm = np.clip(fm, 1.0, None)
    ffm = np.clip(ffm, 50.0, None)
    total = fm + ffm

    out = daily.copy()
    out["fm_lbs_anchor"] = np.round(fm, 2)
    out["ffm_lbs_anchor"] = np.round(ffm, 2)
    out["weight_pred_lbs"] = np.round(total, 2)
    out["fat_pct_anchor"] = np.round(100 * fm / total, 2)
    out_path = ROOT / "analysis" / "AV_anchor_daily_composition.csv"
    out.to_csv(out_path, index=False)

    print("=" * 70)
    print("ANCHOR-ONLY BODY COMPOSITION INTERPOLATION")
    print("=" * 70)
    print(f"\nDays: {n}")
    print(f"Weight observations: {len(weight_idx)}")
    print(f"Composition anchors: {len(comp_idx)}")
    print(f"LSQR iterations: {sol[2]}")
    print(f"LSQR residual norm: {sol[3]:.3f}")

    # Validation on observed weight days.
    weight_sub = out.dropna(subset=["smoothed_weight_lbs"]).copy()
    weight_err = weight_sub["weight_pred_lbs"] - weight_sub["smoothed_weight_lbs"]
    print("\nObserved corrected weight fit:")
    print(f"  RMSE: {np.sqrt(np.mean(weight_err**2)):.3f} lbs")
    print(f"  MAE:  {np.mean(np.abs(weight_err)):.3f} lbs")
    print(f"  Bias: {np.mean(weight_err):+.3f} lbs")

    # Validation on composition anchors.
    anchor_eval = out.merge(
        comp[["date", "fat_mass_lbs", "lean_mass_lbs", "weight_lbs"]],
        on="date",
        how="inner",
    )
    fm_err = anchor_eval["fm_lbs_anchor"] - anchor_eval["fat_mass_lbs"]
    ffm_err = anchor_eval["ffm_lbs_anchor"] - anchor_eval["lean_mass_lbs"]
    wt_err = anchor_eval["weight_pred_lbs"] - anchor_eval["weight_lbs"]
    print("\nComposition anchor fit:")
    print(f"  FM RMSE:   {np.sqrt(np.mean(fm_err**2)):.3f} lbs")
    print(f"  FFM RMSE:  {np.sqrt(np.mean(ffm_err**2)):.3f} lbs")
    print(f"  Wt RMSE:   {np.sqrt(np.mean(wt_err**2)):.3f} lbs")
    print(f"  FM bias:   {np.mean(fm_err):+.3f} lbs")
    print(f"  FFM bias:  {np.mean(ffm_err):+.3f} lbs")

    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
