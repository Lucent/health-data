#!/usr/bin/env python3
"""AX. Anchor-driven daily FM/FFM with weak intake plausibility prior.

Primary constraints:
  - corrected scale weight observations
  - direct FM / FFM composition anchors

Secondary weak priors:
  - FM/FFM smoothness
  - latent TDEE smoothness
  - weak energy balance consistency: intake - TDEE ~= 3500 * dFM
  - broad TDEE prior from lean mass (Cunningham-style)

Intake is used only as a soft plausibility term, not as the dominant state
update. This is a compromise between the fully non-intake AV interpolator and
the stronger intake-driven Kalman pipeline.
"""

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import lsqr

ROOT = Path(__file__).resolve().parent.parent
LBS_TO_KG = 0.453592

# Primary observation weights
W_WEIGHT = 4.0
W_COMP_FM = 18.0
W_COMP_FFM = 18.0

# Shape priors
L2_FM_SECOND = 0.14
L2_FFM_SECOND = 1.00
L2_FM_FIRST = 0.010
L2_FFM_FIRST = 0.060
L2_PRIOR_FM = 0.04
L2_PRIOR_FFM = 0.08

# Weak intake / TDEE priors
W_ENERGY = 0.0015
W_TDEE_PRIOR = 0.0008
L2_TDEE_FIRST = 0.0020
L2_TDEE_SECOND = 0.0100


def add_row_triplets(rows, cols, data, r, cvals):
    for c, v in cvals:
        rows.append(r)
        cols.append(c)
        data.append(v)


def interp_prior(daily_dates, comp_dates, values):
    x = np.array([(d - daily_dates[0]) / np.timedelta64(1, "D") for d in daily_dates], dtype=float)
    xc = np.array([(d - daily_dates[0]) / np.timedelta64(1, "D") for d in comp_dates], dtype=float)
    return np.interp(x, xc, values)


def load_data():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    weight = pd.read_csv(ROOT / "analysis" / "P1_smoothed_weight.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "composition" / "composition.csv", parse_dates=["date"])
    return intake, weight, comp


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--w-energy", type=float, default=W_ENERGY)
    parser.add_argument("--w-tdee-prior", type=float, default=W_TDEE_PRIOR)
    parser.add_argument("--l2-tdee-first", type=float, default=L2_TDEE_FIRST)
    parser.add_argument("--l2-tdee-second", type=float, default=L2_TDEE_SECOND)
    parser.add_argument(
        "--out-file",
        default=str(ROOT / "analysis" / "AX_weakintake_daily_composition.csv"),
        help="Output CSV path",
    )
    args = parser.parse_args()

    intake, weight, comp = load_data()
    daily = pd.DataFrame({"date": pd.date_range(intake["date"].min(), intake["date"].max(), freq="D")})
    daily = daily.merge(intake[["date", "calories"]], on="date", how="left")
    daily["calories"] = daily["calories"].fillna(0)
    daily = daily.merge(weight[["date", "smoothed_weight_lbs"]], on="date", how="left")

    comp = comp[["date", "fat_mass_lbs", "lean_mass_lbs", "weight_lbs"]].dropna(
        subset=["fat_mass_lbs", "lean_mass_lbs"]
    ).copy()
    comp = comp.sort_values("date").drop_duplicates("date", keep="last").reset_index(drop=True)

    n = len(daily)
    fm_offset = 0
    ffm_offset = n
    tdee_offset = 2 * n

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

    # Weight observations: FM + FFM = corrected weight.
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

    # Weak priors to composition interpolation.
    for i in range(n):
        s_fm = np.sqrt(L2_PRIOR_FM)
        s_ffm = np.sqrt(L2_PRIOR_FFM)
        add_row_triplets(rows, cols, data, r, [(fm_offset + i, s_fm)])
        b.append(s_fm * fm_prior[i])
        r += 1
        add_row_triplets(rows, cols, data, r, [(ffm_offset + i, s_ffm)])
        b.append(s_ffm * ffm_prior[i])
        r += 1

    # Smoothness on FM and FFM.
    for i in range(1, n):
        s_fm = np.sqrt(L2_FM_FIRST)
        s_ffm = np.sqrt(L2_FFM_FIRST)
        add_row_triplets(rows, cols, data, r, [(fm_offset + i - 1, -s_fm), (fm_offset + i, s_fm)])
        b.append(0.0)
        r += 1
        add_row_triplets(rows, cols, data, r, [(ffm_offset + i - 1, -s_ffm), (ffm_offset + i, s_ffm)])
        b.append(0.0)
        r += 1

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

    # Weak energy-balance consistency:
    # calories_t - TDEE_t = 3500 * (FM_{t+1} - FM_t)
    for i in range(n - 1):
        s = np.sqrt(args.w_energy)
        add_row_triplets(
            rows, cols, data, r,
            [(fm_offset + i, 3500 * s), (fm_offset + i + 1, -3500 * s), (tdee_offset + i, -s)]
        )
        b.append(-s * daily["calories"].iloc[i])
        r += 1

    # Broad TDEE prior from lean mass only.
    # Cunningham: 500 + 22 * FFM_kg = 500 + 9.98 * FFM_lbs.
    for i in range(n):
        s = np.sqrt(args.w_tdee_prior)
        add_row_triplets(rows, cols, data, r, [(tdee_offset + i, s), (ffm_offset + i, -22 * LBS_TO_KG * s)])
        b.append(500 * s)
        r += 1

    # TDEE smoothness.
    for i in range(1, n):
        s = np.sqrt(args.l2_tdee_first)
        add_row_triplets(rows, cols, data, r, [(tdee_offset + i - 1, -s), (tdee_offset + i, s)])
        b.append(0.0)
        r += 1

    for i in range(1, n - 1):
        s = np.sqrt(args.l2_tdee_second)
        add_row_triplets(
            rows, cols, data, r,
            [(tdee_offset + i - 1, s), (tdee_offset + i, -2 * s), (tdee_offset + i + 1, s)]
        )
        b.append(0.0)
        r += 1

    A = sparse.csr_matrix((data, (rows, cols)), shape=(r, 3 * n))
    b = np.asarray(b)

    sol = lsqr(A, b, atol=1e-8, btol=1e-8, iter_lim=30000)
    x = sol[0]
    fm = np.clip(x[fm_offset:ffm_offset], 1.0, None)
    ffm = np.clip(x[ffm_offset:tdee_offset], 50.0, None)
    tdee = x[tdee_offset:tdee_offset + n]

    out = daily.copy()
    out["fm_lbs_weak"] = np.round(fm, 2)
    out["ffm_lbs_weak"] = np.round(ffm, 2)
    out["weight_pred_lbs"] = np.round(fm + ffm, 2)
    out["fat_pct_weak"] = np.round(100 * fm / (fm + ffm), 2)
    out["tdee_weak"] = np.round(tdee, 1)
    out_path = Path(args.out_file)
    out.to_csv(out_path, index=False)

    print("=" * 70)
    print("WEAK-INTAKE BODY COMPOSITION INTERPOLATION")
    print("=" * 70)
    print(
        f"\nParams: w_energy={args.w_energy}  w_tdee_prior={args.w_tdee_prior}  "
        f"l2_tdee_first={args.l2_tdee_first}  l2_tdee_second={args.l2_tdee_second}"
    )
    print(f"\nDays: {n}")
    print(f"Weight observations: {len(weight_idx)}")
    print(f"Composition anchors: {len(comp_idx)}")
    print(f"LSQR iterations: {sol[2]}")
    print(f"LSQR residual norm: {sol[3]:.3f}")

    weight_sub = out.dropna(subset=["smoothed_weight_lbs"]).copy()
    weight_err = weight_sub["weight_pred_lbs"] - weight_sub["smoothed_weight_lbs"]
    print("\nObserved corrected weight fit:")
    print(f"  RMSE: {np.sqrt(np.mean(weight_err**2)):.3f} lbs")
    print(f"  MAE:  {np.mean(np.abs(weight_err)):.3f} lbs")
    print(f"  Bias: {np.mean(weight_err):+.3f} lbs")

    anchor_eval = out.merge(
        comp[["date", "fat_mass_lbs", "lean_mass_lbs", "weight_lbs"]],
        on="date",
        how="inner",
    )
    fm_err = anchor_eval["fm_lbs_weak"] - anchor_eval["fat_mass_lbs"]
    ffm_err = anchor_eval["ffm_lbs_weak"] - anchor_eval["lean_mass_lbs"]
    wt_err = anchor_eval["weight_pred_lbs"] - anchor_eval["weight_lbs"]
    print("\nComposition anchor fit:")
    print(f"  FM RMSE:   {np.sqrt(np.mean(fm_err**2)):.3f} lbs")
    print(f"  FFM RMSE:  {np.sqrt(np.mean(ffm_err**2)):.3f} lbs")
    print(f"  Wt RMSE:   {np.sqrt(np.mean(wt_err**2)):.3f} lbs")

    implied = out["calories"].values[:-1] - 3500 * np.diff(out["fm_lbs_weak"].values)
    print("\nImplied TDEE from intake and dFM:")
    print(f"  Mean: {np.nanmean(implied):.1f}")
    print(f"  SD:   {np.nanstd(implied):.1f}")
    print(f"  1st-99th pct: {np.nanpercentile(implied, 1):.0f} to {np.nanpercentile(implied, 99):.0f}")

    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
