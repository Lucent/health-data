"""Composition-aware RMR model.

1. Fits individual-specific RMR coefficients from 21 calorimetry measurements
   against body composition (FM/FFM).
2. Interpolates FM/FFM to every day from the composition anchor points.
3. Computes expected RMR at every day for use in energy balance analysis.

Inputs:
    composition/composition.csv — body composition measurements
    RMR/rmr.csv — 21 calorimetry measurements
    analysis/daily_weight.csv — complete daily weight series

Outputs:
    analysis/daily_composition.csv — daily FM, FFM, fat%, expected RMR
    Prints fitted coefficients and validation report
"""

import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

LBS_TO_KG = 0.453592
KG_TO_LBS = 2.20462
FORBES_C_KG = 10.4  # Forbes curve constant for males
RMR_PRIOR = np.array([22.0, 3.2, 500.0])
RMR_RIDGE = 500.0


def load_data():
    comp = pd.read_csv(ROOT / "composition" / "composition.csv", parse_dates=["date"])
    rmr = pd.read_csv(ROOT / "RMR" / "rmr.csv", parse_dates=["date"])
    daily = pd.read_csv(ROOT / "analysis" / "P3_daily_weight.csv", parse_dates=["date"])
    return comp, rmr, daily


def interpolate_composition(comp, daily):
    """Interpolate FM/FFM to every day by linear interpolation between measurements.

    Simple and honest: linearly interpolate fat mass and lean mass between
    the composition anchor points. No Forbes curve, no dependency on
    the daily weight series (which has its own interpolation issues during gaps).

    Before first anchor and after last anchor: hold constant (flat extrapolation).
    """
    comp = comp.sort_values("date").reset_index(drop=True)
    dates = daily["date"].values
    n = len(daily)

    # Build anchor arrays
    date_to_idx = {d: i for i, d in enumerate(dates)}
    anchor_idx = []
    anchor_fm = []
    anchor_ffm = []

    for _, row in comp.iterrows():
        idx = date_to_idx.get(row["date"])
        if idx is None:
            target = np.datetime64(row["date"])
            diffs = np.abs(dates - target)
            nearest = np.argmin(diffs)
            gap = int(diffs[nearest] / np.timedelta64(1, "D"))
            if gap <= 7:
                idx = nearest
        if idx is not None:
            anchor_idx.append(idx)
            anchor_fm.append(row["fat_mass_lbs"])
            anchor_ffm.append(row["lean_mass_lbs"])

    anchor_idx = np.array(anchor_idx)
    anchor_fm = np.array(anchor_fm)
    anchor_ffm = np.array(anchor_ffm)

    # Linear interpolation
    all_idx = np.arange(n)
    fm = np.interp(all_idx, anchor_idx, anchor_fm)
    ffm = np.interp(all_idx, anchor_idx, anchor_ffm)

    return fm, ffm


def fit_rmr_coefficients(comp, rmr, daily, fm_arr, ffm_arr):
    """Fit RMR = a * FFM_kg + b * FM_kg + c using least squares.

    For each RMR measurement, find the nearest composition-interpolated
    FM/FFM values.
    """
    dates = daily["date"].values
    date_to_idx = {d: i for i, d in enumerate(dates)}

    X = []
    y = []
    labels = []
    for _, row in rmr.iterrows():
        idx = date_to_idx.get(row["date"])
        if idx is None:
            target = np.datetime64(row["date"])
            diffs = np.abs(dates - target)
            idx = np.argmin(diffs)
            gap = int(diffs[idx] / np.timedelta64(1, "D"))
            if gap > 14:
                continue

        if np.isnan(fm_arr[idx]) or np.isnan(ffm_arr[idx]):
            continue

        fm_kg = fm_arr[idx] * LBS_TO_KG
        ffm_kg = ffm_arr[idx] * LBS_TO_KG
        X.append([ffm_kg, fm_kg, 1.0])
        y.append(row["rmr_kcal"])
        labels.append((row["date"].strftime("%Y-%m-%d"), row["rmr_kcal"],
                        fm_arr[idx], ffm_arr[idx], row["device"]))

    X = np.array(X)
    y = np.array(y)

    # Ordinary least squares is unstable here because FM and FFM are highly
    # collinear across only 21 measurements. Use a physiologic ridge prior so
    # the fitted coefficients stay plausible and the Kalman prior is not driven
    # by sign-flipped slopes.
    ols_coeffs, *_ = np.linalg.lstsq(X, y, rcond=None)
    ridge_penalty = RMR_RIDGE * np.eye(3)
    coeffs = np.linalg.solve(X.T @ X + ridge_penalty, X.T @ y + ridge_penalty @ RMR_PRIOR)
    a, b, c = coeffs

    print(f"\n=== Fitted RMR Model ===")
    print(f"  RMR = {a:.1f} * FFM_kg + {b:.1f} * FM_kg + {c:.0f}")
    print(f"  (Cunningham: 22.0 * FFM_kg + 0.0 * FM_kg + 500)")
    print(f"  (Hall:       22.0 * FFM_kg + 3.2 * FM_kg + 0)")
    print(f"  Ridge prior: {RMR_PRIOR[0]:.1f} * FFM_kg + {RMR_PRIOR[1]:.1f} * FM_kg + {RMR_PRIOR[2]:.0f}")
    print(f"  OLS (unstable): {ols_coeffs[0]:.1f} * FFM_kg + {ols_coeffs[1]:.1f} * FM_kg + {ols_coeffs[2]:.0f}")
    print(f"  Fitted on {len(y)} measurements")

    # Per-measurement validation
    print(f"\n{'Date':>12} {'Meas':>5} {'Pred':>5} {'Err':>5} {'FM':>5} {'FFM':>5} {'Device':>15}")
    errors = []
    for i, (date, meas, fm_v, ffm_v, device) in enumerate(labels):
        pred = X[i] @ coeffs
        err = pred - meas
        errors.append(err)
        print(f"{date:>12} {meas:5.0f} {pred:5.0f} {err:+5.0f} "
              f"{fm_v:5.1f} {ffm_v:5.1f} {device:>15}")

    errors = np.array(errors)
    print(f"\n  RMSE: {np.sqrt(np.mean(errors**2)):.0f} kcal/day")
    print(f"  MAE:  {np.mean(np.abs(errors)):.0f} kcal/day")
    print(f"  Bias: {np.mean(errors):+.0f} kcal/day")

    return a, b, c


def compute_daily_rmr(fm_arr, ffm_arr, a, b, c):
    """Compute expected RMR at every day using fitted coefficients."""
    fm_kg = fm_arr * LBS_TO_KG
    ffm_kg = ffm_arr * LBS_TO_KG
    return a * ffm_kg + b * fm_kg + c


def save_output(daily, fm_arr, ffm_arr, rmr_arr):
    """Save daily composition CSV."""
    out = daily[["date"]].copy()
    out["weight_lbs"] = daily["smoothed_weight_lbs"]
    out["fm_lbs"] = np.round(fm_arr, 1)
    out["ffm_lbs"] = np.round(ffm_arr, 1)
    out["fat_pct"] = np.round(fm_arr / (fm_arr + ffm_arr) * 100, 1)
    out["expected_rmr"] = np.round(rmr_arr, 0)
    out["date"] = out["date"].dt.strftime("%Y-%m-%d")

    path = ROOT / "analysis" / "P2_daily_composition.csv"
    out.to_csv(path, index=False)
    print(f"\nWrote {len(out)} rows to {path}")


def validate_against_formulas(daily, fm_arr, ffm_arr, a, b, c):
    """Compare fitted model vs MSJ vs Cunningham across weight range."""
    height_cm = 5 * 30.48 + 11.75 * 2.54

    print(f"\n=== Model comparison by weight bin ===")
    print(f"{'Wt bin':>7} {'Fitted':>7} {'Cunn':>7} {'MSJ':>7} {'FM':>6} {'FFM':>6} {'Fat%':>5}")

    valid = ~np.isnan(fm_arr) & ~np.isnan(ffm_arr)
    weights = daily["smoothed_weight_lbs"].values
    dates = daily["date"].values

    for wt_bin in [170, 180, 190, 200, 210, 220, 230]:
        mask = valid & (weights >= wt_bin - 5) & (weights < wt_bin + 5)
        if mask.sum() < 10:
            continue
        fm_med = np.median(fm_arr[mask])
        ffm_med = np.median(ffm_arr[mask])
        fm_kg = fm_med * LBS_TO_KG
        ffm_kg = ffm_med * LBS_TO_KG
        fitted = a * ffm_kg + b * fm_kg + c
        cunningham = 500 + 22 * ffm_kg
        age_med = np.median([(d - np.datetime64("1982-10-15")) / np.timedelta64(1, "D") / 365.25
                             for d in dates[mask]])
        msj = 10 * (wt_bin * LBS_TO_KG) + 6.25 * height_cm - 5 * age_med + 5
        fat_pct = fm_med / (fm_med + ffm_med) * 100
        print(f"{wt_bin:>7} {fitted:7.0f} {cunningham:7.0f} {msj:7.0f} "
              f"{fm_med:6.1f} {ffm_med:6.1f} {fat_pct:5.1f}")


def main():
    print("Loading data...")
    comp, rmr, daily = load_data()
    print(f"  Composition: {len(comp)} measurements")
    print(f"  RMR: {len(rmr)} measurements")
    print(f"  Daily weight: {len(daily)} days")

    print("\nInterpolating FM/FFM to every day...")
    fm_arr, ffm_arr = interpolate_composition(comp, daily)
    filled = (~np.isnan(fm_arr)).sum()
    print(f"  Filled {filled}/{len(daily)} days ({filled/len(daily)*100:.1f}%)")

    # Validate interpolation at composition anchors
    print(f"\n=== Composition interpolation validation ===")
    dates = daily["date"].values
    max_err = 0
    for _, row in comp.iterrows():
        target = np.datetime64(row["date"])
        diffs = np.abs(dates - target)
        idx = np.argmin(diffs)
        gap = int(diffs[idx] / np.timedelta64(1, "D"))
        if gap > 7 or np.isnan(fm_arr[idx]):
            continue
        fm_err = abs(fm_arr[idx] - row["fat_mass_lbs"])
        ffm_err = abs(ffm_arr[idx] - row["lean_mass_lbs"])
        max_err = max(max_err, fm_err, ffm_err)
    print(f"  Max error at composition anchors: {max_err:.2f} lbs")

    print("\nFitting RMR coefficients...")
    a, b, c = fit_rmr_coefficients(comp, rmr, daily, fm_arr, ffm_arr)

    rmr_arr = compute_daily_rmr(fm_arr, ffm_arr, a, b, c)
    validate_against_formulas(daily, fm_arr, ffm_arr, a, b, c)
    save_output(daily, fm_arr, ffm_arr, rmr_arr)


if __name__ == "__main__":
    main()
