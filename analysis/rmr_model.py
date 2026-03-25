"""Composition-aware RMR model.

1. Fits individual-specific RMR coefficients from 21 calorimetry measurements
   against body composition (FM/FFM).
2. Interpolates FM/FFM to every day using the Forbes curve for fat/lean
   partitioning during weight change.
3. Computes expected RMR at every day for use in energy balance analysis.

Inputs:
    composition/composition.csv — 49 body composition measurements
    RMR/rmr.csv — 21 calorimetry measurements
    analysis/daily_weight.csv — complete daily weight series

Outputs:
    analysis/daily_composition.csv — daily FM, FFM, fat%, expected RMR
    Prints fitted coefficients and validation report
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

LBS_TO_KG = 0.453592
KG_TO_LBS = 2.20462
FORBES_C_KG = 10.4  # Forbes curve constant for males


def load_data():
    comp = pd.read_csv(ROOT / "composition" / "composition.csv", parse_dates=["date"])
    rmr = pd.read_csv(ROOT / "RMR" / "rmr.csv", parse_dates=["date"])
    daily = pd.read_csv(ROOT / "analysis" / "daily_weight.csv", parse_dates=["date"])
    return comp, rmr, daily


def interpolate_composition(comp, daily):
    """Interpolate FM/FFM to every day using Forbes curve between measurements.

    Between each pair of composition measurements, daily weight changes are
    partitioned into FM/FFM using the Forbes equation:
        dFFM/dBW = C² / (C + FM_kg)²
    where C = 10.4 kg for males.

    This means at high FM, nearly all weight change is fat.
    At low FM, a larger fraction is lean mass.

    At each composition measurement, we snap to the measured values.
    """
    comp = comp.sort_values("date").reset_index(drop=True)

    # Get daily weight series (smoothed = underlying fat mass)
    dates = daily["date"].values
    weights = daily["smoothed_weight_lbs"].values

    n = len(daily)
    fm = np.full(n, np.nan)
    ffm = np.full(n, np.nan)

    # Map composition dates to daily indices
    date_to_idx = {d: i for i, d in enumerate(dates)}
    comp_anchors = []
    for _, row in comp.iterrows():
        idx = date_to_idx.get(row["date"])
        if idx is not None:
            comp_anchors.append({
                "idx": idx,
                "date": row["date"],
                "fm_lbs": row["fat_mass_lbs"],
                "ffm_lbs": row["lean_mass_lbs"],
            })
        else:
            # Find nearest date
            target = np.datetime64(row["date"])
            diffs = np.abs(dates - target)
            nearest = np.argmin(diffs)
            gap = int(diffs[nearest] / np.timedelta64(1, "D"))
            if gap <= 7:
                comp_anchors.append({
                    "idx": nearest,
                    "date": dates[nearest],
                    "fm_lbs": row["fat_mass_lbs"],
                    "ffm_lbs": row["lean_mass_lbs"],
                })

    comp_anchors.sort(key=lambda x: x["idx"])

    # Snap at each anchor
    for a in comp_anchors:
        fm[a["idx"]] = a["fm_lbs"]
        ffm[a["idx"]] = a["ffm_lbs"]

    # Forward-simulate between anchors using Forbes curve
    for i in range(len(comp_anchors) - 1):
        start = comp_anchors[i]
        end = comp_anchors[i + 1]
        s_idx, e_idx = start["idx"], end["idx"]

        # Simulate forward from start anchor
        current_fm = start["fm_lbs"]
        current_ffm = start["ffm_lbs"]

        for j in range(s_idx + 1, e_idx):
            if np.isnan(weights[j]) or np.isnan(weights[j - 1]):
                fm[j] = current_fm
                ffm[j] = current_ffm
                continue

            dw = weights[j] - weights[j - 1]  # daily weight change in lbs
            dw_kg = dw * LBS_TO_KG
            fm_kg = current_fm * LBS_TO_KG

            # Forbes partition: fraction of weight change that is FFM
            ffm_fraction = FORBES_C_KG ** 2 / (FORBES_C_KG + fm_kg) ** 2
            d_ffm = dw * ffm_fraction
            d_fm = dw - d_ffm

            current_fm += d_fm
            current_ffm += d_ffm
            fm[j] = current_fm
            ffm[j] = current_ffm

        # Snap at end anchor
        fm[e_idx] = end["fm_lbs"]
        ffm[e_idx] = end["ffm_lbs"]

    # Extrapolate before first anchor and after last anchor
    first_a = comp_anchors[0]
    last_a = comp_anchors[-1]

    # Before first anchor: simulate backward
    current_fm = first_a["fm_lbs"]
    current_ffm = first_a["ffm_lbs"]
    for j in range(first_a["idx"] - 1, -1, -1):
        if np.isnan(weights[j]) or np.isnan(weights[j + 1]):
            fm[j] = current_fm
            ffm[j] = current_ffm
            continue
        dw = weights[j] - weights[j + 1]  # going backward
        fm_kg = current_fm * LBS_TO_KG
        ffm_fraction = FORBES_C_KG ** 2 / (FORBES_C_KG + fm_kg) ** 2
        d_ffm = dw * ffm_fraction
        d_fm = dw - d_ffm
        current_fm += d_fm
        current_ffm += d_ffm
        fm[j] = current_fm
        ffm[j] = current_ffm

    # After last anchor: simulate forward
    current_fm = last_a["fm_lbs"]
    current_ffm = last_a["ffm_lbs"]
    for j in range(last_a["idx"] + 1, n):
        if np.isnan(weights[j]) or np.isnan(weights[j - 1]):
            fm[j] = current_fm
            ffm[j] = current_ffm
            continue
        dw = weights[j] - weights[j - 1]
        fm_kg = current_fm * LBS_TO_KG
        ffm_fraction = FORBES_C_KG ** 2 / (FORBES_C_KG + fm_kg) ** 2
        d_ffm = dw * ffm_fraction
        d_fm = dw - d_ffm
        current_fm += d_fm
        current_ffm += d_ffm
        fm[j] = current_fm
        ffm[j] = current_ffm

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

    # Least squares fit
    coeffs, residuals, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    a, b, c = coeffs

    print(f"\n=== Fitted RMR Model ===")
    print(f"  RMR = {a:.1f} * FFM_kg + {b:.1f} * FM_kg + {c:.0f}")
    print(f"  (Cunningham: 22.0 * FFM_kg + 0.0 * FM_kg + 500)")
    print(f"  (Hall:       22.0 * FFM_kg + 3.2 * FM_kg + 0)")
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

    path = ROOT / "analysis" / "daily_composition.csv"
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

    print("\nInterpolating FM/FFM to every day (Forbes curve)...")
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
