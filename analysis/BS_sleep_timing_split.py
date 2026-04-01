#!/usr/bin/env python3
"""BS. Early vs late bedtime split and cutoff sweep.

Question:
  If sleep timing is split into earlier vs later nights, does "going to bed
  early" correlate with any metabolic or behavioral variable, and where is the
  most informative cutoff?

Method:
  - Use Samsung Health sleep records, one primary sleep per wake date.
  - Convert bedtime to a monotonic night clock:
      00:30 -> 24.5, 03:00 -> 27.0, 23:30 -> 23.5
  - Primary analysis: median split (half early, half late as closely as
    possible).
  - Secondary analysis: sweep candidate cutoffs across the middle 60% of the
    bedtime distribution.
  - Report both raw point-biserial correlation and partial correlation after
    controlling for sleep duration, fat mass, weekend, and year. Sleep duration
    itself is reported without controlling for sleep duration.
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT_CSV = ROOT / "analysis" / "BS_sleep_timing_cutoff_sweep.csv"

TARGETS = [
    ("calories", "Same-day calories"),
    ("surplus", "Same-day surplus"),
    ("steps", "Same-day steps"),
    ("tdee", "Same-day TDEE"),
    ("cal_next", "Next-day calories"),
    ("steps_next", "Next-day steps"),
    ("weight_change", "Weight change"),
    ("sleep_hours", "Sleep duration"),
]

PROXY_TARGETS = [
    ("tdee_resid", "TDEE residual"),
    ("tdee_ratio", "TDEE/RMR ratio"),
]

BED_PROXY_WINDOWS = [
    ("bed_hour", "Same-day bedtime"),
    ("bed_3d", "3d bedtime"),
    ("bed_7d", "7d bedtime"),
    ("bed_14d", "14d bedtime"),
    ("bed_30d", "30d bedtime"),
]


def bedtime_to_night_clock(series: pd.Series) -> pd.Series:
    hour = pd.to_timedelta(series + ":00").dt.total_seconds() / 3600
    return pd.Series(np.where(hour < 12, hour + 24, hour), index=series.index)


def fmt_hour(hour: float) -> str:
    h = int(np.floor(hour)) % 24
    m = int(round((hour - np.floor(hour)) * 60))
    if m == 60:
        h = (h + 1) % 24
        m = 0
    return f"{h:02d}:{m:02d}"


def partial_corr(y: np.ndarray, x: np.ndarray, controls: np.ndarray) -> tuple[float, int]:
    mask = np.isfinite(y) & np.isfinite(x)
    if controls.size:
        for col in controls.T:
            mask &= np.isfinite(col)
    y = y[mask]
    x = x[mask]
    c = controls[mask] if controls.size else np.empty((len(y), 0))
    if len(y) < 50:
        return np.nan, len(y)

    xmat = np.column_stack([np.ones(len(y)), c])
    beta_y = np.linalg.lstsq(xmat, y, rcond=None)[0]
    beta_x = np.linalg.lstsq(xmat, x, rcond=None)[0]
    res_y = y - xmat @ beta_y
    res_x = x - xmat @ beta_x
    return np.corrcoef(res_x, res_y)[0, 1], len(y)


def point_biserial(x: np.ndarray, y: np.ndarray) -> tuple[float, int]:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 50:
        return np.nan, int(mask.sum())
    return np.corrcoef(x[mask], y[mask])[0, 1], int(mask.sum())


def build_daily() -> pd.DataFrame:
    sleep = pd.read_csv(ROOT / "steps-sleep" / "sleep.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    steps = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    weight = pd.read_csv(ROOT / "weight" / "weight.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    temp = pd.read_csv(ROOT / "analysis" / "AS_temperature_baseline.csv", parse_dates=["date"])

    sl = sleep[["date", "sleep_start", "sleep_end", "sleep_hours"]].drop_duplicates("date", keep="last").copy()
    sl["bed_hour"] = bedtime_to_night_clock(sl["sleep_start"])
    sl["wake_hour"] = pd.to_timedelta(sl["sleep_end"] + ":00").dt.total_seconds() / 3600
    sl = sl.sort_values("date").reset_index(drop=True)
    for win in [1, 3, 7, 14, 30]:
        sl[f"bed_{win}d"] = sl["bed_hour"].rolling(win, min_periods=max(1, win // 2)).mean()

    df = sl.merge(intake[["date", "calories", "protein_g", "carbs_g", "fat_g"]], on="date", how="left")
    df = df.merge(steps[["date", "steps"]], on="date", how="left")
    df = df.merge(kalman[["date", "tdee", "fat_mass_lbs"]], on="date", how="left")
    df = df.merge(weight[["date", "weight_lbs"]].drop_duplicates("date", keep="first"), on="date", how="left")
    df = df.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    df = df.merge(temp[["date", "baseline_mean", "n_readings"]], on="date", how="left")
    df = df.sort_values("date").reset_index(drop=True)

    df["surplus"] = df["calories"] - df["tdee"]
    df["tdee_resid"] = df["tdee"] - df["expected_rmr"]
    df["tdee_ratio"] = df["tdee"] / df["expected_rmr"]
    df["cal_next"] = df["calories"].shift(-1)
    df["steps_next"] = df["steps"].shift(-1)
    df["weight_change"] = df["weight_lbs"] - df["weight_lbs"].shift(1)
    df["weekend"] = (df["date"].dt.dayofweek >= 5).astype(float)
    df["year"] = df["date"].dt.year.astype(float)
    return df


def control_matrix(df: pd.DataFrame, include_sleep_hours: bool) -> np.ndarray:
    cols = []
    if include_sleep_hours:
        cols.append("sleep_hours")
    cols.extend(["fat_mass_lbs", "weekend", "year"])
    return df[cols].values.astype(float)


def summarize_split(df: pd.DataFrame, cutoff: float, label: str) -> pd.DataFrame:
    out = []
    early = (df["bed_hour"] <= cutoff).astype(float).values

    for target, target_label in TARGETS:
        y = df[target].values.astype(float)
        raw_r, raw_n = point_biserial(early, y)
        valid = np.isfinite(early) & np.isfinite(y)
        if valid.sum() < 50:
            continue
        early_vals = y[(early == 1) & np.isfinite(y)]
        late_vals = y[(early == 0) & np.isfinite(y)]
        controls = control_matrix(df, include_sleep_hours=(target != "sleep_hours"))
        partial_r, partial_n = partial_corr(y, early, controls)
        out.append({
            "analysis": label,
            "cutoff_hour": cutoff,
            "cutoff_hhmm": fmt_hour(cutoff),
            "target": target,
            "target_label": target_label,
            "early_n": int(len(early_vals)),
            "late_n": int(len(late_vals)),
            "early_mean": float(np.nanmean(early_vals)),
            "late_mean": float(np.nanmean(late_vals)),
            "diff_early_minus_late": float(np.nanmean(early_vals) - np.nanmean(late_vals)),
            "raw_r": raw_r,
            "raw_n": raw_n,
            "partial_r": partial_r,
            "partial_n": partial_n,
        })

    return pd.DataFrame(out)


def summarize_debt_and_interaction(df: pd.DataFrame) -> None:
    sleep_ref = float(df["sleep_hours"].median())
    for win in [3, 7, 14, 30]:
        avg = df["sleep_hours"].rolling(win, min_periods=max(1, win // 2)).mean()
        df[f"sleep_{win}d"] = avg
        df[f"sleep_debt_{win}d"] = sleep_ref - avg

    df["late_bed"] = (df["bed_hour"] > df["bed_hour"].median()).astype(float)
    df["short_sleep"] = (df["sleep_hours"] < sleep_ref).astype(float)
    df["late_short"] = df["late_bed"] * df["short_sleep"]

    print("\nSleep debt correlations:")
    print(f"{'Metric':>18} {'Target':>16} {'raw_r':>8} {'partial_r':>10} {'n':>6}")
    debt_targets = [
        ("calories", "Same-day calories"),
        ("surplus", "Same-day surplus"),
        ("steps", "Same-day steps"),
        ("tdee_resid", "TDEE residual"),
        ("tdee_ratio", "TDEE/RMR ratio"),
        ("cal_next", "Next-day calories"),
        ("steps_next", "Next-day steps"),
    ]
    for metric in ["sleep_debt_3d", "sleep_debt_7d", "sleep_debt_14d", "sleep_debt_30d"]:
        for target, target_label in debt_targets:
            raw_r, raw_n = point_biserial(df[metric].values.astype(float), df[target].values.astype(float))
            controls = df[["fat_mass_lbs", "year", "weekend", "bed_hour"]].values.astype(float)
            partial_r, partial_n = partial_corr(
                df[target].values.astype(float),
                df[metric].values.astype(float),
                controls,
            )
            if partial_n < 50:
                continue
            print(f"{metric:>18} {target_label:>16} {raw_r:+8.3f} {partial_r:+10.3f} {partial_n:6d}")

    print("\nLate-bed x short-sleep interaction:")
    print(f"Reference median sleep = {sleep_ref:.2f}h; late bed = after {fmt_hour(float(df['bed_hour'].median()))}")
    print(f"{'Target':>18} {'late_short diff':>16} {'raw_r':>8} {'partial_r':>10} {'n':>6}")
    for target, target_label in debt_targets:
        y = df[target].values.astype(float)
        x = df["late_short"].values.astype(float)
        valid = np.isfinite(y) & np.isfinite(x)
        if valid.sum() < 100:
            continue
        group = y[(x == 1) & np.isfinite(y)]
        rest = y[(x == 0) & np.isfinite(y)]
        diff = float(np.nanmean(group) - np.nanmean(rest)) if len(group) and len(rest) else np.nan
        raw_r, _ = point_biserial(x, y)
        controls = df[["fat_mass_lbs", "year", "weekend", "bed_hour", "sleep_hours"]].values.astype(float)
        partial_r, partial_n = partial_corr(y, x, controls)
        if partial_n < 50:
            continue
        print(f"{target_label:>18} {diff:+16.1f} {raw_r:+8.3f} {partial_r:+10.3f} {partial_n:6d}")


def main() -> None:
    df = build_daily()
    median_cutoff = float(df["bed_hour"].median())
    split = summarize_split(df, median_cutoff, "median_split")

    candidates = sorted(df["bed_hour"].dropna().quantile(np.linspace(0.2, 0.8, 25)).unique())
    sweep_rows = []
    for cutoff in candidates:
        early = (df["bed_hour"] <= cutoff).astype(float).values
        for target, target_label in TARGETS:
            y = df[target].values.astype(float)
            valid = np.isfinite(early) & np.isfinite(y)
            if valid.sum() < 100:
                continue
            early_n = int((early[valid] == 1).sum())
            late_n = int((early[valid] == 0).sum())
            if min(early_n, late_n) < 200:
                continue
            raw_r, raw_n = point_biserial(early, y)
            controls = control_matrix(df, include_sleep_hours=(target != "sleep_hours"))
            partial_r, partial_n = partial_corr(y, early, controls)
            sweep_rows.append({
                "analysis": "cutoff_sweep",
                "cutoff_hour": cutoff,
                "cutoff_hhmm": fmt_hour(cutoff),
                "target": target,
                "target_label": target_label,
                "early_n": early_n,
                "late_n": late_n,
                "early_mean": float(np.nanmean(y[(early == 1) & np.isfinite(y)])),
                "late_mean": float(np.nanmean(y[(early == 0) & np.isfinite(y)])),
                "diff_early_minus_late": float(np.nanmean(y[(early == 1) & np.isfinite(y)]) - np.nanmean(y[(early == 0) & np.isfinite(y)])),
                "raw_r": raw_r,
                "raw_n": raw_n,
                "partial_r": partial_r,
                "partial_n": partial_n,
            })

    sweep = pd.DataFrame(sweep_rows)
    out = pd.concat([split, sweep], ignore_index=True)
    out.to_csv(OUT_CSV, index=False)

    proxy_rows = []
    for proxy_col, proxy_label in [("baseline_mean", "Baseline temp")] + BED_PROXY_WINDOWS:
        for target, target_label in PROXY_TARGETS:
            raw_r, raw_n = point_biserial(df[proxy_col].values.astype(float), df[target].values.astype(float))
            if "temp" in proxy_col:
                controls = df[["fat_mass_lbs", "year", "weekend", "n_readings"]].values.astype(float)
            else:
                controls = df[["fat_mass_lbs", "year", "weekend", "sleep_hours"]].values.astype(float)
            partial_r, partial_n = partial_corr(
                df[target].values.astype(float),
                df[proxy_col].values.astype(float),
                controls,
            )
            proxy_rows.append({
                "proxy": proxy_label,
                "target": target_label,
                "raw_r": raw_r,
                "partial_r": partial_r,
                "n": partial_n if partial_n else raw_n,
            })
    proxy_df = pd.DataFrame(proxy_rows)

    print("=" * 84)
    print("SLEEP TIMING: EARLY VS LATE BEDTIME")
    print("=" * 84)
    print(f"Sleep rows: {len(df)}")
    print(f"Bedtime median split cutoff: {fmt_hour(median_cutoff)} ({median_cutoff:.2f} on night clock)")
    print(f"Bedtime quartiles: Q1={fmt_hour(df['bed_hour'].quantile(0.25))}, "
          f"median={fmt_hour(median_cutoff)}, Q3={fmt_hour(df['bed_hour'].quantile(0.75))}")

    split_show = split.sort_values("target")[[
        "target_label", "early_n", "late_n", "diff_early_minus_late", "raw_r", "partial_r"
    ]]
    print("\nMedian split results:")
    print(split_show.to_string(index=False, formatters={
        "diff_early_minus_late": lambda x: f"{x:+.1f}",
        "raw_r": lambda x: f"{x:+.3f}",
        "partial_r": lambda x: f"{x:+.3f}",
    }))

    print("\nBest cutoff by |partial_r| for each target:")
    for target, target_label in TARGETS:
        sub = sweep[sweep["target"] == target].copy()
        if len(sub) == 0:
            continue
        best = sub.iloc[sub["partial_r"].abs().argmax()]
        print(
            f"  {target_label:>18}: cutoff {best['cutoff_hhmm']}  "
            f"partial_r={best['partial_r']:+.3f}  "
            f"diff={best['diff_early_minus_late']:+.1f}  "
            f"early_n={int(best['early_n'])} late_n={int(best['late_n'])}"
        )

    bed = df["bed_hour"].values.astype(float)
    print("\nContinuous bedtime correlations (later bedtime = larger number):")
    for target, target_label in TARGETS:
        y = df[target].values.astype(float)
        raw_r, raw_n = point_biserial(bed, y)
        controls = control_matrix(df, include_sleep_hours=(target != "sleep_hours"))
        partial_r, partial_n = partial_corr(y, bed, controls)
        if raw_n < 50:
            continue
        print(f"  {target_label:>18}: raw_r={raw_r:+.3f}  partial_r={partial_r:+.3f}  n={partial_n}")

    print("\nMetabolic proxy correlations:")
    print(f"{'Proxy':>18} {'Target':>16} {'raw_r':>8} {'partial_r':>10} {'n':>6}")
    for _, row in proxy_df.iterrows():
        print(f"{row['proxy']:>18} {row['target']:>16} {row['raw_r']:+8.3f} {row['partial_r']:+10.3f} {int(row['n']):6d}")

    summarize_debt_and_interaction(df)

    print(f"\nArtifact: {OUT_CSV}")


if __name__ == "__main__":
    main()
