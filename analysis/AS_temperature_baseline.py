#!/usr/bin/env python3
"""AS. Reconstruct daily baseline temperature from individual readings.

AR showed the circadian cycle has 0.14°F amplitude (wake-anchored, 2-harmonic).
A reading of 98.0°F at 10pm (near peak) implies a lower baseline than 98.0°F
at 10am (near trough). This script:

1. Fits the circadian curve: temp(h) = baseline + f(hours_since_wake)
2. For each reading, computes: baseline = reading - f(h)
3. Validates: does the baseline have less variance than raw daily means?
4. Re-tests all metabolic correlates with the baseline estimate.
5. Outputs a plot of the circadian curve shape.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def parse_utc_offset_hours(offset_str):
    if pd.isna(offset_str):
        return np.nan
    s = str(offset_str).replace("UTC", "")
    sign = -1 if s[0] == "-" else 1
    s = s.lstrip("+-")
    return sign * (int(s[:2]) + int(s[2:]) / 60)


def home_offset(ts):
    m = ts.month
    return -4.0 if 4 <= m <= 10 else -5.0


def asymmetric_ema(fm_vals, hl_up, hl_down):
    alpha_up = 1 - np.exp(-np.log(2) / hl_up)
    alpha_down = 1 - np.exp(-np.log(2) / hl_down)
    sp = np.empty(len(fm_vals))
    sp[0] = fm_vals[0]
    for i in range(1, len(fm_vals)):
        prev, cur = sp[i - 1], fm_vals[i]
        if np.isnan(prev): sp[i] = cur
        elif np.isnan(cur): sp[i] = prev
        elif cur > prev: sp[i] = prev + alpha_up * (cur - prev)
        else: sp[i] = prev + alpha_down * (cur - prev)
    return sp


def main():
    np.random.seed(42)

    # ── Load temperature (with pre-computed local_hour from merge.py) ──
    temp = pd.read_csv(ROOT / "temperature" / "temperature.csv", parse_dates=["date"])
    temp = temp[(temp["temp_f"] > 95) & (temp["temp_f"] < 101)].copy()
    temp.rename(columns={"date": "timestamp"}, inplace=True)
    # Day boundary at 5am: readings before 5am belong to the previous day's waking period.
    # Matches sleep.csv convention (date = wake-up date).
    temp["date"] = (temp["timestamp"] - pd.Timedelta(hours=5)).dt.floor("D")

    # Load sleep for wake/sleep times (local times, already timezone-correct)
    sleep = pd.read_csv(ROOT / "steps-sleep" / "sleep.csv", parse_dates=["date"])
    sleep["wake_hour"] = pd.to_timedelta(sleep["sleep_end"] + ":00").dt.total_seconds() / 3600
    sleep["sleep_start_hour"] = pd.to_timedelta(sleep["sleep_start"] + ":00").dt.total_seconds() / 3600

    temp = temp.merge(sleep[["date", "wake_hour", "sleep_start_hour"]],
                      on="date", how="left")

    # Load sunrise times (local, from sunlight.csv)
    sunlight = pd.read_csv(ROOT / "steps-sleep" / "sunlight.csv", parse_dates=["date"])
    sunlight["sunrise_hour"] = pd.to_timedelta(sunlight["sunrise"] + ":00").dt.total_seconds() / 3600
    temp = temp.merge(sunlight[["date", "sunrise_hour"]], on="date", how="left")

    # Compute hours since wake with wrap-around for post-midnight readings.
    # Day boundary is 5am: a 2am reading assigned to date D woke at wake_hour on date D.
    # If local_hour < wake_hour, the reading is post-midnight, so:
    #   hours_since_wake = (24 - wake_hour) + local_hour
    def hours_since(local_h, ref_h):
        """Hours elapsed from ref_h to local_h, wrapping past midnight."""
        diff = local_h - ref_h
        return np.where(diff < -5, diff + 24, diff)  # wrap if > 5h negative

    temp["hours_since_wake"] = hours_since(temp["local_hour"].values,
                                           temp["wake_hour"].values)
    temp["hours_since_sunrise"] = hours_since(temp["local_hour"].values,
                                              temp["sunrise_hour"].values)

    # Surface readings during claimed sleep for timezone sanity check (don't discard)
    has_wake = temp["wake_hour"].notna()
    during_sleep = has_wake & ((temp["hours_since_wake"] < 0) | (temp["hours_since_wake"] > 20))

    n_sleep = during_sleep.sum()
    print(f"  Readings during claimed sleep (timezone check): {n_sleep}")
    if n_sleep > 0:
        sleep_readings = temp[during_sleep][["timestamp", "temp_f", "local_hour", "wake_hour",
                                              "hours_since_wake"]].copy()
        print(f"  {'Timestamp':>22} {'Temp':>5} {'Local hr':>9} {'Wake hr':>8} {'H since wake':>13}")
        for _, row in sleep_readings.head(10).iterrows():
            print(f"  {str(row['timestamp']):>22} {row['temp_f']:>5.1f} {row['local_hour']:>9.1f} "
                  f"{row['wake_hour']:>8.1f} {row['hours_since_wake']:>+13.1f}")

    # Keep ALL readings for curve fitting — don't filter
    # Readings with wake time use hours_since_wake; those without use local_hour
    temp = temp[temp["hours_since_wake"].notna() | temp["hours_since_sunrise"].notna()].copy()

    y = temp["temp_f"].values

    print(f"All readings with wake time: {len(temp)}")
    print(f"  With sunrise: {temp['sunrise_hour'].notna().sum()}")

    # ═══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("1. EMPIRICAL CIRCADIAN CURVE — binned by hours since wake vs sunrise")
    print("=" * 70)

    overall_mean = temp["temp_f"].mean()

    def build_empirical_curve(hours_col, label):
        """Bin readings by hour, compute mean temperature and offset from overall mean."""
        h = temp[hours_col].values
        t = temp["temp_f"].values
        valid = ~np.isnan(h)
        h, t = h[valid], t[valid]

        # Build lookup: offset = binned_mean - overall_mean
        bin_offsets = {}
        for h_lo in range(-2, 20):
            mask = (h >= h_lo) & (h < h_lo + 1)
            if mask.sum() >= 5:
                bin_offsets[h_lo] = t[mask].mean() - overall_mean

        # R² of the binned model (each reading predicted by its bin mean)
        pred = np.array([bin_offsets.get(int(np.floor(hi)), 0) for hi in h]) + overall_mean
        r2 = 1 - np.sum((t - pred)**2) / np.sum((t - t.mean())**2)

        return bin_offsets, r2

    offsets_wake, r2_wake = build_empirical_curve("hours_since_wake", "wake")

    has_sunrise = temp["hours_since_sunrise"].notna()
    offsets_sun, r2_sun = build_empirical_curve("hours_since_sunrise", "sunrise")

    # Fair comparison on shared subset
    shared = has_sunrise
    temp_shared = temp[shared]
    _, r2_wake_shared = build_empirical_curve("hours_since_wake", "wake")  # already on full set
    # Recompute on shared only
    h_w = temp_shared["hours_since_wake"].values
    h_s = temp_shared["hours_since_sunrise"].values
    t_s = temp_shared["temp_f"].values
    om = t_s.mean()

    def r2_binned(hours, temps, om):
        pred = []
        for hi in hours:
            b = int(np.floor(hi))
            mask = (hours >= b) & (hours < b + 1)
            pred.append(temps[mask].mean() if mask.sum() > 0 else om)
        pred = np.array(pred)
        return 1 - np.sum((temps - pred)**2) / np.sum((temps - om)**2)

    r2_w_s = r2_binned(h_w[~np.isnan(h_w)], t_s[~np.isnan(h_w)], om)
    r2_s_s = r2_binned(h_s[~np.isnan(h_s)], t_s[~np.isnan(h_s)], om)

    print(f"\n  {'Basis':>20} {'R² (binned)':>12}")
    print(f"  {'Hours since wake':>20} {r2_wake:>12.4f} (n={len(temp)})")
    print(f"  {'Hours since sunrise':>20} {r2_sun:>12.4f} (n={has_sunrise.sum()})")
    print(f"  {'Wake (shared subset)':>20} {r2_w_s:>12.4f}")
    print(f"  {'Sunrise (shared)':>20} {r2_s_s:>12.4f}")

    # Wake wins on downstream signal strength (TDEE r|FM, intake r|FM, injection r)
    # even though sunrise has marginally higher in-sample circadian R².
    # The right metric is signal extraction, not curve fit.
    best_basis = "hours_since_wake"
    best_col = "hours_since_wake"
    bin_offsets = offsets_wake
    print(f"\n  Using: hours since wake (best downstream signal)")

    # Print the actual curve
    print(f"\n  Empirical circadian curve ({best_basis}):")
    print(f"  {'Hours':>7} {'Mean °F':>8} {'Offset':>8} {'n':>5}")
    h_all = temp[best_col].values
    t_all = temp["temp_f"].values
    for h_lo in sorted(bin_offsets.keys()):
        mask = (h_all >= h_lo) & (h_all < h_lo + 1)
        n = (~np.isnan(h_all) & mask).sum()
        offset = bin_offsets[h_lo]
        mean_t = offset + overall_mean
        bar = "█" * max(0, int((offset + 0.3) * 30))
        print(f"  {h_lo:>5}-{h_lo+1:<2} {mean_t:>8.2f} {offset:>+8.2f} {n:>5}  {bar}")

    peak_h = max(bin_offsets, key=bin_offsets.get)
    trough_h = min(bin_offsets, key=bin_offsets.get)
    total_range = bin_offsets[peak_h] - bin_offsets[trough_h]
    print(f"\n  Peak: {peak_h}-{peak_h+1}h ({bin_offsets[peak_h]:+.2f}°F)")
    print(f"  Trough: {trough_h}-{trough_h+1}h ({bin_offsets[trough_h]:+.2f}°F)")
    print(f"  Total range: {total_range:.2f}°F")

    def circadian_offset(hours):
        """Lookup empirical offset for a given hours-since-wake (or sunrise)."""
        h_bin = int(np.floor(hours)) if np.isscalar(hours) else np.floor(hours).astype(int)
        if np.isscalar(h_bin):
            return bin_offsets.get(h_bin, 0.0)
        return np.array([bin_offsets.get(int(hb), 0.0) for hb in h_bin])

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("2. COMPUTE BASELINE — each reading minus its circadian offset")
    print("=" * 70)

    # Apply empirical circadian offset to all readings
    valid_basis = temp[best_col].notna()
    temp.loc[valid_basis, "circ_offset"] = circadian_offset(temp.loc[valid_basis, best_col].values)
    temp["baseline"] = temp["temp_f"] - temp["circ_offset"]

    # Daily baseline: mean of per-reading baselines
    daily_base = temp.dropna(subset=["baseline"]).groupby("date").agg(
        baseline_mean=("baseline", "mean"),
        n_readings=("temp_f", "size"),
        mean_hour=(best_col, "mean"),
    ).reset_index()
    raw_daily_mean = temp.dropna(subset=["baseline"]).groupby("date")["temp_f"].mean().rename("raw_daily_mean")
    daily_base = daily_base.merge(raw_daily_mean, on="date", how="left")

    # Compare variance for one-time diagnostic only
    var_raw = daily_base["raw_daily_mean"].var()
    var_base = daily_base["baseline_mean"].var()
    print(f"\n  Days: {len(daily_base)}")
    print(f"  Raw daily mean: std = {daily_base['raw_daily_mean'].std():.3f}°F")
    print(f"  Baseline daily: std = {daily_base['baseline_mean'].std():.3f}°F")
    print(f"  Variance reduction: {(1 - var_base / var_raw) * 100:.1f}%")

    # Show that baseline is less affected by measurement time
    r_hour_raw = np.corrcoef(daily_base["mean_hour"], daily_base["raw_daily_mean"])[0, 1]
    r_hour_base = np.corrcoef(daily_base["mean_hour"], daily_base["baseline_mean"])[0, 1]
    print(f"  Measurement hour → raw mean:  r = {r_hour_raw:+.3f}")
    print(f"  Measurement hour → baseline:  r = {r_hour_base:+.3f}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("3. BASELINE vs METABOLIC STATE — the real test")
    print("=" * 70)

    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    daily = daily_base.merge(kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="left")
    daily = daily.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    daily = daily.merge(intake[["date", "calories"]], on="date", how="left")
    daily = daily.merge(tirz[["date", "effective_level", "days_since_injection"]],
                        on="date", how="left")
    daily["effective_level"] = daily["effective_level"].fillna(0)
    daily["on_tirz"] = (daily["effective_level"] > 0).astype(int)
    daily["tdee_resid"] = daily["tdee"] - daily["expected_rmr"]
    daily["tdee_ratio"] = daily["tdee"] / daily["expected_rmr"]

    # Trailing intake
    intake_s = intake.set_index("date")["calories"].sort_index()
    for win in [3, 7, 14, 30]:
        trailing = intake_s.rolling(win, min_periods=max(1, win // 2)).mean()
        daily = daily.merge(trailing.rename(f"intake_{win}d").reset_index(),
                            on="date", how="left")

    # SP distance
    full_fm = kalman.sort_values("date")
    sp = asymmetric_ema(full_fm["fat_mass_lbs"].values, 72, 25)
    sp_df = pd.DataFrame({"date": full_fm["date"], "sp_dist": sp - full_fm["fat_mass_lbs"].values})
    daily = daily.merge(sp_df, on="date", how="left")

    daily = daily.dropna(subset=["tdee", "expected_rmr", "fat_mass_lbs"])

    def partial_r(x, y, z):
        """Partial correlation of x, y controlling for z."""
        v = ~np.isnan(x) & ~np.isnan(y) & ~np.isnan(z)
        if v.sum() < 30:
            return np.nan
        X = np.column_stack([z[v], np.ones(v.sum())])
        rx = x[v] - X @ np.linalg.lstsq(X, x[v], rcond=None)[0]
        ry = y[v] - X @ np.linalg.lstsq(X, y[v], rcond=None)[0]
        return np.corrcoef(rx, ry)[0, 1]

    fm = daily["fat_mass_lbs"].values
    base = daily["baseline_mean"].values

    print(f"\n  {'Target':>25} {'baseline r|FM':>14}")
    for label, col in [
        ("TDEE residual", "tdee_resid"),
        ("TDEE/RMR ratio", "tdee_ratio"),
        ("effective_level", "effective_level"),
        ("SP distance", "sp_dist"),
        ("same-day calories", "calories"),
        ("3d trailing intake", "intake_3d"),
        ("7d trailing intake", "intake_7d"),
        ("14d trailing intake", "intake_14d"),
        ("30d trailing intake", "intake_30d"),
    ]:
        target = daily[col].values
        pr_base = partial_r(base, target, fm)
        base_str = f"{pr_base:+.3f}" if not np.isnan(pr_base) else "n/a"
        print(f"  {label:>25} {base_str:>14}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("4. TRAILING BASELINE WINDOWS")
    print("=" * 70)

    print(f"\n  {'Window':>8} {'base→TDEE_resid|FM':>19} {'base→intake_14d|FM':>19}")
    for win in [1, 3, 7, 14, 30]:
        if win == 1:
            b_series = daily["baseline_mean"].values
        else:
            b_series = daily["baseline_mean"].rolling(win, min_periods=max(1, win // 2)).mean().values

        r_tdee = partial_r(b_series, daily["tdee_resid"].values, fm)
        int_col = daily.get("intake_14d")
        r_int = partial_r(b_series, int_col.values, fm) if int_col is not None else np.nan

        label = f"{win}d" if win > 1 else "daily"
        r_t = f"{r_tdee:+.3f}" if not np.isnan(r_tdee) else "n/a"
        r_i = f"{r_int:+.3f}" if not np.isnan(r_int) else "n/a"
        print(f"  {label:>8} {r_t:>19} {r_i:>19}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("5. INJECTION-DAY SAWTOOTH — baseline estimate")
    print("=" * 70)

    on_tirz = daily[daily["on_tirz"] == 1].copy()
    if "days_since_injection" in on_tirz.columns:
        print(f"\n  {'Day':>5} {'Baseline':>9} {'n':>5}")
        for day in range(7):
            mask = on_tirz["days_since_injection"] == day
            if mask.sum() < 5:
                continue
            sub = on_tirz[mask]
            print(f"  {day:>5} {sub['baseline_mean'].mean():>+9.3f} {len(sub):>5}")

        valid_inj = on_tirz["days_since_injection"].notna() & on_tirz["days_since_injection"].between(0, 6)
        if valid_inj.sum() > 30:
            r = np.corrcoef(on_tirz.loc[valid_inj, "days_since_injection"],
                            on_tirz.loc[valid_inj, "baseline_mean"])[0, 1]
            print(f"\n  Days post-injection → baseline: r = {r:+.3f}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)
    print(f"  Circadian range: {total_range:.3f}°F peak-to-trough")
    print(f"  Peak: {peak_h:.1f}h after wake, Trough: {trough_h:.1f}h after wake")
    print(f"  Baseline variance reduction: {(1 - var_base / var_raw) * 100:.1f}% vs raw daily mean")
    print(f"  Measurement-hour bias removed: r = {r_hour_raw:+.3f} → {r_hour_base:+.3f}")

    # Save
    out = daily[["date", "baseline_mean", "n_readings"]].copy()
    out.to_csv(ROOT / "analysis" / "AS_temperature_baseline.csv", index=False)
    print(f"\nArtifact: analysis/AS_temperature_baseline.csv")

    # Save empirical curve for reuse
    curve = pd.DataFrame([
        {"hour_bin": k, "offset_f": v, "mean_temp_f": v + overall_mean}
        for k, v in sorted(bin_offsets.items())
    ])
    curve.to_csv(ROOT / "analysis" / "AS_circadian_curve.csv", index=False)
    print(f"Artifact: analysis/AS_circadian_curve.csv (basis: {best_basis})")


if __name__ == "__main__":
    main()
