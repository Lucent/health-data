#!/usr/bin/env python3
"""AR. Temperature circadian correction and metabolic correlates.

Previous approach: average to daily means, correlate with daily TDEE.
Problem: 2-3 readings/day averaged into one number buries signal in noise.

Better approach:
1. Model individual readings: temp ~ circadian(hour) + metabolic_state
2. Use trailing windows of circadian-corrected residuals (not daily means)
3. Exploit the injection-day sawtooth as a natural experiment
4. Test lagged intake → temperature (thermic effect of food)

Timezone: timestamps are in Eastern US time. Corrected to local using
sleep time_offset. Sleep-period readings excluded.
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


def main():
    np.random.seed(42)

    # ── Load all data ──
    temp = pd.read_csv(ROOT / "temperature" / "temperature.csv", parse_dates=["date"])
    temp = temp[(temp["temp_f"] > 95) & (temp["temp_f"] < 101)].copy()
    temp.rename(columns={"date": "timestamp"}, inplace=True)
    # Day boundary at 5am: readings before 5am belong to the previous day's waking period.
    temp["date"] = (temp["timestamp"] - pd.Timedelta(hours=5)).dt.floor("D")

    sleep = pd.read_csv(ROOT / "steps-sleep" / "sleep.csv", parse_dates=["date"])
    sleep["wake_hour"] = pd.to_timedelta(sleep["sleep_end"] + ":00").dt.total_seconds() / 3600
    sleep["sleep_start_hour"] = pd.to_timedelta(sleep["sleep_start"] + ":00").dt.total_seconds() / 3600
    sleep["local_offset_h"] = sleep["time_offset"].apply(parse_utc_offset_hours)

    sunlight = pd.read_csv(ROOT / "steps-sleep" / "sunlight.csv", parse_dates=["date"])
    sunlight["sunrise_hour"] = pd.to_timedelta(sunlight["sunrise"] + ":00").dt.total_seconds() / 3600

    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    # ── Timezone correction ──
    temp = temp.merge(sleep[["date", "wake_hour", "sleep_start_hour", "local_offset_h"]],
                      on="date", how="left")
    raw_hour = temp["timestamp"].dt.hour + temp["timestamp"].dt.minute / 60
    temp["home_offset_h"] = temp["timestamp"].apply(home_offset)
    temp["tz_correction"] = temp["local_offset_h"].fillna(temp["home_offset_h"]) - temp["home_offset_h"]
    temp["local_hour"] = (raw_hour + temp["tz_correction"]) % 24

    # ── Filter sleep-period readings ──
    awake = (temp["local_hour"] >= temp["wake_hour"]) & (
        (temp["sleep_start_hour"].isna()) |
        (temp["local_hour"] < temp["sleep_start_hour"]) |
        (temp["sleep_start_hour"] < temp["wake_hour"])
    )
    n_asleep = (~awake & temp["wake_hour"].notna()).sum()
    temp = temp[awake | temp["wake_hour"].isna()].copy()
    print(f"Readings: {len(temp)} awake ({n_asleep} sleep-period filtered)")

    # ── Merge metabolic state onto each reading ──
    temp = temp.merge(kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="left")
    temp = temp.merge(comp[["date", "expected_rmr", "ffm_lbs"]], on="date", how="left")
    temp = temp.merge(intake[["date", "calories"]], on="date", how="left")
    temp = temp.merge(tirz[["date", "effective_level", "days_since_injection"]],
                      on="date", how="left")
    temp["effective_level"] = temp["effective_level"].fillna(0)
    temp["on_tirz"] = (temp["effective_level"] > 0).astype(int)
    temp["tdee_resid"] = temp["tdee"] - temp["expected_rmr"]
    temp["tdee_ratio"] = temp["tdee"] / temp["expected_rmr"]

    # Trailing intake windows
    intake_s = intake.set_index("date")["calories"].sort_index()
    for win in [1, 3, 7, 14]:
        trailing = intake_s.rolling(win, min_periods=win).mean()
        trailing.name = f"intake_{win}d"
        temp = temp.merge(trailing.reset_index(), on="date", how="left")

    # Hours since wake
    temp["hours_since_wake"] = temp["local_hour"] - temp["wake_hour"]

    valid = (temp["tdee"].notna() & temp["expected_rmr"].notna() &
             temp["hours_since_wake"].notna() & temp["fat_mass_lbs"].notna())
    temp_v = temp[valid].copy()
    print(f"Readings with full metabolic data + wake time: {len(temp_v)}")

    # ═══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("1. JOINT MODEL — circadian + metabolic predictors on individual readings")
    print("=" * 70)

    # Build design matrix: temp ~ sin(h) + cos(h) + sin(h/12) + cos(h/12) + metabolic
    h = temp_v["hours_since_wake"].values
    y = temp_v["temp_f"].values

    # Circadian-only model
    X_circ = np.column_stack([
        np.sin(2 * np.pi * h / 24),
        np.cos(2 * np.pi * h / 24),
        np.sin(2 * np.pi * h / 12),
        np.cos(2 * np.pi * h / 12),
        np.ones(len(h)),
    ])
    beta_circ = np.linalg.lstsq(X_circ, y, rcond=None)[0]
    r2_circ = 1 - np.sum((y - X_circ @ beta_circ)**2) / np.sum((y - y.mean())**2)

    # Circadian + FM + TDEE residual
    X_met1 = np.column_stack([
        X_circ[:, :-1],  # sin/cos terms without intercept
        temp_v["fat_mass_lbs"].values,
        temp_v["tdee_resid"].values,
        np.ones(len(h)),
    ])
    beta_met1 = np.linalg.lstsq(X_met1, y, rcond=None)[0]
    r2_met1 = 1 - np.sum((y - X_met1 @ beta_met1)**2) / np.sum((y - y.mean())**2)

    # Circadian + FM + TDEE residual + effective_level
    X_met2 = np.column_stack([
        X_circ[:, :-1],
        temp_v["fat_mass_lbs"].values,
        temp_v["tdee_resid"].values,
        temp_v["effective_level"].values,
        np.ones(len(h)),
    ])
    beta_met2 = np.linalg.lstsq(X_met2, y, rcond=None)[0]
    r2_met2 = 1 - np.sum((y - X_met2 @ beta_met2)**2) / np.sum((y - y.mean())**2)

    # Circadian + FM + TDEE_ratio
    X_met3 = np.column_stack([
        X_circ[:, :-1],
        temp_v["fat_mass_lbs"].values,
        temp_v["tdee_ratio"].values,
        np.ones(len(h)),
    ])
    beta_met3 = np.linalg.lstsq(X_met3, y, rcond=None)[0]
    r2_met3 = 1 - np.sum((y - X_met3 @ beta_met3)**2) / np.sum((y - y.mean())**2)

    # Circadian + trailing intake
    for win in [1, 3, 7, 14]:
        col = f"intake_{win}d"
        valid_i = temp_v[col].notna()
        if valid_i.sum() < 200:
            continue
        sub = temp_v[valid_i]
        h_sub = sub["hours_since_wake"].values
        X_int = np.column_stack([
            np.sin(2 * np.pi * h_sub / 24),
            np.cos(2 * np.pi * h_sub / 24),
            np.sin(2 * np.pi * h_sub / 12),
            np.cos(2 * np.pi * h_sub / 12),
            sub["fat_mass_lbs"].values,
            sub[col].values,
            np.ones(len(h_sub)),
        ])
        y_sub = sub["temp_f"].values
        beta_int = np.linalg.lstsq(X_int, y_sub, rcond=None)[0]
        r2_int = 1 - np.sum((y_sub - X_int @ beta_int)**2) / np.sum((y_sub - y_sub.mean())**2)
        coef_intake = beta_int[5]
        # How many °F per 100 cal?
        print(f"  + {win}d intake:  R² = {r2_int:.4f} (Δ = +{r2_int - r2_circ:.4f}), "
              f"coef = {coef_intake * 100:+.3f}°F per 100 cal")

    print(f"\n  Circadian only:            R² = {r2_circ:.4f}")
    print(f"  + FM + TDEE residual:      R² = {r2_met1:.4f} (Δ = +{r2_met1 - r2_circ:.4f})")
    print(f"    TDEE_resid coef: {beta_met1[5]:+.5f}°F per cal")
    print(f"    FM coef:         {beta_met1[4]:+.5f}°F per lb")
    print(f"  + FM + TDEE resid + drug:  R² = {r2_met2:.4f} (Δ = +{r2_met2 - r2_circ:.4f})")
    print(f"    TDEE_resid coef: {beta_met2[5]:+.5f}°F per cal")
    print(f"    eff_level coef:  {beta_met2[6]:+.5f}°F per unit")
    print(f"  + FM + TDEE/RMR ratio:     R² = {r2_met3:.4f} (Δ = +{r2_met3 - r2_circ:.4f})")
    print(f"    TDEE_ratio coef: {beta_met3[5]:+.4f}°F per 0.01 ratio")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("2. TRAILING WINDOWS — circadian-corrected residual over multiple days")
    print("=" * 70)

    # Compute per-reading circadian residual
    temp_v["circ_resid"] = y - X_circ @ beta_circ

    # Aggregate to daily, then compute trailing windows
    daily = temp_v.groupby("date").agg(
        temp_resid=("circ_resid", "mean"),
        temp_raw=("temp_f", "mean"),
        n_readings=("temp_f", "size"),
    ).reset_index()

    daily = daily.merge(kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="left")
    daily = daily.merge(comp[["date", "expected_rmr"]], on="date", how="left")
    daily = daily.merge(intake[["date", "calories"]], on="date", how="left")
    daily = daily.merge(tirz[["date", "effective_level"]], on="date", how="left")
    daily["effective_level"] = daily["effective_level"].fillna(0)
    daily["tdee_resid"] = daily["tdee"] - daily["expected_rmr"]
    daily["tdee_ratio"] = daily["tdee"] / daily["expected_rmr"]
    daily = daily.dropna(subset=["tdee", "expected_rmr"])

    print(f"\n  {'Window':>8} {'temp→TDEE_resid':>16} {'temp→TDEE_ratio':>16} {'intake→temp':>16}")
    for win in [1, 3, 7, 14, 30]:
        if win == 1:
            t_series = daily["temp_resid"]
        else:
            t_series = daily["temp_resid"].rolling(win, min_periods=max(1, win // 2)).mean()

        tdee_r = daily["tdee_resid"]
        tdee_ratio = daily["tdee_ratio"]
        fm = daily["fat_mass_lbs"]

        valid_w = t_series.notna() & tdee_r.notna() & fm.notna()
        if valid_w.sum() < 50:
            continue

        # Partial r: temp → TDEE_resid | FM
        X = np.column_stack([fm[valid_w], np.ones(valid_w.sum())])
        res_t = t_series[valid_w].values - X @ np.linalg.lstsq(X, t_series[valid_w].values, rcond=None)[0]
        res_tdee = tdee_r[valid_w].values - X @ np.linalg.lstsq(X, tdee_r[valid_w].values, rcond=None)[0]
        res_ratio = tdee_ratio[valid_w].values - X @ np.linalg.lstsq(X, tdee_ratio[valid_w].values, rcond=None)[0]
        r_tdee = np.corrcoef(res_t, res_tdee)[0, 1]
        r_ratio = np.corrcoef(res_t, res_ratio)[0, 1]

        # Trailing intake → temp (partial | FM)
        cal_trail = daily["calories"].rolling(win, min_periods=max(1, win // 2)).mean()
        valid_c = valid_w & cal_trail.notna()
        if valid_c.sum() > 50:
            X2 = np.column_stack([fm[valid_c], np.ones(valid_c.sum())])
            res_t2 = t_series[valid_c].values - X2 @ np.linalg.lstsq(X2, t_series[valid_c].values, rcond=None)[0]
            res_c = cal_trail[valid_c].values - X2 @ np.linalg.lstsq(X2, cal_trail[valid_c].values, rcond=None)[0]
            r_intake = np.corrcoef(res_t2, res_c)[0, 1]
        else:
            r_intake = np.nan

        label = f"{win}d" if win > 1 else "daily"
        r_i_str = f"{r_intake:+.3f}" if not np.isnan(r_intake) else "n/a"
        print(f"  {label:>8} {r_tdee:+16.3f} {r_ratio:+16.3f} {r_i_str:>16}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("3. INJECTION-DAY SAWTOOTH — temperature by day post-injection")
    print("=" * 70)

    on_tirz = temp_v[temp_v["on_tirz"] == 1].copy()
    if "days_since_injection" in on_tirz.columns and on_tirz["days_since_injection"].notna().sum() > 50:
        print(f"\n  {'Day':>5} {'Mean °F':>8} {'Corrected':>10} {'n':>5} {'Mean intake':>12}")
        for day in range(8):
            mask = on_tirz["days_since_injection"] == day
            if mask.sum() < 10:
                continue
            sub = on_tirz[mask]
            print(f"  {day:>5} {sub['temp_f'].mean():>8.2f} {sub['circ_resid'].mean():>+10.3f} "
                  f"{len(sub):>5} {sub['calories'].mean():>12.0f}")

        # Correlation: injection day → circadian-corrected temp
        valid_inj = on_tirz["days_since_injection"].notna() & on_tirz["days_since_injection"].between(0, 6)
        if valid_inj.sum() > 50:
            r_inj = np.corrcoef(on_tirz.loc[valid_inj, "days_since_injection"],
                                on_tirz.loc[valid_inj, "circ_resid"])[0, 1]
            print(f"\n  Days post-injection → corrected temp: r = {r_inj:+.3f} (n={valid_inj.sum()})")

            # Same but controlling for same-day calories
            d_inj = on_tirz.loc[valid_inj, "days_since_injection"].values
            t_inj = on_tirz.loc[valid_inj, "circ_resid"].values
            c_inj = on_tirz.loc[valid_inj, "calories"].values
            valid_c = ~np.isnan(c_inj)
            if valid_c.sum() > 50:
                X = np.column_stack([c_inj[valid_c], np.ones(valid_c.sum())])
                res_d = d_inj[valid_c] - X @ np.linalg.lstsq(X, d_inj[valid_c], rcond=None)[0]
                res_t = t_inj[valid_c] - X @ np.linalg.lstsq(X, t_inj[valid_c], rcond=None)[0]
                r_inj_partial = np.corrcoef(res_d, res_t)[0, 1]
                print(f"  Partial (| same-day calories): r = {r_inj_partial:+.3f}")
    else:
        print("  Insufficient injection-day data")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("4. LAGGED INTAKE → TEMPERATURE — thermic effect timescale")
    print("=" * 70)

    # For each reading, compute trailing intake at various lags
    # This uses individual readings, not daily means
    print(f"\n  {'Lag':>8} {'r (raw)':>8} {'r|FM':>8} {'coef (°F/100cal)':>18}")
    for lag_name, col in [("same day", "intake_1d"), ("3d trail", "intake_3d"),
                          ("7d trail", "intake_7d"), ("14d trail", "intake_14d")]:
        valid_i = temp_v[col].notna() & temp_v["fat_mass_lbs"].notna()
        if valid_i.sum() < 200:
            continue
        sub = temp_v[valid_i]
        # Raw correlation with circadian residual
        r_raw = np.corrcoef(sub[col], sub["circ_resid"])[0, 1]
        # Partial | FM
        fm_s = sub["fat_mass_lbs"].values
        X = np.column_stack([fm_s, np.ones(len(fm_s))])
        res_t = sub["circ_resid"].values - X @ np.linalg.lstsq(X, sub["circ_resid"].values, rcond=None)[0]
        res_c = sub[col].values - X @ np.linalg.lstsq(X, sub[col].values, rcond=None)[0]
        r_partial = np.corrcoef(res_t, res_c)[0, 1]
        # Coefficient
        slope = np.polyfit(res_c, res_t, 1)[0]
        print(f"  {lag_name:>8} {r_raw:+8.3f} {r_partial:+8.3f} {slope * 100:+18.4f}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("5. SET POINT DISTANCE → TEMPERATURE")
    print("=" * 70)

    def asym_ema(fm_vals, hl_up, hl_down):
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

    full_fm = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    full_fm = full_fm.sort_values("date")
    sp = asym_ema(full_fm["fat_mass_lbs"].values, 72, 25)
    sp_df = pd.DataFrame({"date": full_fm["date"], "sp_dist": sp - full_fm["fat_mass_lbs"].values})

    temp_sp = temp_v.merge(sp_df, on="date", how="left")
    valid_sp = temp_sp["sp_dist"].notna() & temp_sp["fat_mass_lbs"].notna()
    if valid_sp.sum() > 200:
        sub = temp_sp[valid_sp]
        fm_s = sub["fat_mass_lbs"].values
        X = np.column_stack([fm_s, np.ones(len(fm_s))])
        res_t = sub["circ_resid"].values - X @ np.linalg.lstsq(X, sub["circ_resid"].values, rcond=None)[0]
        res_sp = sub["sp_dist"].values - X @ np.linalg.lstsq(X, sub["sp_dist"].values, rcond=None)[0]
        r_sp = np.corrcoef(res_t, res_sp)[0, 1]
        print(f"\n  SP distance → corrected temp | FM: r = {r_sp:+.3f} (n={valid_sp.sum()})")

        # Binned
        print(f"\n  {'SP distance':>15} {'Mean temp resid':>16} {'n':>5}")
        for lo, hi in [(-5, -2.5), (-2.5, 0), (0, 2.5), (2.5, 5), (5, 10)]:
            mask = (sub["sp_dist"] > lo) & (sub["sp_dist"] <= hi)
            if mask.sum() >= 10:
                print(f"  {lo:+.0f} to {hi:+.0f} lbs {sub.loc[mask, 'circ_resid'].mean():>+16.3f} {mask.sum():>5}")

    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)

    # Save
    out = daily[["date", "temp_resid", "temp_raw", "n_readings"]].copy()
    out.to_csv(ROOT / "analysis" / "AR_temperature_corrected.csv", index=False)
    print(f"\nArtifact: analysis/AR_temperature_corrected.csv")


if __name__ == "__main__":
    main()
