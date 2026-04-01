#!/usr/bin/env python3
"""AZ. Exhaustive search for a set point model that fits both this subject
and SURMOUNT-4 regain.

Search space:
  SP functional forms: EMA, rolling mean, linear approach, log, two-phase
  Outcome measures: daily surplus, 30/60/90/120d rolling surplus, FM velocity
  Pressure functions: linear, sqrt, quadratic, sigmoid

Targets:
  1. This subject 2014+ pre-tirz: r(SP_distance, outcome)
  2. SURMOUNT-4 regain: RMSE of predicted vs actual surplus at each timepoint
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

# SURMOUNT-4 regain
S4_WEEKS = np.array([0, 4, 8, 12, 16, 20, 24, 32, 40, 52])
S4_PCT = np.array([0, 2.5, 5.0, 7.0, 8.5, 9.8, 11.0, 12.5, 13.3, 14.0])
S4_W36_LBS = 180.0
S4_FM36 = 48.0
S4_SP_ORIGINAL = 90.0
S4_FM = S4_FM36 + (S4_W36_LBS * S4_PCT / 100) * 0.93

# Actual regain surplus at each interval (cal/day)
S4_DAYS = S4_WEEKS * 7
S4_FM_RATES = np.diff(S4_FM) / np.diff(S4_DAYS)  # lbs/day
S4_SURPLUS = S4_FM_RATES * 3500 / 0.93  # cal/day in each interval


def load_subject():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])

    df = intake[["date", "calories"]].merge(
        kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="inner")
    df = df.merge(tirz[["date", "effective_level"]], on="date", how="left")
    df["effective_level"] = df["effective_level"].fillna(0)
    df = df.sort_values("date").reset_index(drop=True)
    df["surplus"] = df["calories"] - df["tdee"]
    post = df[(df["date"] >= "2014-01-01") & (df["effective_level"] == 0)].copy()
    return post.reset_index(drop=True)


# ═══════════════════════════════════════════════════════════════════
# SP FUNCTIONAL FORMS
# ═══════════════════════════════════════════════════════════════════

def sp_ema(fm, hl):
    """Exponential moving average."""
    alpha = 1 - np.exp(-np.log(2) / hl)
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        sp[i] = sp[i - 1] + alpha * (fm[i] - sp[i - 1])
    return sp


def sp_rolling_mean(fm, window):
    """Rolling mean of FM."""
    sp = pd.Series(fm).rolling(int(window), min_periods=1, center=False).mean().values
    return sp


def sp_rolling_median(fm, window):
    """Rolling median of FM."""
    sp = pd.Series(fm).rolling(int(window), min_periods=1, center=False).median().values
    return sp


def sp_linear_approach(fm, rate_per_day):
    """SP moves toward FM at a fixed rate (lbs/day), never overshooting."""
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    for i in range(1, len(fm)):
        diff = fm[i] - sp[i - 1]
        step = np.sign(diff) * min(abs(diff), rate_per_day)
        sp[i] = sp[i - 1] + step
    return sp


def sp_ema_with_floor(fm, hl, floor_fraction):
    """EMA that can't adapt below floor_fraction of its all-time max."""
    alpha = 1 - np.exp(-np.log(2) / hl)
    sp = np.empty(len(fm))
    sp[0] = fm[0]
    peak = fm[0]
    for i in range(1, len(fm)):
        sp[i] = sp[i - 1] + alpha * (fm[i] - sp[i - 1])
        peak = max(peak, sp[i])
        sp[i] = max(sp[i], peak * floor_fraction)
    return sp


def sp_two_phase(fm, hl_fast, hl_slow, fast_weight):
    """Two-phase: SP = w*EMA(fast) + (1-w)*EMA(slow)."""
    sp_f = sp_ema(fm, hl_fast)
    sp_s = sp_ema(fm, hl_slow)
    return fast_weight * sp_f + (1 - fast_weight) * sp_s


# ═══════════════════════════════════════════════════════════════════
# PRESSURE FUNCTIONS
# ═══════════════════════════════════════════════════════════════════

def pressure_linear(gap, cal_per_lb):
    return cal_per_lb * gap


def pressure_sqrt(gap, cal_per_sqrtlb):
    return cal_per_sqrtlb * np.sign(gap) * np.sqrt(np.abs(gap))


def pressure_quadratic(gap, cal_per_lb2):
    return cal_per_lb2 * np.sign(gap) * gap ** 2


# ═══════════════════════════════════════════════════════════════════
# EVALUATION
# ═══════════════════════════════════════════════════════════════════

def eval_subject(fm, sp, surplus, windows=[1, 30, 60, 90, 120]):
    """Best r(distance, surplus) across outcome windows."""
    dist = sp - fm
    best_r = 0
    best_win = 90
    for win in windows:
        if win == 1:
            target = surplus
        else:
            target = pd.Series(surplus).rolling(win, min_periods=win).mean().values
        v = ~np.isnan(dist) & ~np.isnan(target)
        if v.sum() < 200:
            continue
        r = np.corrcoef(dist[v], target[v])[0, 1]
        if abs(r) > abs(best_r):
            best_r = r
            best_win = win
    return best_r, best_win


def eval_surmount(sp_func, sp_params, pressure_func, pressure_params):
    """Predict SURMOUNT-4 regain and compute RMSE against actual surplus."""
    # Treatment phase: FM drops from 90 to 48 over 252 days
    fm_treatment = np.linspace(S4_SP_ORIGINAL, S4_FM36, 252)
    sp_treat = sp_func(fm_treatment, *sp_params)
    sp_at_stop = sp_treat[-1]

    # Regain phase: FM rises per published data
    fm_daily = np.interp(np.arange(S4_DAYS[-1] + 1), S4_DAYS, S4_FM)
    sp_regain = np.empty(len(fm_daily))
    sp_regain[0] = sp_at_stop

    # Propagate SP during regain
    if sp_func == sp_ema:
        hl = sp_params[0]
        alpha = 1 - np.exp(-np.log(2) / hl)
        for d in range(1, len(fm_daily)):
            sp_regain[d] = sp_regain[d - 1] + alpha * (fm_daily[d] - sp_regain[d - 1])
    elif sp_func == sp_linear_approach:
        rate = sp_params[0]
        for d in range(1, len(fm_daily)):
            diff = fm_daily[d] - sp_regain[d - 1]
            step = np.sign(diff) * min(abs(diff), rate)
            sp_regain[d] = sp_regain[d - 1] + step
    elif sp_func == sp_ema_with_floor:
        hl, floor_frac = sp_params
        alpha = 1 - np.exp(-np.log(2) / hl)
        peak = sp_at_stop
        for d in range(1, len(fm_daily)):
            sp_regain[d] = sp_regain[d - 1] + alpha * (fm_daily[d] - sp_regain[d - 1])
            peak = max(peak, sp_regain[d])
            sp_regain[d] = max(sp_regain[d], peak * floor_frac)
    elif sp_func == sp_two_phase:
        hl_f, hl_s, w = sp_params
        alpha_f = 1 - np.exp(-np.log(2) / hl_f)
        alpha_s = 1 - np.exp(-np.log(2) / hl_s)
        # Need to decompose sp_at_stop into fast/slow components
        sp_f_treat = sp_ema(fm_treatment, hl_f)
        sp_s_treat = sp_ema(fm_treatment, hl_s)
        sp_f = np.empty(len(fm_daily))
        sp_s = np.empty(len(fm_daily))
        sp_f[0] = sp_f_treat[-1]
        sp_s[0] = sp_s_treat[-1]
        for d in range(1, len(fm_daily)):
            sp_f[d] = sp_f[d - 1] + alpha_f * (fm_daily[d] - sp_f[d - 1])
            sp_s[d] = sp_s[d - 1] + alpha_s * (fm_daily[d] - sp_s[d - 1])
        sp_regain = w * sp_f + (1 - w) * sp_s
    else:
        return np.inf  # unsupported for forward propagation

    # Predicted surplus at midpoints of each interval
    pred_surplus = []
    for i in range(len(S4_DAYS) - 1):
        d_mid = (S4_DAYS[i] + S4_DAYS[i + 1]) // 2
        gap = sp_regain[d_mid] - fm_daily[d_mid]
        pred_surplus.append(pressure_func(gap, *pressure_params))

    pred_surplus = np.array(pred_surplus)
    rmse = np.sqrt(np.mean((pred_surplus - S4_SURPLUS) ** 2))
    return rmse


def main():
    post = load_subject()
    fm = post["fat_mass_lbs"].values.copy()
    surplus = post["surplus"].values.copy()

    # Fill NaN
    first_valid = np.where(~np.isnan(fm))[0][0]
    fm[:first_valid] = fm[first_valid]

    # Fixed outcome: 90d rolling surplus (the known best for this subject)
    surplus_90 = pd.Series(surplus).rolling(90, min_periods=90).mean().values

    results = []

    def eval_r90(fm_arr, sp_arr):
        """r(distance, 90d surplus) — the benchmark metric."""
        dist = sp_arr - fm_arr
        v = ~np.isnan(dist) & ~np.isnan(surplus_90)
        if v.sum() < 200:
            return 0
        return np.corrcoef(dist[v], surplus_90[v])[0, 1]

    # ═══════════════════════════════════════════════════════════════
    print("=" * 70)
    print("SEARCH: what SP functions match this subject AND SURMOUNT-4 regain?")
    print("Outcome fixed: 90-day rolling surplus (known best)")
    print("=" * 70)

    # 1. EMA sweep
    for hl in [30, 40, 45, 50, 60, 80, 100, 120, 150, 200, 250, 300]:
        sp = sp_ema(fm, hl)
        r = eval_r90(fm, sp)
        for cal_lb in [15, 20, 27, 35, 45, 60, 80]:
            rmse = eval_surmount(sp_ema, (hl,), pressure_linear, (cal_lb,))
            results.append(("EMA", f"hl={hl}", "linear", f"{cal_lb}cal/lb",
                            r, 90, rmse))

    # 2. Rolling mean sweep
    for win_sp in [60, 90, 120, 180, 270, 365]:
        sp = sp_rolling_mean(fm, win_sp)
        r = eval_r90(fm, sp)
        for cal_lb in [20, 27, 35, 45]:
            rmse = eval_surmount(sp_rolling_mean, (win_sp,), pressure_linear, (cal_lb,))
            results.append(("RollingMean", f"win={win_sp}", "linear", f"{cal_lb}cal/lb",
                            r, 90, rmse))

    # 3. Linear approach sweep
    for rate in [0.05, 0.1, 0.2, 0.3, 0.5, 1.0]:
        sp = sp_linear_approach(fm, rate)
        r = eval_r90(fm, sp)
        for cal_lb in [20, 27, 35, 45]:
            rmse = eval_surmount(sp_linear_approach, (rate,), pressure_linear, (cal_lb,))
            results.append(("LinearAppr", f"rate={rate}", "linear", f"{cal_lb}cal/lb",
                            r, 90, rmse))

    # 4. EMA with floor
    for hl in [45, 80, 120]:
        for floor_f in [0.5, 0.6, 0.7, 0.8]:
            sp = sp_ema_with_floor(fm, hl, floor_f)
            r = eval_r90(fm, sp)
            rmse = eval_surmount(sp_ema_with_floor, (hl, floor_f), pressure_linear, (27,))
            results.append(("EMA+floor", f"hl={hl},f={floor_f}", "linear", "27cal/lb",
                            r, 90, rmse))

    # 5. Two-phase
    for hl_f in [30, 45, 60]:
        for hl_s in [150, 250, 400]:
            for w in [0.3, 0.5, 0.7, 0.9]:
                sp = sp_two_phase(fm, hl_f, hl_s, w)
                r = eval_r90(fm, sp)
                for cal_lb in [27, 45, 60]:
                    rmse = eval_surmount(sp_two_phase, (hl_f, hl_s, w),
                                         pressure_linear, (cal_lb,))
                    results.append(("TwoPhase", f"f={hl_f},s={hl_s},w={w}",
                                    "linear", f"{cal_lb}cal/lb", r, 90, rmse))

    # 6. EMA with sqrt pressure
    for hl in [45, 80, 120, 200]:
        sp = sp_ema(fm, hl)
        r = eval_r90(fm, sp)
        for cal_sqrt in [30, 50, 80, 120]:
            rmse = eval_surmount(sp_ema, (hl,), pressure_sqrt, (cal_sqrt,))
            results.append(("EMA", f"hl={hl}", "sqrt", f"{cal_sqrt}cal/√lb",
                            r, 90, rmse))

    # Convert to DataFrame and rank
    rdf = pd.DataFrame(results, columns=["SP_form", "SP_params", "Pressure", "P_params",
                                          "r_subject", "best_window", "regain_RMSE"])
    rdf["abs_r"] = rdf["r_subject"].abs()

    # Pareto front: best tradeoffs between subject r and regain RMSE
    print(f"\n  Total configurations tested: {len(rdf)}")

    # Top 10 by combined score (normalize both to 0-1, then add)
    rdf["r_norm"] = (rdf["abs_r"] - rdf["abs_r"].min()) / (rdf["abs_r"].max() - rdf["abs_r"].min())
    rdf["rmse_norm"] = 1 - (rdf["regain_RMSE"] - rdf["regain_RMSE"].min()) / (rdf["regain_RMSE"].max() - rdf["regain_RMSE"].min())
    rdf["score"] = rdf["r_norm"] + rdf["rmse_norm"]

    print(f"\n  Top 15 by combined score (subject r + regain fit):")
    print(f"  {'SP form':>12} {'SP params':>25} {'Pressure':>10} {'P params':>12} "
          f"{'r':>7} {'win':>4} {'RMSE':>7} {'score':>6}")
    top = rdf.nlargest(15, "score")
    for _, row in top.iterrows():
        print(f"  {row['SP_form']:>12} {row['SP_params']:>25} {row['Pressure']:>10} {row['P_params']:>12} "
              f"{row['r_subject']:>+7.3f} {row['best_window']:>4}d {row['regain_RMSE']:>7.0f} {row['score']:>6.3f}")

    # Best subject-only
    print(f"\n  Best subject r (ignoring regain):")
    best_subj = rdf.nlargest(5, "abs_r")
    for _, row in best_subj.iterrows():
        print(f"  {row['SP_form']:>12} {row['SP_params']:>25} r={row['r_subject']:+.3f} RMSE={row['regain_RMSE']:.0f}")

    # Best regain-only
    print(f"\n  Best regain RMSE (ignoring subject):")
    best_regain = rdf.nsmallest(5, "regain_RMSE")
    for _, row in best_regain.iterrows():
        print(f"  {row['SP_form']:>12} {row['SP_params']:>25} {row['P_params']:>12} r={row['r_subject']:+.3f} RMSE={row['regain_RMSE']:.0f}")


if __name__ == "__main__":
    main()
