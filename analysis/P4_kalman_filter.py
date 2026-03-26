"""Kalman filter for fat mass and TDEE estimation.

State: [fat_mass_lbs, tdee_cal]
Process: fat_mass(t+1) = fat_mass(t) + (intake(t) - tdee(t)) / 3500
         tdee(t+1) = tdee(t) + noise  (random walk)
Observation: smoothed_weight(t) - lean_mass(t) = fat_mass(t) + noise
             (only on days with weigh-ins)

Includes Rauch-Tung-Striebel backward smoother for optimal estimates.

Inputs:
    analysis/smoothed_weight.csv — glycogen+sodium corrected weight (sparse)
    analysis/daily_composition.csv — daily lean mass + expected RMR
    intake/intake_daily.csv — daily calories

Outputs:
    analysis/kalman_daily.csv — daily fat mass, TDEE, uncertainties
    analysis/plot_kalman.png — diagnostic plots
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

CAL_PER_LB = 3500

# Noise parameters (tuned below)
R_DEFAULT = 0.97       # observation noise variance (lbs²)
Q_FAT_DEFAULT = 0.0025 # fat mass process noise (lbs²/day)
Q_TDEE_DEFAULT = 25    # TDEE process noise (cal²/day)

CALORIMETRY = {
    "2011-04-30": 2415,
    "2012-05-03": 1956,
    "2016-01-30": 1700,
}


def load_data():
    smooth = pd.read_csv(ROOT / "analysis" / "P1_smoothed_weight.csv", parse_dates=["date"])
    comp = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    return smooth, comp, intake


def build_inputs(smooth, comp, intake):
    """Build aligned arrays for the filter."""
    # Base: all intake days
    dates = intake["date"].values
    calories = intake["calories"].values
    n = len(dates)

    # Lean mass from composition
    lean = comp.set_index("date")["ffm_lbs"].reindex(intake["date"]).values

    # Observations: smoothed weight on days with weigh-ins
    obs_map = smooth.set_index("date")["smoothed_weight_lbs"]
    observations = np.full(n, np.nan)
    for i, d in enumerate(dates):
        if d in obs_map.index:
            v = obs_map[d]
            if not np.isnan(v) and not np.isnan(lean[i]):
                observations[i] = v - lean[i]  # fat mass observation

    # Expected RMR for initialization
    expected_rmr = comp.set_index("date")["expected_rmr"].reindex(intake["date"]).values

    return dates, calories, observations, lean, expected_rmr


def kalman_forward(calories, observations, expected_rmr, q_tdee,
                   q_fat=Q_FAT_DEFAULT, r=R_DEFAULT, mean_revert=0.02):
    """Forward Kalman filter pass with mean-reverting TDEE.

    TDEE process: tdee(t+1) = tdee(t) + alpha * (target(t) - tdee(t)) + noise
    where target = expected_rmr * 1.2. This prevents unbounded drift during
    long gaps while still allowing the filter to discover TDEE from observations.

    alpha=0: pure random walk. alpha=1: snaps to expected each day.
    Default 0.02: pulls ~2% toward expected per day → 50% correction over 35 days.

    Returns:
        x_filt, P_filt, x_pred, P_pred, innovations
    """
    n = len(calories)

    H = np.array([[1.0, 0.0]])
    Q = np.array([[q_fat, 0.0],
                  [0.0,   q_tdee]])
    R_mat = np.array([[r]])

    first_obs_idx = np.where(~np.isnan(observations))[0]
    if len(first_obs_idx) == 0:
        raise ValueError("No observations")
    i0 = first_obs_idx[0]

    init_tdee = expected_rmr[i0] * 1.2 if not np.isnan(expected_rmr[i0]) else 2100.0
    x = np.array([observations[i0], init_tdee])
    P = np.array([[25.0,     0.0],
                  [0.0,  250000.0]])

    x_filt = np.full((n, 2), np.nan)
    P_filt = np.full((n, 2, 2), np.nan)
    x_pred_arr = np.full((n, 2), np.nan)
    P_pred_arr = np.full((n, 2, 2), np.nan)
    innovations = np.full(n, np.nan)

    x_filt[i0] = x
    P_filt[i0] = P
    x_pred_arr[i0] = x
    P_pred_arr[i0] = P

    for i in range(i0, n - 1):
        # TDEE target: composition-aware expected RMR × activity factor
        target_tdee = expected_rmr[i] * 1.2 if not np.isnan(expected_rmr[i]) else x[1]

        # Process model with mean-reverting TDEE
        # fat(t+1) = fat(t) + (intake - tdee) / 3500
        # tdee(t+1) = tdee(t) + alpha * (target - tdee(t))
        F = np.array([[1.0, -1.0 / CAL_PER_LB],
                      [0.0,  1.0 - mean_revert]])
        B_vec = np.array([calories[i] / CAL_PER_LB,
                          mean_revert * target_tdee])

        x_pred = F @ x + B_vec
        P_pred = F @ P @ F.T + Q

        x_pred_arr[i + 1] = x_pred
        P_pred_arr[i + 1] = P_pred

        if not np.isnan(observations[i + 1]):
            y = observations[i + 1]
            y_pred = (H @ x_pred)[0]
            innovation = y - y_pred
            S = (H @ P_pred @ H.T + R_mat)[0, 0]
            K = (P_pred @ H.T) / S
            x = x_pred + (K * innovation).flatten()
            P = P_pred - np.outer(K.flatten(), H @ P_pred)
            innovations[i + 1] = innovation
        else:
            x = x_pred
            P = P_pred

        x_filt[i + 1] = x
        P_filt[i + 1] = P

    return x_filt, P_filt, x_pred_arr, P_pred_arr, innovations


def rts_smoother(x_filt, P_filt, x_pred, P_pred, calories, mean_revert=0.02):
    """Rauch-Tung-Striebel backward smoother.

    Produces optimal estimates using all observations (past and future).
    """
    n = len(x_filt)
    F = np.array([[1.0, -1.0 / CAL_PER_LB],
                  [0.0,  1.0 - mean_revert]])

    x_smooth = np.copy(x_filt)
    P_smooth = np.copy(P_filt)

    # Find last valid index
    valid = ~np.isnan(x_filt[:, 0])
    last_valid = np.where(valid)[0][-1]

    for i in range(last_valid - 1, -1, -1):
        if np.isnan(x_filt[i, 0]) or np.isnan(P_pred[i + 1, 0, 0]):
            continue

        P_pred_inv = np.linalg.inv(P_pred[i + 1])
        G = P_filt[i] @ F.T @ P_pred_inv  # smoother gain

        x_smooth[i] = x_filt[i] + G @ (x_smooth[i + 1] - x_pred[i + 1])
        P_smooth[i] = P_filt[i] + G @ (P_smooth[i + 1] - P_pred[i + 1]) @ G.T

    return x_smooth, P_smooth


def tune_params(calories, observations, expected_rmr):
    """Sweep Q_tdee and mean_revert, select by innovation whiteness + calorimetry."""
    print("\n=== Tuning Q_tdee and mean_revert ===")
    print(f"{'Q':>5} {'alpha':>6} {'acf1':>7} {'innov_std':>10} {'tdee_dd':>8} {'cal_err':>8} {'score':>7}")

    dates_arr = pd.read_csv(ROOT / "intake" / "intake_daily.csv")["date"].values
    best = None
    best_score = 999
    best_params = (25, 0.02)

    for q_tdee in [10, 25, 50, 100, 200, 500]:
        for mr in [0.005, 0.01, 0.02, 0.05, 0.10]:
            x_f, P_f, x_p, P_p, innov = kalman_forward(
                calories, observations, expected_rmr, q_tdee, mean_revert=mr)

            valid_innov = innov[~np.isnan(innov)]
            if len(valid_innov) < 50:
                continue
            innov_std = np.std(valid_innov)

            centered = valid_innov - valid_innov.mean()
            acf1 = np.corrcoef(centered[:-1], centered[1:])[0, 1]

            valid_tdee = x_f[~np.isnan(x_f[:, 1]), 1]
            tdee_dd = np.std(np.diff(valid_tdee))

            cal_errors = []
            for date_str, rmr in CALORIMETRY.items():
                idx = np.searchsorted(dates_arr, date_str)
                if idx < len(x_f) and not np.isnan(x_f[idx, 1]):
                    cal_errors.append(abs(x_f[idx, 1] / rmr - 1.2))
            cal_err = np.mean(cal_errors) if cal_errors else 999

            score = abs(acf1) + cal_err * 2
            print(f"{q_tdee:5.0f} {mr:6.3f} {acf1:7.4f} {innov_std:10.3f} "
                  f"{tdee_dd:8.1f} {cal_err:8.4f} {score:7.4f}")

            if score < best_score:
                best_score = score
                best_params = (q_tdee, mr)

    q, mr = best_params
    print(f"\nSelected: Q_tdee={q}, mean_revert={mr}")
    return q, mr


def validate(dates, x_smooth, P_smooth, innovations, observations, lean, expected_rmr):
    """Validation report."""
    print("\n=== Validation ===\n")

    # Innovation stats
    valid_innov = innovations[~np.isnan(innovations)]
    centered = valid_innov - valid_innov.mean()
    acf1 = np.corrcoef(centered[:-1], centered[1:])[0, 1]
    print(f"Innovations: n={len(valid_innov)}  mean={valid_innov.mean():.3f}  "
          f"std={valid_innov.std():.3f}  lag1_acf={acf1:.4f}")

    # Fat mass at composition dates
    comp = pd.read_csv(ROOT / "composition" / "composition.csv", parse_dates=["date"])
    intake_dates = pd.read_csv(ROOT / "intake" / "intake_daily.csv")["date"].values
    print(f"\nFat mass at composition measurement dates:")
    print(f"{'Date':>12} {'Meas FM':>8} {'Filter FM':>10} {'±':>5} {'Error':>6}")
    for _, row in comp.iterrows():
        date_str = row["date"].strftime("%Y-%m-%d")
        idx = np.searchsorted(intake_dates, date_str)
        if idx < len(x_smooth) and not np.isnan(x_smooth[idx, 0]):
            fm_filter = x_smooth[idx, 0]
            fm_std = np.sqrt(P_smooth[idx, 0, 0])
            fm_meas = row["fat_mass_lbs"]
            err = fm_filter - fm_meas
            print(f"{date_str:>12} {fm_meas:8.1f} {fm_filter:10.1f} {fm_std:5.1f} {err:+6.1f}")

    # Calorimetry
    print(f"\nCalorimetry anchors:")
    for date_str, rmr in CALORIMETRY.items():
        idx = np.searchsorted(intake_dates, date_str)
        if idx < len(x_smooth) and not np.isnan(x_smooth[idx, 1]):
            tdee = x_smooth[idx, 1]
            tdee_std = np.sqrt(P_smooth[idx, 1, 1])
            print(f"  {date_str}: RMR={rmr}  TDEE={tdee:.0f}±{tdee_std:.0f}  "
                  f"ratio={tdee/rmr:.2f}")

    # TDEE smoothness vs window method
    valid_tdee = x_smooth[~np.isnan(x_smooth[:, 1]), 1]
    tdee_diff_std = np.std(np.diff(valid_tdee))
    print(f"\nTDEE day-to-day std: {tdee_diff_std:.1f} cal/day "
          f"(window method: ~640 cal/day)")

    # Uncertainty bands
    valid_std = np.sqrt(P_smooth[~np.isnan(P_smooth[:, 0, 0]), 0, 0])
    print(f"\nFat mass uncertainty (std):")
    print(f"  Median: {np.median(valid_std):.2f} lbs")
    print(f"  P5:     {np.percentile(valid_std, 5):.2f} lbs")
    print(f"  P95:    {np.percentile(valid_std, 95):.2f} lbs")
    print(f"  Max:    {np.max(valid_std):.2f} lbs (during longest gap)")


def save_output(dates, x_filt, P_filt, x_smooth, P_smooth, innovations):
    """Save daily Kalman output.

    The filtered columns are causal estimates available at that date.
    The smoothed columns are retrospective estimates that use future weigh-ins.
    """
    n = len(dates)
    out = pd.DataFrame({
        "date": dates,
        "fat_mass_lbs_filtered": np.round(x_filt[:, 0], 2),
        "fat_mass_std_filtered": np.round(np.sqrt(np.maximum(P_filt[:, 0, 0], 0)), 2),
        "tdee_filtered": np.round(x_filt[:, 1], 0),
        "tdee_std_filtered": np.round(np.sqrt(np.maximum(P_filt[:, 1, 1], 0)), 0),
        "fat_mass_lbs": np.round(x_smooth[:, 0], 2),
        "fat_mass_std": np.round(np.sqrt(np.maximum(P_smooth[:, 0, 0], 0)), 2),
        "tdee": np.round(x_smooth[:, 1], 0),
        "tdee_std": np.round(np.sqrt(np.maximum(P_smooth[:, 1, 1], 0)), 0),
        "innovation": np.round(innovations, 3),
    })
    out["date"] = pd.to_datetime(out["date"]).dt.strftime("%Y-%m-%d")
    path = ROOT / "analysis" / "P4_kalman_daily.csv"
    out.to_csv(path, index=False)
    print(f"\nWrote {len(out)} rows to {path}")


def plot_kalman(dates, x_smooth, P_smooth, innovations, observations, lean,
                expected_rmr, calories):
    """Diagnostic plots."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    dates_dt = pd.to_datetime(dates)
    fat = x_smooth[:, 0]
    tdee = x_smooth[:, 1]
    fat_std = np.sqrt(np.maximum(P_smooth[:, 0, 0], 0))
    tdee_std = np.sqrt(np.maximum(P_smooth[:, 1, 1], 0))

    fig, axes = plt.subplots(4, 1, figsize=(18, 16), sharex=True)

    # Panel 1: Fat mass with uncertainty
    ax = axes[0]
    obs_mask = ~np.isnan(observations)
    ax.scatter(dates_dt[obs_mask], observations[obs_mask], s=2, color="tab:blue",
               alpha=0.3, label="Observed (corrected wt - lean)", zorder=3)
    ax.plot(dates_dt, fat, color="tab:red", lw=1, label="Kalman fat mass")
    ax.fill_between(dates_dt, fat - 2 * fat_std, fat + 2 * fat_std,
                    alpha=0.15, color="tab:red", label="±2σ")
    ax.set_ylabel("Fat mass (lbs)")
    ax.legend(fontsize=8)
    ax.set_title("Kalman Filter: Fat Mass Estimate with Uncertainty")
    ax.grid(True, alpha=0.2)

    # Panel 2: TDEE with uncertainty
    ax = axes[1]
    ax.plot(dates_dt, tdee, color="tab:purple", lw=1, label="Kalman TDEE")
    ax.fill_between(dates_dt, tdee - 2 * tdee_std, tdee + 2 * tdee_std,
                    alpha=0.15, color="tab:purple", label="±2σ")

    # 30-day rolling intake
    intake_30 = pd.Series(calories).rolling(30, center=True, min_periods=10).mean()
    ax.plot(dates_dt, intake_30, color="tab:green", lw=0.8, alpha=0.5,
            label="30-day intake")

    # Expected RMR
    rmr_30 = pd.Series(expected_rmr).rolling(30, center=True, min_periods=10).mean()
    ax.plot(dates_dt, rmr_30, color="tab:orange", lw=0.8, ls="--",
            label="Expected RMR (composition)")

    # Calorimetry
    for date_str, rmr in CALORIMETRY.items():
        ax.plot(pd.Timestamp(date_str), rmr, "D", color="black", ms=8, zorder=5)

    ax.set_ylabel("Cal/day")
    ax.set_ylim(1000, 3500)
    ax.legend(fontsize=8, ncol=3)
    ax.set_title("Kalman Filter: TDEE Estimate")
    ax.grid(True, alpha=0.2)

    # Panel 3: Innovations
    ax = axes[2]
    valid_mask = ~np.isnan(innovations)
    ax.scatter(dates_dt[valid_mask], innovations[valid_mask], s=1, color="tab:purple",
               alpha=0.3)
    innov_30 = pd.Series(innovations).rolling(30, center=True, min_periods=5).mean()
    ax.plot(dates_dt, innov_30, color="tab:purple", lw=1.5, label="30-day avg")
    ax.axhline(y=0, color="gray", ls="--")
    ax.set_ylabel("Innovation (lbs)")
    ax.legend(fontsize=8)
    ax.set_title("Innovation Sequence (should be white noise around zero)")
    ax.grid(True, alpha=0.2)

    # Panel 4: Uncertainty over time
    ax = axes[3]
    ax.plot(dates_dt, fat_std, color="tab:red", lw=1, label="Fat mass σ")
    ax2 = ax.twinx()
    ax2.plot(dates_dt, tdee_std, color="tab:purple", lw=1, label="TDEE σ")
    ax.set_ylabel("Fat mass σ (lbs)", color="tab:red")
    ax2.set_ylabel("TDEE σ (cal)", color="tab:purple")
    ax.set_title("Uncertainty: Wide during gaps, narrow during dense weighing")
    ax.grid(True, alpha=0.2)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

    plt.tight_layout()
    out = ROOT / "analysis" / "P4_plot_kalman.png"
    plt.savefig(out, dpi=150)
    print(f"Saved {out}")
    plt.close()


def main():
    print("Loading data...")
    smooth, comp, intake = load_data()

    print("Building filter inputs...")
    dates, calories, observations, lean, expected_rmr = build_inputs(smooth, comp, intake)
    n_obs = (~np.isnan(observations)).sum()
    print(f"  {len(dates)} days, {n_obs} observations")

    # Tune parameters
    q_tdee, mean_revert = tune_params(calories, observations, expected_rmr)

    # Run filter with selected parameters
    print(f"\nRunning Kalman filter (Q_tdee={q_tdee}, mean_revert={mean_revert})...")
    x_filt, P_filt, x_pred, P_pred, innovations = kalman_forward(
        calories, observations, expected_rmr, q_tdee, mean_revert=mean_revert)

    print("Running RTS smoother...")
    x_smooth, P_smooth = rts_smoother(x_filt, P_filt, x_pred, P_pred, calories,
                                       mean_revert=mean_revert)

    validate(dates, x_smooth, P_smooth, innovations, observations, lean, expected_rmr)
    save_output(dates, x_filt, P_filt, x_smooth, P_smooth, innovations)
    plot_kalman(dates, x_smooth, P_smooth, innovations, observations, lean,
                expected_rmr, calories)


if __name__ == "__main__":
    main()
