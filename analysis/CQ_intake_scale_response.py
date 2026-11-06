"""CQ. Non-circular estimate of the intake scale factor from weight response.

Why this exists
---------------
Every existing estimate of the intake undercount is anchored on an assumption
rather than a measurement:

  - A_energy_balance_quality.py computes undercount as
    `expected_rmr × 1.2 − kalman_tdee`, so the answer is a linear function of
    the stipulated activity factor 1.2.
  - CI/CJ solve for the multiplier by driving the Kalman's cumulative
    eb_residual_lb to zero. But TDEE is a free random walk, so scaling intake by
    m raises TDEE by (m−1)×intake almost day for day. Measured leakage:
    dEB/dm is 24.9 lb per unit m against the 3330.7 lb/unit it would be if TDEE
    were fixed — the filter absorbs 99.25% of any level change. CJ's protein
    sweep contains two configs differing by 8.4% of total lifetime intake that
    both close EB to within 0.002 lb.

So the level is not identified by anything currently in the pipeline.

What IS identifiable
--------------------
A constant absolute undercount is absorbed 1:1 into TDEE and is genuinely
unidentifiable from weight + intake. A PROPORTIONAL undercount is identifiable,
because intake varies by hundreds of calories/day across windows and TDEE cannot
follow that variation (rate cap 100 cal/day, and real adaptation is far slower).

Over a window between two weigh-ins (t1, t2), n = t2 − t1 days:

    rho × (W[t2] − W[t1]) = k × sum(intake[t1 .. t2−1]) − TDEE × n
    =>  y = rho×dW/n = k × mean_intake − TDEE

Block-demean y and mean_intake to difference out slow TDEE drift. Slope = k.

Inputs are Kalman-free: P1 smoothed weight (glycogen+sodium-corrected scale
readings), raw logged intake, and P2 scan-interpolated FM for the Forbes energy
density. Nothing here has seen INTAKE_MULTIPLIER.

Windows are one per weigh-in, paired with the weigh-in nearest +TARGET days
(median consecutive weigh-in gap is 1 day, so restricting to consecutive pairs
would discard ~90% of the data). Windows overlap, so inference is by block
bootstrap over year-length blocks.

Bias is measured, not assumed: simulate_recovery() generates trajectories from
the real intake series under a known k with realistic AR(1) water noise and
measurement noise, subsamples to the real weigh-in dates, and runs the same
estimator. The correction applied to the headline number comes from that.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
import sys
sys.path.insert(0, str(ROOT / "analysis"))
from intake_model import INTAKE_MULTIPLIER, describe as describe_intake_model  # noqa: E402
from results import save as save_result, load as load_result  # noqa: E402

E_FAT, E_LEAN = 3500.0, 1100.0
FORBES_C_KG, LBS_TO_KG = 10.4, 0.453592

TARGET, TOL = 21, 4          # headline window: 3 weeks — the middle of DM's 14-28 d criterion; shorter windows are water-inflated (see table)
BLOCK_DAYS = 180             # TDEE-drift removal
BOOT_BLOCK = 365             # bootstrap resampling unit
BUF_HL, BUF_SD, MEAS_SD = 5.0, 0.55, 0.6   # noise matched to filter innov std ~1.15 lb


def forbes_rho(fm_lbs):
    """Effective energy density of weight change (cal/lb) under Forbes partition."""
    fm_kg = np.maximum(np.asarray(fm_lbs) * LBS_TO_KG, 1.0)
    ff = 1.0 - (FORBES_C_KG / (FORBES_C_KG + fm_kg)) ** 2
    return ff * E_FAT + (1.0 - ff) * E_LEAN


def load():
    p1 = (pd.read_csv(ROOT / "analysis/P1_smoothed_weight.csv", parse_dates=["date"])
          .dropna(subset=["smoothed_weight_lbs"]).sort_values("date").reset_index(drop=True))
    intake = pd.read_csv(ROOT / "intake/intake_daily.csv", parse_dates=["date"])
    p2 = pd.read_csv(ROOT / "analysis/P2_daily_composition.csv", parse_dates=["date"])
    fm = p2.set_index("date")["fm_lbs"].reindex(intake.date).ffill().bfill()
    return p1, intake, fm


def build_windows(p1, intake, fm, target=TARGET, tol=TOL):
    cal = intake.set_index("date")["calories"]
    rho_daily = pd.Series(forbes_rho(fm.values), index=intake.date)
    d = p1.date.values
    w = p1.smoothed_weight_lbs.values
    rows = []
    for a in range(len(d)):
        lo = d[a] + np.timedelta64(target - tol, "D")
        hi = d[a] + np.timedelta64(target + tol, "D")
        cand = np.where((d >= lo) & (d <= hi))[0]
        if len(cand) == 0:
            continue
        b = cand[np.argmin(np.abs(d[cand] - (d[a] + np.timedelta64(target, "D"))))]
        n = int((d[b] - d[a]) / np.timedelta64(1, "D"))
        if n <= 0:
            continue
        # off-by-one: morning weight reflects the PREVIOUS day's intake, so the
        # window covers intake on days t1 .. t2-1
        days = pd.date_range(pd.Timestamp(d[a]), pd.Timestamp(d[b]) - pd.Timedelta(days=1))
        c = cal.reindex(days)
        if c.isna().any():
            continue
        rows.append(dict(t=pd.Timestamp(d[a]), n=n, x=c.mean(),
                         y=rho_daily.reindex(days).mean() * (w[b] - w[a]) / n))
    return pd.DataFrame(rows)


def demean(df, block_days=BLOCK_DAYS):
    df = df.copy()
    df["block"] = ((df.t - df.t.min()).dt.days // block_days).values
    df = df[df.groupby("block")["y"].transform("size") >= 5].copy()
    for c in ("x", "y"):
        df[c + "d"] = df[c] - df.groupby("block")[c].transform("mean")
    return df


def slope(x, y):
    return float(np.sum(x * y) / np.sum(x * x))


def bootstrap(df, boot_block=BOOT_BLOCK, n_boot=3000, seed=0):
    rng = np.random.default_rng(seed)
    bb = ((df.t - df.t.min()).dt.days // boot_block).values
    uniq = np.unique(bb)
    by = {b: np.where(bb == b)[0] for b in uniq}
    xd, yd = df.xd.values, df.yd.values
    out = [slope(xd[i], yd[i]) for i in
           (np.concatenate([by[b] for b in rng.choice(uniq, len(uniq), True)])
            for _ in range(n_boot))]
    return tuple(np.percentile(out, [2.5, 97.5]))


def simulate_recovery(p1, intake, fm, k_true, n_seeds=200, tdee_mode="kalman"):
    """Measure estimator bias: simulate from a known k and re-estimate."""
    k4 = pd.read_csv(ROOT / "analysis/P4_kalman_daily.csv", parse_dates=["date"])
    cal = intake.calories.values.astype(float)
    rho_path = forbes_rho(fm.values)
    tdee_ref = k4.set_index("date")["tdee"].reindex(intake.date).ffill().bfill().values
    if tdee_mode == "constant":
        tdee_ref = np.full(len(cal), float(np.mean(tdee_ref)))
    tdee = tdee_ref * (k_true / 1.09)
    dW = (k_true * cal - tdee) / rho_path
    base = 200.0 + np.concatenate([[0.0], np.cumsum(dW)[:-1]])   # w[i] reflects intake thru i-1
    decay = np.exp(-np.log(2) / BUF_HL)

    dates = intake.date.values
    weigh = intake.date.isin(p1.date).values
    idx = np.where(weigh)[0]
    dd = pd.to_datetime(dates[idx])
    A, B, X, R, T = [], [], [], [], []
    for a in range(len(idx)):
        lo, hi = dd[a] + pd.Timedelta(days=TARGET - TOL), dd[a] + pd.Timedelta(days=TARGET + TOL)
        c = np.where((dd >= lo) & (dd <= hi))[0]
        if len(c) == 0:
            continue
        b = c[np.argmin(np.abs((dd[c] - (dd[a] + pd.Timedelta(days=TARGET))).values))]
        i0, i1 = idx[a], idx[b]
        if i1 <= i0:
            continue
        A.append(i0); B.append(i1); X.append(cal[i0:i1].mean())
        R.append(rho_path[i0:i1].mean() / (i1 - i0))
        T.append(dd[a])          # keep dates aligned with the windows actually kept
    A, B, X, R = map(np.array, (A, B, X, R))
    t = pd.Series(pd.to_datetime(T))
    blk = ((t - t.min()).dt.days // BLOCK_DAYS).values
    keep = pd.Series(blk).groupby(blk).transform("size").values >= 5
    A, B, X, R, blk = A[keep], B[keep], X[keep], R[keep], blk[keep]
    xd = X - pd.Series(X).groupby(blk).transform("mean").values

    ests = []
    for s in range(n_seeds):
        rng = np.random.default_rng(10_000 + s)
        buf = np.zeros(len(cal))
        for i in range(1, len(cal)):
            buf[i] = buf[i - 1] * decay + rng.normal(0, BUF_SD)
        w = base + buf + rng.normal(0, MEAS_SD, len(cal))
        y = R * (w[B] - w[A])
        yd = y - pd.Series(y).groupby(blk).transform("mean").values
        ests.append(slope(xd, yd))
    return float(np.mean(ests)), float(np.std(ests, ddof=1) / np.sqrt(n_seeds))


def main():
    p1, intake, fm = load()

    print("=== CQ. Intake scale factor k from weight response to intake variation ===")
    print("k = 1.00 means logged intake is unbiased in SCALE.")
    print("A constant absolute offset is NOT identifiable here (TDEE absorbs it 1:1).")
    print(f"Canonical intake model: {describe_intake_model()}\n")

    print("Window-length sensitivity (180d TDEE blocks):")
    raw = {}
    for tgt, tol in [(10, 3), (14, 3), (21, 4), (28, 5), (42, 7), (60, 10)]:
        d = demean(build_windows(p1, intake, fm, tgt, tol))
        k = slope(d.xd.values, d.yd.values)
        lo, hi = bootstrap(d)
        r = float(np.corrcoef(d.xd, d.yd)[0, 1])
        raw[tgt] = k
        print(f"  +{tgt:>2}d  n={len(d):>5}  r={r:+.3f}  k_raw={k:.3f}  95%CI[{lo:.3f}, {hi:.3f}]")

    print("\nEstimator bias (simulate known k, re-estimate; 200 seeds/cell):")
    bias = {}
    for mode in ("kalman", "constant"):
        for kt in (1.00, 1.20, 1.40):
            m, se = simulate_recovery(p1, intake, fm, kt, 200, mode)
            bias[(mode, kt)] = m - kt
            print(f"  TDEE={mode:<9} k_true={kt:.2f}  k_hat={m:.3f}±{se:.3f}  "
                  f"bias={m - kt:+.3f}")

    corr = -np.mean([bias[("kalman", k)] for k in (1.00, 1.20, 1.40)])
    d = demean(build_windows(p1, intake, fm))
    k_raw = slope(d.xd.values, d.yd.values)
    lo, hi = bootstrap(d)
    print(f"\n--- Headline ({TARGET}d windows, bias-corrected by {corr:+.3f}; DM's criterion spans 14-28 d) ---")
    print(f"  k = {k_raw + corr:.2f}   95% CI [{lo + corr:.2f}, {hi + corr:.2f}]")

    # implied activity factor: TDEE/RMR on calorimetry days scales ~linearly with k.
    # The reference ratio comes from P4's own results file at the canonical multiplier.
    p4 = load_result("P4_kalman_filter", required=False)
    ratio_canon = (p4["diagnostics"]["tdee_rmr_mean"] if p4 else float("nan"))
    print(f"  implied TDEE/RMR on calorimetry days: "
          f"{ratio_canon * (k_raw + corr) / INTAKE_MULTIPLIER:.2f} "
          f"[{ratio_canon * (lo + corr) / INTAKE_MULTIPLIER:.2f}, {ratio_canon * (hi + corr) / INTAKE_MULTIPLIER:.2f}]")
    print(f"  (canonical multiplier {INTAKE_MULTIPLIER:.3f} gives {ratio_canon:.2f}; free-living adult PAL "
          f"is conventionally >=1.40)")
    save_result("CQ_intake_scale_response", k=k_raw + corr, k_ci=[lo + corr, hi + corr], k_raw=k_raw,
                bias_correction=corr, window_days=TARGET, by_window={str(t): v for t, v in raw.items()},
                implied_tdee_rmr=ratio_canon * (k_raw + corr) / INTAKE_MULTIPLIER)
    print("Wrote analysis/results/CQ_intake_scale_response.json")


if __name__ == "__main__":
    main()
