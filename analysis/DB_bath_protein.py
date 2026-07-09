"""Test the Lamming-Lab protein/isoleucine → body temp depression hypothesis.

Mechanism: an inability to fully process dietary isoleucine causes downstream
body-temperature depression. Higher protein (proxy for isoleucine) →
cooler body → cold feeling → warming bath behavior. Predicts a POSITIVE
correlation between protein and (a) bath occurrence and (b) lower baseline
body temperature, possibly with a 1-2 day lag for the metabolic processing
to land.

Two endpoints:
  1. Bath occurrence (behavioral proxy for feeling cold).
     baths/baths.csv since Sept 2020 (~440 events in this window).
  2. Circadian-corrected daily basal body temperature.
     analysis/AS_temperature_baseline.csv since Dec 2023 (~360 days).

Intake from intake/intake_daily.csv (food-scale logged, ~10% known undercount
uniform across phases — see CLAUDE.md).

Reported metrics:
  - Bath endpoint: point-biserial r and logistic OR per +10 g protein for
    same-day, lag-1, lag-2, and 3-day trailing mean. Seasonal subsets.
  - Body-temp endpoint: Pearson r between protein (same-day, lag-1, lag-2,
    3-day mean) and baseline_mean.

Artifacts:
  - analysis/DB_bath_protein_summary.csv (bath endpoint)
  - analysis/DB_temp_protein_summary.csv (body-temp endpoint)
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
BATHS = ROOT / "baths" / "baths.csv"
INTAKE = ROOT / "intake" / "intake_daily.csv"
TEMP = ROOT / "analysis" / "AS_temperature_baseline.csv"
OUT_BATH = ROOT / "analysis" / "DB_bath_protein_summary.csv"
OUT_TEMP = ROOT / "analysis" / "DB_temp_protein_summary.csv"


def load_daily():
    intake = pd.read_csv(INTAKE, parse_dates=["date"])
    bath_dates = pd.read_csv(BATHS, parse_dates=["date"])["date"]

    # Build a daily frame over the bath-observation window only (Sept 2020 → today)
    start = bath_dates.min()
    end = bath_dates.max()
    days = pd.DataFrame({"date": pd.date_range(start, end, freq="D")})
    d = days.merge(intake[["date", "protein_g", "calories"]], on="date", how="left")
    d["bath"] = d["date"].isin(set(bath_dates)).astype(int)

    # Lagged & trailing protein
    d["protein_lag1"] = d["protein_g"].shift(1)
    d["protein_lag2"] = d["protein_g"].shift(2)
    d["protein_mean3"] = d["protein_g"].rolling(3).mean().shift(1)  # prior 3 days
    d["month"] = d["date"].dt.month
    return d.dropna(subset=["protein_g", "protein_lag1", "protein_lag2", "protein_mean3"])


def pb(d, x_col, label):
    """Point-biserial correlation + group means + bootstrap CI."""
    x = d[x_col].values
    y = d["bath"].values
    r, p = stats.pointbiserialr(y, x)
    bath_mean = x[y == 1].mean()
    nobath_mean = x[y == 0].mean()
    # Bootstrap CI for the difference in means
    diffs = []
    rng = np.random.default_rng(42)
    for _ in range(2000):
        idx = rng.integers(0, len(d), len(d))
        bx, by = x[idx], y[idx]
        if by.sum() == 0 or by.sum() == len(by): continue
        diffs.append(bx[by == 1].mean() - bx[by == 0].mean())
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return dict(label=label, n=len(d), n_bath=int(y.sum()),
                r=r, p=p, bath_mean=bath_mean, nobath_mean=nobath_mean,
                diff=bath_mean - nobath_mean, diff_lo=lo, diff_hi=hi)


def main():
    d = load_daily()

    # Subsets
    subsets = {
        "all": d,
        "winter (Dec/Jan/Feb)": d[d["month"].isin([12, 1, 2])],
        "spring (Mar/Apr/May)": d[d["month"].isin([3, 4, 5])],
        "summer (Jun/Jul/Aug)": d[d["month"].isin([6, 7, 8])],
        "fall (Sep/Oct/Nov)": d[d["month"].isin([9, 10, 11])],
    }
    cols = [("protein_g",   "protein same day"),
            ("protein_lag1","protein lag-1 day"),
            ("protein_lag2","protein lag-2 days"),
            ("protein_mean3","protein mean prior 3 days")]

    rows = []
    print(f"\n=== Bath ~ protein correlation tests ===\n")
    for sub_name, sub in subsets.items():
        print(f"--- {sub_name}  (n={len(sub)}, baths={int(sub['bath'].sum())}) ---")
        for col, lbl in cols:
            res = pb(sub, col, lbl)
            res["subset"] = sub_name
            rows.append(res)
            sig = "*" if res["p"] < 0.05 else " "
            print(f"  {lbl:<26} r={res['r']:+.3f}  p={res['p']:.3f}{sig}  "
                  f"bath_mean={res['bath_mean']:5.1f}g  "
                  f"nobath_mean={res['nobath_mean']:5.1f}g  "
                  f"Δ={res['diff']:+.1f}g  CI=[{res['diff_lo']:+.1f},{res['diff_hi']:+.1f}]")
        print()

    # Logistic regression on "all": log-odds per +10g protein
    from numpy.linalg import LinAlgError
    print("=== Logistic regression (all days) ===")
    for col, lbl in cols:
        x = d[col].values
        y = d["bath"].values
        # Simple Newton-Raphson logistic
        X = np.column_stack([np.ones(len(x)), x / 10.0])  # scale to +10 g
        beta = np.zeros(2)
        for _ in range(30):
            p = 1 / (1 + np.exp(-X @ beta))
            W = np.diag(p * (1 - p))
            try:
                beta_new = beta + np.linalg.solve(X.T @ W @ X, X.T @ (y - p))
            except LinAlgError:
                break
            if np.allclose(beta_new, beta, atol=1e-7): break
            beta = beta_new
        # std errors
        p = 1 / (1 + np.exp(-X @ beta))
        W = np.diag(p * (1 - p))
        try:
            cov = np.linalg.inv(X.T @ W @ X)
            se = np.sqrt(np.diag(cov))
            z = beta[1] / se[1]
            p_val = 2 * (1 - stats.norm.cdf(abs(z)))
            print(f"  {lbl:<26} log-OR per +10g: {beta[1]:+.4f}  "
                  f"OR={np.exp(beta[1]):.3f}  z={z:+.2f}  p={p_val:.3f}")
        except LinAlgError:
            print(f"  {lbl:<26} (singular)")

    pd.DataFrame(rows).to_csv(OUT_BATH, index=False)
    print(f"\nArtifact: {OUT_BATH.relative_to(ROOT)}")

    # ---- Body-temperature endpoint ----
    print("\n\n=== Body-temperature endpoint (Dec 2023 onward) ===\n")
    temp = pd.read_csv(TEMP, parse_dates=["date"])
    intake = pd.read_csv(INTAKE, parse_dates=["date"])[["date","protein_g","calories"]]
    tdf = temp[["date","baseline_mean","n_readings"]].merge(intake, on="date", how="left")
    tdf = tdf.sort_values("date").reset_index(drop=True)
    tdf["protein_lag1"]  = tdf["protein_g"].shift(1)
    tdf["protein_lag2"]  = tdf["protein_g"].shift(2)
    tdf["protein_mean3"] = tdf["protein_g"].rolling(3).mean().shift(1)
    tdf["month"] = tdf["date"].dt.month
    tdf = tdf.dropna(subset=["baseline_mean","protein_g","protein_lag1",
                             "protein_lag2","protein_mean3"])

    temp_rows = []
    for sub_name, sub in [
            ("all", tdf),
            ("winter (Dec/Jan/Feb)", tdf[tdf.month.isin([12,1,2])]),
            ("summer (Jun/Jul/Aug)", tdf[tdf.month.isin([6,7,8])]),
    ]:
        if len(sub) < 30:
            print(f"--- {sub_name}: only {len(sub)} days, skipping ---\n")
            continue
        print(f"--- {sub_name}  (n={len(sub)}  temp mean={sub.baseline_mean.mean():.2f}°F) ---")
        for col, lbl in cols:
            r, p = stats.pearsonr(sub[col], sub["baseline_mean"])
            slope, intercept, _, _, _ = stats.linregress(sub[col], sub["baseline_mean"])
            sig = "*" if p < 0.05 else " "
            print(f"  {lbl:<26} r={r:+.3f}  p={p:.3f}{sig}  "
                  f"slope={slope*10:+.3f} °F per +10 g protein")
            temp_rows.append(dict(subset=sub_name, predictor=lbl, n=len(sub),
                                  r=r, p=p, slope_per_10g=slope*10))
        print()

    pd.DataFrame(temp_rows).to_csv(OUT_TEMP, index=False)
    print(f"Artifact: {OUT_TEMP.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
