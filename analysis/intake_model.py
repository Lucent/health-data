"""Canonical logged-intake → consumed-energy correction. ONE source of truth.

Every script that needs consumed calories imports from here:

    from intake_model import true_intake, INTAKE_MULTIPLIER, CALIBRATION
    cal_true = true_intake(intake_df)        # array, same length as intake_df

The parameters are NOT typed into this file. They are read from
analysis/results/DM_intake_calibration.json, written by
analysis/DM_intake_calibration.py, which compares candidate correction forms
(uniform multiplier, affine, intake-level-dependent, macro-specific,
calendar-era, FM-direction) by the Kalman innovation likelihood and writes the
winner. P4 consumes `true_intake` and writes the corrected series into
P4_kalman_daily.csv as `intake_true`, so scripts that merge the Kalman output
can use that column directly instead of re-applying anything.

`INTAKE_MULTIPLIER` is the record-wide effective ratio Σ true / Σ logged. It is
exact for the multiplicative form and an approximation for the others; prefer
`true_intake` or P4's `intake_true` column wherever a daily series exists.

Bootstrap: if the calibration file is absent (fresh clone, DM not yet run), the
module falls back to a uniform multiplier of 1.0 and prints a warning. That is
deliberate — the fallback must be visibly *uncalibrated*, not a guess.
"""
from pathlib import Path
import sys
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analysis"))
from results import load as _load_result  # noqa: E402

# Calendar eras used by the "era" form (from DL's year-error pattern; the
# direction form below is the data-driven alternative that needs no dates).
ERA_BOUNDS = [
    ("loss",    "2011-01-01", "2016-12-31"),
    ("regain",  "2017-01-01", "2021-12-31"),
    ("plateau", "2022-01-01", "2024-09-16"),
    ("drug",    "2024-09-17", "2099-12-31"),
]
DIRECTION_SLOPE_LB_PER_DAY = 0.01   # |90-d scan-interpolated FM slope| below this = stable
REF_CAL = 2000.0                    # centring for the level-dependent form

MODEL_PARAM_NAMES = {
    "multiplicative": ["m"],
    "affine":         ["a", "m"],
    "level":          ["m", "m1"],
    "macro":          ["m", "a_fiber", "a_protein"],
    "era":            [f"m_{e[0]}" for e in ERA_BOUNDS],
    "direction":      ["m_losing", "m_stable", "m_gaining"],
}


def _fm_direction(dates):
    """losing / stable / gaining from the 90-day slope of P2's scan-interpolated
    fat mass (no Kalman state involved)."""
    p2 = pd.read_csv(ROOT / "analysis" / "P2_daily_composition.csv", parse_dates=["date"])
    fm = p2.set_index("date")["fm_lbs"].reindex(pd.DatetimeIndex(dates)).ffill().bfill()
    slope = (fm.rolling(90, center=True, min_periods=30).mean().diff(30) / 30).ffill().bfill().values
    out = np.full(len(dates), "stable", dtype=object)
    out[slope < -DIRECTION_SLOPE_LB_PER_DAY] = "losing"
    out[slope > DIRECTION_SLOPE_LB_PER_DAY] = "gaining"
    return out


def era_of(dates):
    d = pd.DatetimeIndex(dates)
    out = np.full(len(d), ERA_BOUNDS[-1][0], dtype=object)
    for name, lo, hi in ERA_BOUNDS:
        out[(d >= lo) & (d <= hi)] = name
    return out


def apply_model(model, params, df):
    """Consumed energy for each row of `df` (needs `calories`; `date`,
    `fiber_g`, `protein_g` for the forms that use them)."""
    cal = df["calories"].values.astype(float)
    if model == "multiplicative":
        return params["m"] * cal
    if model == "affine":
        return params["a"] + params["m"] * cal
    if model == "level":
        return cal * (params["m"] + params["m1"] * (cal - REF_CAL) / 1000.0)
    if model == "macro":
        fib = df["fiber_g"].fillna(0).values.astype(float)
        pro = df["protein_g"].fillna(0).values.astype(float)
        return params["m"] * cal + params["a_fiber"] * fib + params["a_protein"] * pro
    if model == "era":
        e = era_of(df["date"].values)
        m = np.array([params[f"m_{x}"] for x in e])
        return m * cal
    if model == "direction":
        dirn = _fm_direction(df["date"].values)
        m = np.array([params[f"m_{x}"] for x in dirn])
        return m * cal
    raise ValueError(f"unknown intake model {model!r}")


_cal = _load_result("DM_intake_calibration", required=False)
if _cal is None:
    print("[intake_model] WARNING: analysis/results/DM_intake_calibration.json not found; "
          "using an UNCALIBRATED multiplier of 1.0. Run analysis/DM_intake_calibration.py.",
          file=sys.stderr)
    CALIBRATION = {"model": "multiplicative", "params": {"m": 1.0},
                   "effective_multiplier": 1.0, "calibrated": False}
else:
    CALIBRATION = dict(_cal)
    CALIBRATION.setdefault("calibrated", True)

MODEL = CALIBRATION["model"]
PARAMS = CALIBRATION["params"]
INTAKE_MULTIPLIER = float(CALIBRATION["effective_multiplier"])
# Filter noise hyperparameters selected by the same likelihood (DM). P4 reads
# these for its canonical run; the fallback is P4's hand-set tuning.
FILTER = dict(sigma_tdee=50.0, q_buf=0.32, r_weight=0.97, max_dtdee=100, beta_buf=0.0, buf_hl=5.0)
FILTER.update(CALIBRATION.get("filter") or {})


def true_intake(df):
    return apply_model(MODEL, PARAMS, df)


def load_intake():
    """intake_daily.csv with a `true_cal` column attached."""
    df = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    df["true_cal"] = true_intake(df)
    return df


def describe():
    p = ", ".join(f"{k}={v:.4g}" for k, v in PARAMS.items())
    return f"{MODEL} ({p}); effective multiplier {INTAKE_MULTIPLIER:.3f}"


def describe_filter():
    return (f"σ_TDEE {FILTER['sigma_tdee']:.1f} cal/d, Q_BUF {FILTER['q_buf']:.3f} lb², "
            f"R_WEIGHT {FILTER['r_weight']:.2f} lb², BETA_BUF {FILTER['beta_buf']:.3f} lb/1000 cal, "
            f"BUF_HL {FILTER['buf_hl']:.1f} d, max_dTDEE {FILTER['max_dtdee']:.0f} cal/d")
