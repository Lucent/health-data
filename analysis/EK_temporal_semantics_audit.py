"""Audit temporal semantics in legacy analysis scripts.

This is a narrow, deterministic audit. It does not decide whether every
association in the repository is causal. It inventories code paths where the
implemented time window differs from the prose, where a sparse-table row shift
is called a calendar-day shift, or where a retrospective feature is presented
as forward-looking.

Outputs
-------
analysis/EK_temporal_semantics_inventory.csv
analysis/EK_temporal_semantics_checks.csv
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
ANALYSIS_DIR = ROOT / "analysis"
INVENTORY_PATH = ANALYSIS_DIR / "EK_temporal_semantics_inventory.csv"
CHECKS_PATH = ANALYSIS_DIR / "EK_temporal_semantics_checks.csv"


# Each needle is matched against the current source so the exported line number
# cannot silently drift when a file is edited.
AUDIT_ITEMS = [
    {
        "script": "U_steps_compensation.py",
        "needle": 'daily["calories"].shift(-1).rolling(horizon',
        "pattern": "shift(-1).rolling(h)",
        "implemented_semantics": "At t, h=7 uses t-5..t+1 and h=14 uses t-12..t+1.",
        "claimed_semantics": "Mean of the following 7 or 14 calendar days.",
        "claim_impact": "Invalidates the promoted U future-intake window and the analogous future TDEE/RMR window.",
        "severity": "high",
        "disposition": "Already corrected by EI direct re-audit; do not use U future outcomes.",
    },
    {
        "script": "V_exercise_walk_analysis.py",
        "needle": 'daily["calories"].shift(-1).rolling(',
        "pattern": "shift(-1).rolling(h)",
        "implemented_semantics": "The nominal future-14 window is mostly trailing and is also indexed on available step rows.",
        "claimed_semantics": "Mean of the following 14 calendar days.",
        "claim_impact": "Invalidates the promoted paired-walk future-14 intake contrast.",
        "severity": "high",
        "disposition": "Already corrected by EI direct re-audit; do not use V future outcomes.",
    },
    {
        "script": "AE_sleep_null.py",
        "needle": 'daily["cal_next"] = daily["calories"].shift(-1)',
        "pattern": "shift(-1) after anchoring on sleep rows",
        "implemented_semantics": "Next recorded sleep row, which is not always the next calendar date.",
        "claimed_semantics": "Next calendar-day intake.",
        "claim_impact": "Invalidates AE's exact next-day numbers, especially in sparse early sleep years.",
        "severity": "high",
        "disposition": "Already corrected by EI exact-calendar re-audit.",
    },
    {
        "script": "BS_sleep_timing_split.py",
        "needle": 'df["cal_next"] = df["calories"].shift(-1)',
        "pattern": "shift(-1) after anchoring on sleep rows",
        "implemented_semantics": "Next recorded sleep row, not necessarily the next calendar date.",
        "claimed_semantics": "Next calendar-day intake.",
        "claim_impact": "Makes the legacy next-day bedtime/sleep-debt endpoints unreliable.",
        "severity": "high",
        "disposition": "Already corrected by EI exact-calendar re-audit.",
    },
    {
        "script": "AF_intake_variance.py",
        "needle": 'df[f"fm_change_{w}d"] = df["fat_mass_lbs"].shift(-w)',
        "pattern": "trailing exposure paired to non-overlapping future outcome",
        "implemented_semantics": "Intake mean/variance/TDEE use t-w+1..t; FM change uses t..t+w.",
        "claimed_semantics": "Fat gain caused during a window by that same window's mean intake and variance.",
        "claim_impact": "Does not establish that variable overeating is less fattening at the same contemporaneous surplus.",
        "severity": "high",
        "disposition": "Quarantine the causal/protective interpretation; both intervals also use legacy P4 latent states.",
    },
    {
        "script": "CL_revealed_sp_from_binges.py",
        "needle": 'pre[col] = pre["binge"].rolling(window',
        "pattern": "centered predictor overlaps future outcome",
        "implemented_semantics": "The 30d/90d centered binge-rate predictor contains future binges.",
        "claimed_semantics": "A revealed current gap that predicts later binges.",
        "claim_impact": "The predictor contains all t+1..t+7 outcome days; the NEXT-7d diagnostic is target leakage.",
        "severity": "high",
        "disposition": "Do not use CL's next-7-day prediction diagnostic.",
    },
    {
        "script": "AB_running_vs_walking.py",
        "needle": 'daily["next_3d_cal"] = daily["calories"].rolling(3).mean().shift(-1)',
        "pattern": "rolling(3).shift(-1)",
        "implemented_semantics": "At t, averages t-1, t, and t+1.",
        "claimed_semantics": "Mean intake on t+1, t+2, and t+3.",
        "claim_impact": "The documented next-3d run-versus-walk contrast has the wrong horizon.",
        "severity": "medium",
        "disposition": "Use the exact-calendar recomputation in EK checks; same/next-day results are separate.",
    },
    {
        "script": "N_dietary_predictors.py",
        "needle": 'v["next_cal"] = v["calories"].shift(-1)',
        "occurrence": 2,
        "pattern": "shift(-1) after dropna on steps and modeled FM",
        "implemented_semantics": "Next row with both steps and modeled FM.",
        "claimed_semantics": "Next calendar-day intake.",
        "claim_impact": "The exact gravitostat coefficient is not an exact next-day estimate.",
        "severity": "medium",
        "disposition": "Use the exact-calendar comparison in EK checks; the qualitative null is re-evaluated there.",
    },
    {
        "script": "AK_sunlight_exposure.py",
        "needle": 'df["next_cal"] = df["calories"].shift(-1)',
        "pattern": "shift(-1) on sunlight/sleep-anchored rows",
        "implemented_semantics": "Next available sunlight row, not guaranteed to be the next calendar day.",
        "claimed_semantics": "Next calendar-day intake.",
        "claim_impact": "The legacy next-day sunlight endpoints are mislabeled; same-day/RMR analyses are separate.",
        "severity": "medium",
        "disposition": "Do not quote AK next-day coefficients without an exact date join.",
    },
    {
        "script": "CA_stress_daily.py",
        "needle": 'daily["cal_next"] = daily["calories"].shift(-1)',
        "pattern": "shift(-1) on 208 sparse stress days",
        "implemented_semantics": "Next recorded stress day, potentially many calendar days later.",
        "claimed_semantics": "Next calendar-day intake.",
        "claim_impact": "Invalidates CA's exported next-day stress endpoint, but not its promoted same-day sparse association.",
        "severity": "medium",
        "disposition": "Retain same-day result as exploratory; rebuild next-day result by exact date if stress data are restored.",
    },
    {
        "script": "D_food_noise_variance.py",
        "needle": 'rolling(180, center=True',
        "pattern": "centered 180-day latent baseline",
        "implemented_semantics": "The daily 'set point' includes up to about 90 days of future modeled FM.",
        "claimed_semantics": "A contemporaneous defended-weight distance influencing intake.",
        "claim_impact": "The distance-to-intake slope is retrospective detrending, not a forward appetite-pressure test.",
        "severity": "medium",
        "disposition": "Treat the D slope as exploratory and retrospective only.",
    },
    {
        "script": "DA5_intake_pressure_foods.py",
        "needle": 'rolling(30, min_periods=15, center=True)',
        "pattern": "centered same-day outcome baseline",
        "implemented_semantics": "A food day's residual is centered against roughly 15 prior and 14 future intake days.",
        "claimed_semantics": "Food-specific downward or upward intake pressure.",
        "claim_impact": "Useful retrospective de-trending, but not evidence that a food prospectively changes intake.",
        "severity": "medium",
        "disposition": "Label rankings associative; use EB's ordered holdout features for prospective claims.",
    },
    {
        "script": "DA6_food_macros.py",
        "needle": 'rolling(30, min_periods=15, center=True)',
        "pattern": "centered same-day outcome baseline",
        "implemented_semantics": "Trigger/protective labels incorporate future intake through the centered residual.",
        "claimed_semantics": "Macro features distinguish foods that trigger or protect against intake.",
        "claim_impact": "Cannot support prospective trigger/protective food classification.",
        "severity": "medium",
        "disposition": "Keep as exploratory retrospective clustering only.",
    },
    {
        "script": "CR_jump_events.py",
        "needle": 'pre["fm"].shift(-15) - pre["fm"].shift(15)',
        "pattern": "centered event detector",
        "implemented_semantics": "Each event label uses modeled FM 15 records before and after the labeled day.",
        "claimed_semantics": "The module docstring calls this a 30-day forward slope.",
        "claim_impact": "Events are retrospective segments and cannot be interpreted as real-time trigger detection.",
        "severity": "medium",
        "disposition": "Acceptable for descriptive segmentation after relabeling; not predictive.",
    },
    {
        "script": "CS_cycle_peaks_troughs.py",
        "needle": 'rolling(14, min_periods=7, center=True)',
        "pattern": "centered peak/trough detector",
        "implemented_semantics": "Peak/trough dates require future modeled FM observations.",
        "claimed_semantics": "Retrospective cycle decomposition.",
        "claim_impact": "No defect if kept retrospective; unsuitable for online detection or forecast claims.",
        "severity": "low",
        "disposition": "Retrospective descriptive use only.",
    },
    {
        "script": "CT_cycle_pairing.py",
        "needle": 'rolling(14, min_periods=7, center=True)',
        "pattern": "centered peak/trough detector",
        "implemented_semantics": "Cycle boundaries depend on future modeled FM.",
        "claimed_semantics": "Retrospective adjacent-cycle pairing.",
        "claim_impact": "No defect if kept retrospective; subsequent-cycle associations remain post hoc.",
        "severity": "low",
        "disposition": "Retrospective descriptive use only.",
    },
    {
        "script": "A_energy_balance_quality.py",
        "needle": 'rolling(90, center=True',
        "pattern": "centered phase classification",
        "implemented_semantics": "Gaining/losing/stable labels use roughly 45 future days.",
        "claimed_semantics": "Phase-stratified descriptive undercount.",
        "claim_impact": "Not a forward predictor; phase labels are acceptable only as retrospective annotations.",
        "severity": "low",
        "disposition": "Retrospective description only; model-dependence is already separately documented.",
    },
    {
        "script": "Y_set_point_intake_tests.py",
        "needle": 'daily["calories"].rolling(7).mean().shift(-7)',
        "pattern": "rolling(7).shift(-7)",
        "implemented_semantics": "At t, exactly averages t+1 through t+7 on the complete intake calendar.",
        "claimed_semantics": "Following seven calendar days.",
        "claim_impact": "Verified control: this superficially suspicious expression is temporally correct.",
        "severity": "none",
        "disposition": "No temporal correction required; other model-dependence remains separate.",
    },
]


def locate_line(script: str, needle: str, occurrence: int = 1) -> tuple[int, str]:
    path = ANALYSIS_DIR / script
    lines = path.read_text(encoding="utf-8").splitlines()
    matches = [(i, line) for i, line in enumerate(lines, start=1) if needle in line]
    if len(matches) < occurrence:
        raise RuntimeError(
            f"Expected occurrence {occurrence} of {needle!r} in {script}; "
            f"found {len(matches)}"
        )
    return matches[occurrence - 1]


def build_inventory() -> pd.DataFrame:
    rows = []
    for item in AUDIT_ITEMS:
        row = item.copy()
        occurrence = int(row.pop("occurrence", 1))
        line, source = locate_line(row["script"], row.pop("needle"), occurrence)
        row["line"] = line
        row["source"] = source.strip()
        rows.append(row)
    columns = [
        "script",
        "line",
        "severity",
        "pattern",
        "implemented_semantics",
        "claimed_semantics",
        "claim_impact",
        "disposition",
        "source",
    ]
    return pd.DataFrame(rows)[columns].sort_values(
        ["severity", "script", "line"],
        key=lambda s: (
            s.map({"high": 0, "medium": 1, "low": 2, "none": 3})
            if s.name == "severity"
            else s
        ),
    )


def add_check(
    rows: list[dict[str, object]],
    check: str,
    metric: str,
    value: object,
    notes: str,
) -> None:
    rows.append(
        {
            "check": check,
            "metric": metric,
            "value": value,
            "notes": notes,
        }
    )


def partial_corr_controlling_today(
    frame: pd.DataFrame, x_col: str, y_col: str
) -> tuple[int, float]:
    d = frame[[x_col, y_col, "calories"]].dropna()
    x = d[x_col].to_numpy(float)
    y = d[y_col].to_numpy(float)
    controls = np.column_stack([d["calories"].to_numpy(float), np.ones(len(d))])
    x_resid = x - controls @ np.linalg.lstsq(controls, x, rcond=None)[0]
    y_resid = y - controls @ np.linalg.lstsq(controls, y, rcond=None)[0]
    return len(d), float(np.corrcoef(x_resid, y_resid)[0, 1])


def gap_checks(rows: list[dict[str, object]]) -> None:
    sources = {
        "sleep_rows": (
            ROOT / "steps-sleep" / "sleep.csv",
            lambda d: d[["date"]].drop_duplicates("date", keep="last"),
        ),
        "sunlight_rows": (
            ROOT / "steps-sleep" / "sunlight.csv",
            lambda d: d[["date"]].drop_duplicates("date", keep="last"),
        ),
        "step_rows": (
            ROOT / "steps-sleep" / "steps.csv",
            lambda d: d[["date"]].dropna().drop_duplicates("date", keep="last"),
        ),
    }
    for label, (path, transform) in sources.items():
        dates = transform(pd.read_csv(path, parse_dates=["date"])).sort_values("date")[
            "date"
        ]
        gaps = dates.diff().dt.days.dropna()
        add_check(rows, label, "rows", len(dates), "Unique source dates.")
        add_check(
            rows,
            label,
            "consecutive_transition_fraction",
            round(float((gaps == 1).mean()), 6),
            "Fraction of shift(-1) transitions that are exactly one calendar day.",
        )
        add_check(
            rows,
            label,
            "maximum_gap_days",
            int(gaps.max()),
            "Largest next-record date gap.",
        )


def n_gravitostat_check(rows: list[dict[str, object]]) -> None:
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    steps = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])
    kalman = pd.read_csv(
        ANALYSIS_DIR / "P4_kalman_daily.csv", parse_dates=["date"]
    )
    daily = (
        intake.merge(steps[["date", "steps"]], on="date", how="left")
        .merge(kalman[["date", "fat_mass_lbs"]], on="date", how="left")
        .sort_values("date")
        .reset_index(drop=True)
    )
    eligible = daily.dropna(subset=["steps", "fat_mass_lbs"]).copy()
    eligible["foot_lbs"] = eligible["steps"] * eligible["fat_mass_lbs"]
    eligible["legacy_next_cal"] = eligible["calories"].shift(-1)
    calendar_calories = intake.set_index("date")["calories"]
    eligible["exact_next_cal"] = (
        eligible["date"] + pd.Timedelta(days=1)
    ).map(calendar_calories)
    eligible["legacy_gap_days"] = eligible["date"].shift(-1).sub(
        eligible["date"]
    ).dt.days

    for target in ["legacy_next_cal", "exact_next_cal"]:
        n, r = partial_corr_controlling_today(eligible, "foot_lbs", target)
        add_check(
            rows,
            "N_gravitostat",
            f"{target}_partial_r",
            round(r, 6),
            f"n={n}; controls same-day calories, matching legacy N.",
        )
    add_check(
        rows,
        "N_gravitostat",
        "legacy_exact_next_day_fraction",
        round(float((eligible["legacy_gap_days"] == 1).mean()), 6),
        "Share of legacy row shifts that actually land on t+1.",
    )


def ab_running_check(rows: list[dict[str, object]]) -> None:
    exercises = pd.read_csv(
        ROOT / "steps-sleep" / "exercises.csv", parse_dates=["date"]
    )
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    steps = pd.read_csv(ROOT / "steps-sleep" / "steps.csv", parse_dates=["date"])

    daily = intake[["date", "calories"]].copy().sort_values("date")
    daily["legacy_next_3d"] = daily["calories"].rolling(3).mean().shift(-1)
    exact_future = pd.concat(
        [daily["calories"].shift(-lag) for lag in (1, 2, 3)], axis=1
    )
    daily["exact_next_3d"] = exact_future.mean(axis=1)

    runs = (
        exercises[exercises["type"] == "running"]
        .groupby("date")
        .size()
        .rename("run_sessions")
        .reset_index()
    )
    walks = (
        exercises[exercises["type"] == "walking"]
        .groupby("date")
        .size()
        .rename("walk_sessions")
        .reset_index()
    )
    daily = (
        daily.merge(runs, on="date", how="left")
        .merge(walks, on="date", how="left")
        .merge(steps[["date", "steps"]], on="date", how="left")
    )
    daily["ran"] = daily["run_sessions"].fillna(0) > 0
    daily["walked"] = daily["walk_sessions"].fillna(0) > 0
    era = daily[
        daily["date"].between(pd.Timestamp("2014-01-01"), pd.Timestamp("2016-12-31"))
    ].dropna(subset=["steps"])
    run_era = era[era["ran"]]
    walk_era = era[era["walked"] & ~era["ran"]]

    matches = []
    for _, run_row in run_era.iterrows():
        diffs = (walk_era["steps"] - run_row["steps"]).abs()
        if diffs.empty:
            continue
        best = walk_era.loc[diffs.idxmin()]
        if float(diffs.min()) < 2000:
            matches.append(
                {
                    "legacy_diff": run_row["legacy_next_3d"]
                    - best["legacy_next_3d"],
                    "exact_diff": run_row["exact_next_3d"] - best["exact_next_3d"],
                }
            )
    matched = pd.DataFrame(matches)
    add_check(rows, "AB_running", "matched_pairs", len(matched), "Legacy matching rule.")
    add_check(
        rows,
        "AB_running",
        "legacy_next_3d_run_minus_walk_kcal",
        round(float(matched["legacy_diff"].mean()), 3),
        "Actually t-1,t,t+1.",
    )
    add_check(
        rows,
        "AB_running",
        "exact_next_3d_run_minus_walk_kcal",
        round(float(matched["exact_diff"].mean()), 3),
        "Exact t+1,t+2,t+3.",
    )


def af_variance_check(rows: list[dict[str, object]]) -> None:
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    kalman = pd.read_csv(
        ANALYSIS_DIR / "P4_kalman_daily.csv", parse_dates=["date"]
    )
    daily = (
        intake[["date", "calories"]]
        .merge(kalman[["date", "fat_mass_lbs", "tdee"]], on="date", how="inner")
        .sort_values("date")
        .reset_index(drop=True)
    )

    for horizon in [14, 30]:
        mean_col = f"cal_mean_{horizon}"
        std_col = f"cal_std_{horizon}"
        surplus_col = f"surplus_{horizon}"
        daily[mean_col] = daily["calories"].rolling(
            horizon, min_periods=horizon
        ).mean()
        daily[std_col] = daily["calories"].rolling(
            horizon, min_periods=horizon
        ).std()
        tdee_mean = daily["tdee"].rolling(horizon, min_periods=horizon).mean()
        daily[surplus_col] = daily[mean_col] - tdee_mean
        daily["future_delta"] = (
            daily["fat_mass_lbs"].shift(-horizon) - daily["fat_mass_lbs"]
        )
        daily["same_window_delta"] = (
            daily["fat_mass_lbs"] - daily["fat_mass_lbs"].shift(horizon)
        )

        for outcome in ["future_delta", "same_window_delta"]:
            d = daily[[std_col, surplus_col, outcome]].dropna()
            controls = np.column_stack(
                [d[surplus_col].to_numpy(float), np.ones(len(d))]
            )
            x = d[std_col].to_numpy(float)
            y = d[outcome].to_numpy(float)
            x_resid = x - controls @ np.linalg.lstsq(
                controls, x, rcond=None
            )[0]
            y_resid = y - controls @ np.linalg.lstsq(
                controls, y, rcond=None
            )[0]
            partial_r = float(np.corrcoef(x_resid, y_resid)[0, 1])
            regression = np.column_stack(
                [d[surplus_col].to_numpy(float), x, np.ones(len(d))]
            )
            coefficient = float(
                np.linalg.lstsq(regression, y, rcond=None)[0][1]
            )
            add_check(
                rows,
                "AF_variance",
                f"h{horizon}_{outcome}_partial_r",
                round(partial_r, 6),
                f"n={len(d)}; controls trailing modeled surplus.",
            )
            add_check(
                rows,
                "AF_variance",
                f"h{horizon}_{outcome}_std_coefficient_lb_per_kcal_sd",
                round(coefficient, 9),
                "Legacy P4 latent FM/TDEE; comparison diagnoses temporal alignment only.",
            )


def window_identity_checks(rows: list[dict[str, object]]) -> None:
    identities = [
        (
            "shift(-1).rolling(7)",
            "t-5 through t+1",
            "5 prior days, current day, 1 future day",
        ),
        (
            "shift(-1).rolling(14)",
            "t-12 through t+1",
            "12 prior days, current day, 1 future day",
        ),
        (
            "rolling(3).shift(-1)",
            "t-1 through t+1",
            "1 prior day, current day, 1 future day",
        ),
        (
            "rolling(7).shift(-7)",
            "t+1 through t+7",
            "Correct following-seven-day window on a complete daily index",
        ),
    ]
    for expression, span, notes in identities:
        add_check(rows, "window_algebra", expression, span, notes)


def build_checks() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    window_identity_checks(rows)
    gap_checks(rows)
    n_gravitostat_check(rows)
    ab_running_check(rows)
    af_variance_check(rows)
    return pd.DataFrame(rows)


def main() -> None:
    inventory = build_inventory()
    checks = build_checks()
    inventory.to_csv(INVENTORY_PATH, index=False)
    checks.to_csv(CHECKS_PATH, index=False)

    print("Temporal semantics audit")
    print(f"  inventory rows: {len(inventory)}")
    print(inventory.groupby("severity").size().to_string())
    print(f"  quantitative checks: {len(checks)}")
    print(f"  wrote {INVENTORY_PATH}")
    print(f"  wrote {CHECKS_PATH}")


if __name__ == "__main__":
    main()
