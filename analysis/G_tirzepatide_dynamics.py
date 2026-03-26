"""Analyze tirzepatide as a transition modifier, not just a mean-intake shift.

Questions:
1. Does the drug reduce escalation from normal/high days into binges?
2. Does it reduce persistence of high-intake streaks?
3. Does it suppress rebound after multi-day restriction runs?
4. Where does the drug fail: which high-drug days still become high intake?

Outputs:
    analysis/tirzepatide_transition_summary.csv
    analysis/tirzepatide_rebound_summary.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

BINGE_THRESHOLD = 2800
HIGH_THRESHOLD = 2400
RESTRICTION_THRESHOLD = 1800
RUN_MIN_DAYS = 3

STATE_ORDER = ["restriction", "typical", "high", "binge"]


def load_data():
    intake = pd.read_csv(ROOT / "intake" / "intake_daily.csv", parse_dates=["date"])
    tirz = pd.read_csv(ROOT / "drugs" / "tirzepatide.csv", parse_dates=["date"])
    kalman = pd.read_csv(ROOT / "analysis" / "P4_kalman_daily.csv", parse_dates=["date"])

    daily = intake.merge(
        tirz[["date", "effective_level", "dose_mg", "days_since_injection"]],
        on="date",
        how="left",
    )
    daily = daily.merge(
        kalman[["date", "tdee_filtered", "fat_mass_lbs_filtered"]],
        on="date",
        how="left",
    )
    daily["effective_level"] = daily["effective_level"].fillna(0)
    daily["on_tirz"] = daily["effective_level"] > 0
    return daily.sort_values("date").reset_index(drop=True)


def classify_state(calories):
    if pd.isna(calories):
        return np.nan
    if calories < RESTRICTION_THRESHOLD:
        return "restriction"
    if calories < HIGH_THRESHOLD:
        return "typical"
    if calories < BINGE_THRESHOLD:
        return "high"
    return "binge"


def add_states(daily):
    d = daily.copy()
    d["state"] = d["calories"].map(classify_state)
    d["next_state"] = d["state"].shift(-1)
    d["next_calories"] = d["calories"].shift(-1)
    d["next_binge"] = (d["next_calories"] > BINGE_THRESHOLD).astype(float)
    d["next_high_or_binge"] = (d["next_calories"] >= HIGH_THRESHOLD).astype(float)
    d["prior_3d_mean"] = d["calories"].rolling(3, min_periods=3).mean().shift(1)
    d["high_drug"] = d["effective_level"] >= d.loc[d["on_tirz"], "effective_level"].median()
    return d


def transition_summary(daily):
    rows = []
    valid = daily.dropna(subset=["state", "next_state"]).copy()

    for cohort_name, cohort_mask in [
        ("pre_tirzepatide", valid["on_tirz"] == 0),
        ("on_tirzepatide", valid["on_tirz"] == 1),
    ]:
        cohort = valid[cohort_mask]
        for state in STATE_ORDER:
            sub = cohort[cohort["state"] == state]
            row = {
                "cohort": cohort_name,
                "state": state,
                "n_days": len(sub),
                "next_day_calories": round(sub["next_calories"].mean(), 1),
                "p_next_high_or_binge": round(sub["next_high_or_binge"].mean(), 4),
                "p_next_binge": round(sub["next_binge"].mean(), 4),
            }
            for next_state in STATE_ORDER:
                row[f"to_{next_state}"] = round((sub["next_state"] == next_state).mean(), 4)
            rows.append(row)

    return pd.DataFrame(rows)


def detect_restriction_runs(daily):
    d = daily.copy()
    d["is_restriction"] = d["calories"] < RESTRICTION_THRESHOLD
    groups = (d["is_restriction"] != d["is_restriction"].shift()).cumsum()
    runs = []

    for _, grp in d.groupby(groups):
        if not bool(grp["is_restriction"].iloc[0]):
            continue
        if len(grp) < RUN_MIN_DAYS:
            continue
        end_idx = grp.index[-1]
        if end_idx + 7 >= len(d):
            continue
        future = d.loc[end_idx + 1:end_idx + 7]
        run = {
            "end_date": d.loc[end_idx, "date"],
            "run_days": len(grp),
            "run_mean_calories": grp["calories"].mean(),
            "cohort": "on_tirzepatide" if d.loc[end_idx, "on_tirz"] else "pre_tirzepatide",
            "effective_level_end": d.loc[end_idx, "effective_level"],
            "dose_mg_end": d.loc[end_idx, "dose_mg"],
            "days_since_injection_end": d.loc[end_idx, "days_since_injection"],
            "next7_mean_calories": future["calories"].mean(),
            "next7_max_calories": future["calories"].max(),
            "next7_binge": int((future["calories"] > BINGE_THRESHOLD).any()),
            "next7_high_or_binge": int((future["calories"] >= HIGH_THRESHOLD).any()),
        }
        runs.append(run)

    return pd.DataFrame(runs)


def summarize_restriction_runs(runs):
    rows = []
    if runs.empty:
        return pd.DataFrame(rows)

    for cohort, grp in runs.groupby("cohort"):
        rows.append(
            {
                "cohort": cohort,
                "n_runs": len(grp),
                "mean_run_days": round(grp["run_days"].mean(), 2),
                "mean_run_calories": round(grp["run_mean_calories"].mean(), 1),
                "next7_mean_calories": round(grp["next7_mean_calories"].mean(), 1),
                "next7_binge_rate": round(grp["next7_binge"].mean(), 4),
                "next7_high_or_binge_rate": round(grp["next7_high_or_binge"].mean(), 4),
            }
        )
    return pd.DataFrame(rows)


def failure_modes(daily):
    on = daily[daily["on_tirz"]].copy()
    if on.empty:
        return pd.DataFrame()
    fail = on[(on["high_drug"]) & (on["calories"] >= HIGH_THRESHOLD)].copy()
    cols = ["date", "calories", "state", "effective_level", "dose_mg", "days_since_injection", "prior_3d_mean"]
    return fail.sort_values(["calories", "effective_level"], ascending=[False, False])[cols].head(15)


def print_report(transitions, rebound, failures):
    print("\n=== Tirzepatide as Transition Modifier ===")
    for state in ["typical", "high", "binge", "restriction"]:
        pre = transitions[(transitions["cohort"] == "pre_tirzepatide") & (transitions["state"] == state)]
        post = transitions[(transitions["cohort"] == "on_tirzepatide") & (transitions["state"] == state)]
        if pre.empty or post.empty:
            continue
        pre_row = pre.iloc[0]
        post_row = post.iloc[0]
        print(
            f"{state:>11}: next binge {pre_row['p_next_binge']*100:5.1f}% -> "
            f"{post_row['p_next_binge']*100:5.1f}% | next high/binge "
            f"{pre_row['p_next_high_or_binge']*100:5.1f}% -> "
            f"{post_row['p_next_high_or_binge']*100:5.1f}%"
        )

    if not rebound.empty:
        print("\n=== Post-restriction rebound (next 7 days after >=3-day restriction runs) ===")
        for _, row in rebound.iterrows():
            print(
                f"{row['cohort']:>15}: runs={int(row['n_runs'])}  next7 binge="
                f"{row['next7_binge_rate']*100:5.1f}%  next7 high/binge="
                f"{row['next7_high_or_binge_rate']*100:5.1f}%  "
                f"next7 mean cal={row['next7_mean_calories']:.0f}"
            )

    if not failures.empty:
        print("\n=== Highest-intake high-drug days (failure modes) ===")
        for _, row in failures.head(8).iterrows():
            print(
                f"{row['date'].strftime('%Y-%m-%d')}: {int(row['calories'])} cal  "
                f"level={row['effective_level']:.2f}  dose={row['dose_mg']:.1f}  "
                f"day={int(row['days_since_injection']) if not pd.isna(row['days_since_injection']) else 'NA'}"
            )


def save_outputs(transitions, rebound):
    transitions.to_csv(ROOT / "analysis" / "G_tirzepatide_transition_summary.csv", index=False)
    rebound.to_csv(ROOT / "analysis" / "G_tirzepatide_rebound_summary.csv", index=False)


def main():
    daily = add_states(load_data())
    transitions = transition_summary(daily)
    runs = detect_restriction_runs(daily)
    rebound = summarize_restriction_runs(runs)
    failures = failure_modes(daily)
    save_outputs(transitions, rebound)
    print_report(transitions, rebound, failures)


if __name__ == "__main__":
    main()
