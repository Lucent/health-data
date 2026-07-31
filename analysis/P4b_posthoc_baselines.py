"""Post-hoc persistence baselines for the frozen P4b v1 outputs.

This script does not tune or modify P4b. The 2025-2026 interval has already
been observed, so these comparisons are descriptive post-hoc benchmarks, not
a second locked validation.

Run:
    uv run --with numpy --with pandas python analysis/P4b_posthoc_baselines.py
"""

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
STATE_PATH = ROOT / "analysis" / "P4b_state_daily.csv"
OUTPUT_PATH = ROOT / "analysis" / "P4b_posthoc_baselines.csv"


def persistence_predictions(values):
    """Return the most recent strictly prior observation and its age."""
    predictions = np.full(len(values), np.nan)
    ages = np.full(len(values), np.nan)
    last_value = np.nan
    last_index = None
    for index, value in enumerate(values):
        if np.isfinite(last_value):
            predictions[index] = last_value
            ages[index] = index - last_index
        if np.isfinite(value):
            last_value = value
            last_index = index
    return predictions, ages


def fm_weight_adjusted_predictions(fm_values, total_weight_values):
    """Carry prior FM forward plus observed total-weight change since that scan."""
    predictions = np.full(len(fm_values), np.nan)
    last_fm = np.nan
    last_total_weight = np.nan
    for index, (fm, total_weight) in enumerate(
        zip(fm_values, total_weight_values)
    ):
        if np.isfinite(last_fm) and np.isfinite(last_total_weight) and np.isfinite(
            total_weight
        ):
            predictions[index] = last_fm + total_weight - last_total_weight
        if np.isfinite(fm) and np.isfinite(total_weight):
            last_fm = fm
            last_total_weight = total_weight
    return predictions


def metrics(observed, predicted):
    mask = np.isfinite(observed) & np.isfinite(predicted)
    residual = observed[mask] - predicted[mask]
    return {
        "n": int(mask.sum()),
        "mae": float(np.mean(np.abs(residual))),
        "rmse": float(np.sqrt(np.mean(residual**2))),
    }


def main():
    frame = pd.read_csv(STATE_PATH, parse_dates=["date"])

    total_observation = frame["scale_weight_lbs_observed"].to_numpy(
        dtype=float, copy=True
    )
    missing_scale = ~np.isfinite(total_observation)
    total_observation[missing_scale] = frame.loc[
        missing_scale, "composition_weight_lbs_observed"
    ].to_numpy(dtype=float)

    fm_observation = frame["composition_fm_lbs_observed"].to_numpy(dtype=float)
    weight_persistence, weight_age = persistence_predictions(total_observation)
    fm_persistence, fm_age = persistence_predictions(fm_observation)
    fm_weight_adjusted = fm_weight_adjusted_predictions(
        fm_observation, total_observation
    )
    fm_conditional_prediction = (
        fm_observation
        - frame["composition_fm_innovation_lbs"].to_numpy(dtype=float)
    )

    comparisons = (
        (
            "total_weight",
            total_observation,
            frame["predicted_weight_lbs"].to_numpy(dtype=float),
            weight_persistence,
            weight_age,
            None,
        ),
        (
            "composition_fm",
            fm_observation,
            fm_conditional_prediction,
            fm_persistence,
            fm_age,
            fm_weight_adjusted,
        ),
    )

    rows = []
    initialization = (
        frame["update_source"].astype(str).str.startswith("initial_").to_numpy()
    )
    for segment in ("development", "calibration", "locked_test", "all"):
        segment_mask = (
            np.ones(len(frame), dtype=bool)
            if segment == "all"
            else frame["segment"].eq(segment).to_numpy()
        )
        for (
            channel,
            observed,
            p4b_prediction,
            persistence,
            age,
            weight_adjusted,
        ) in comparisons:
            common_mask = (
                segment_mask
                & ~initialization
                & np.isfinite(observed)
                & np.isfinite(p4b_prediction)
                & np.isfinite(persistence)
            )
            common_observed = np.where(common_mask, observed, np.nan)
            p4b = metrics(common_observed, p4b_prediction)
            naive = metrics(common_observed, persistence)
            adjusted = (
                metrics(common_observed, weight_adjusted)
                if weight_adjusted is not None
                else {"n": 0, "mae": np.nan, "rmse": np.nan}
            )
            age_mask = common_mask & np.isfinite(age)
            rows.append(
                {
                    "segment": segment,
                    "channel": channel,
                    "comparison_status": "posthoc_after_locked_test_opened",
                    "p4b_n": p4b["n"],
                    "p4b_mae": p4b["mae"],
                    "p4b_rmse": p4b["rmse"],
                    "persistence_n": naive["n"],
                    "persistence_mae": naive["mae"],
                    "persistence_rmse": naive["rmse"],
                    "weight_adjusted_persistence_n": adjusted["n"],
                    "weight_adjusted_persistence_mae": adjusted["mae"],
                    "weight_adjusted_persistence_rmse": adjusted["rmse"],
                    "p4b_mae_improvement_vs_persistence": (
                        naive["mae"] - p4b["mae"]
                    ),
                    "p4b_mae_improvement_vs_weight_adjusted_persistence": (
                        adjusted["mae"] - p4b["mae"]
                        if np.isfinite(adjusted["mae"])
                        else np.nan
                    ),
                    "median_days_since_persistence_anchor": (
                        float(np.median(age[age_mask]))
                        if age_mask.any()
                        else np.nan
                    ),
                }
            )

    output = pd.DataFrame(rows)
    output.to_csv(OUTPUT_PATH, index=False, float_format="%.8f")
    print(output.to_string(index=False))
    print(f"\nWrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
