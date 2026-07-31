"""Leakage-resistant forward state-space model for body composition and expenditure.

This is an opt-in replacement candidate for P4. It deliberately does not read
P1, P2, P3, P4, RMR, temperature, or any downstream analysis artifact.

Primary model:
    state = [fat mass, fat-free mass, transient mass buffer,
             expenditure conditional on logged intake]

The primary intake multiplier is 1.0. Consequently, the expenditure state is
on the logged-intake scale and must not be presented as independently measured
physiological TDEE. Other intake multipliers belong in explicit sensitivity
analyses, not in model selection.

Outputs:
    analysis/P4b_state_daily.csv
    analysis/P4b_covariance_daily.csv
    analysis/P4b_temporal_validation.csv
    analysis/P4b_config_selection.csv

Run:
    uv run --with numpy --with pandas python analysis/P4b_forward_state_space.py
"""

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
STATE_OUTPUT = ROOT / "analysis" / "P4b_state_daily.csv"
COVARIANCE_OUTPUT = ROOT / "analysis" / "P4b_covariance_daily.csv"
VALIDATION_OUTPUT = ROOT / "analysis" / "P4b_temporal_validation.csv"
CONFIG_OUTPUT = ROOT / "analysis" / "P4b_config_selection.csv"

RAW_INPUT_PATHS = {
    "intake": ROOT / "intake" / "intake_daily.csv",
    "weight": ROOT / "weight" / "weight.csv",
    "composition": ROOT / "composition" / "composition.csv",
}

STATE_NAMES = ("fm_lbs", "ffm_lbs", "buffer_lbs", "tdee_logged_kcal")
STATE_DIM = len(STATE_NAMES)
UPPER_TRIANGLE = tuple(
    (i, j) for i in range(STATE_DIM) for j in range(i, STATE_DIM)
)

E_FAT_KCAL_PER_LB = 3500.0
E_FFM_KCAL_PER_LB = 1100.0
FORBES_C_KG = 10.4
LBS_TO_KG = 0.453592

FM_GAIN_LIMIT_LBS_PER_DAY = 0.5
FM_LOSS_LIMIT_PER_LB_FM = 31.0 / E_FAT_KCAL_PER_LB
FFM_GAIN_LIMIT_LBS_PER_DAY = 0.15
FFM_LOSS_LIMIT_LBS_PER_DAY = 0.5

# Fixed analysis assumptions. These are not estimates of device truth.
# They are intentionally conservative and must be varied in sensitivity work.
COMPOSITION_ERROR_SD_LBS = {
    "BOD POD": {"weight": 1.50, "fm": 2.50},
    "370": {"weight": 0.75, "fm": 3.00},
    "570": {"weight": 0.75, "fm": 3.50},
    "H20": {"weight": 1.00, "fm": 4.00},
    "H20N": {"weight": 1.00, "fm": 4.00},
    "H30": {"weight": 0.75, "fm": 4.00},
}
UNKNOWN_DEVICE_ERROR_SD_LBS = {"weight": 1.50, "fm": 4.50}


@dataclass(frozen=True)
class ModelConfig:
    intake_multiplier: float = 1.0
    initial_tdee_kcal: float = 2200.0
    initial_tdee_sd_kcal: float = 600.0
    initial_buffer_sd_lbs: float = 2.0
    buffer_half_life_days: float = 5.0
    scale_weight_sd_lbs: float = 1.50
    process_fm_sd_lbs: float = np.sqrt(0.005)
    process_ffm_sd_lbs: float = np.sqrt(0.005)
    process_buffer_sd_lbs: float = np.sqrt(0.32)
    process_tdee_sd_kcal: float = 25.0


@dataclass(frozen=True)
class RawInputs:
    dates: np.ndarray
    calories_logged: np.ndarray
    scale_weight_lbs: np.ndarray
    composition_weight_lbs: np.ndarray
    composition_fm_lbs: np.ndarray
    composition_device: np.ndarray
    segment: np.ndarray
    days_since_prior_scale_weight: np.ndarray


@dataclass(frozen=True)
class FilterResult:
    x_predicted: np.ndarray
    p_predicted: np.ndarray
    x_filtered: np.ndarray
    p_filtered: np.ndarray
    transition_jacobian: np.ndarray
    weight_innovation: np.ndarray
    weight_innovation_variance: np.ndarray
    composition_fm_innovation: np.ndarray
    composition_fm_innovation_variance: np.ndarray
    composition_joint_nis: np.ndarray
    scale_weight_used: np.ndarray
    composition_weight_used: np.ndarray
    composition_fm_used: np.ndarray
    update_source: np.ndarray
    device_error_fallback: np.ndarray
    covariance_repaired_predicted: np.ndarray
    covariance_repaired_filtered: np.ndarray
    first_state_index: int


@dataclass(frozen=True)
class SmoothResult:
    x_smoothed: np.ndarray
    p_smoothed: np.ndarray
    covariance_repaired_smoothed: np.ndarray


def _validate_unique_dates(frame, label):
    duplicate = frame["date"].duplicated(keep=False)
    if duplicate.any():
        dates = frame.loc[duplicate, "date"].dt.strftime("%Y-%m-%d").tolist()
        raise ValueError(f"{label} contains duplicate dates: {dates[:10]}")


def _temporal_segment(dates):
    timestamp = pd.to_datetime(dates)
    return np.where(
        timestamp <= pd.Timestamp("2023-12-31"),
        "development",
        np.where(
            timestamp <= pd.Timestamp("2024-12-31"),
            "calibration",
            "locked_test",
        ),
    )


def _days_since_prior_observation(values):
    result = np.full(len(values), np.nan)
    last_index = None
    for i, value in enumerate(values):
        if last_index is not None:
            result[i] = i - last_index
        if np.isfinite(value):
            last_index = i
    return result


def load_raw_inputs():
    """Load only the three canonical raw/minimally transformed CSVs."""
    intake = pd.read_csv(RAW_INPUT_PATHS["intake"], parse_dates=["date"])
    weight = pd.read_csv(RAW_INPUT_PATHS["weight"], parse_dates=["date"])
    composition = pd.read_csv(
        RAW_INPUT_PATHS["composition"], parse_dates=["date"]
    )

    intake = intake[["date", "calories"]].sort_values("date").reset_index(drop=True)
    weight = (
        weight[["date", "weight_lbs"]]
        .dropna(subset=["weight_lbs"])
        .sort_values("date")
        .reset_index(drop=True)
    )
    composition = (
        composition[
            ["date", "device", "weight_lbs", "fat_mass_lbs"]
        ]
        .dropna(subset=["weight_lbs", "fat_mass_lbs"])
        .sort_values("date")
        .reset_index(drop=True)
    )

    _validate_unique_dates(intake, "intake")
    _validate_unique_dates(weight, "weight")
    _validate_unique_dates(composition, "composition")

    expected_dates = pd.date_range(
        intake["date"].min(), intake["date"].max(), freq="D"
    )
    if len(expected_dates) != len(intake) or not np.array_equal(
        expected_dates.values, intake["date"].values
    ):
        raise ValueError("intake dates are not a complete consecutive daily series")
    if intake["calories"].isna().any():
        raise ValueError("intake calories contain missing values")

    daily = intake.merge(weight, on="date", how="left", validate="one_to_one")
    daily = daily.merge(
        composition.rename(
            columns={
                "device": "composition_device",
                "weight_lbs": "composition_weight_lbs",
                "fat_mass_lbs": "composition_fm_lbs",
            }
        ),
        on="date",
        how="left",
        validate="one_to_one",
    )

    scale_weight = daily["weight_lbs"].to_numpy(dtype=float)
    return RawInputs(
        dates=daily["date"].to_numpy(),
        calories_logged=daily["calories"].to_numpy(dtype=float),
        scale_weight_lbs=scale_weight,
        composition_weight_lbs=daily["composition_weight_lbs"].to_numpy(dtype=float),
        composition_fm_lbs=daily["composition_fm_lbs"].to_numpy(dtype=float),
        composition_device=daily["composition_device"].fillna("").to_numpy(dtype=object),
        segment=_temporal_segment(daily["date"]),
        days_since_prior_scale_weight=_days_since_prior_observation(scale_weight),
    )


def forbes_fat_fraction(fm_lbs):
    fm_kg = max(float(fm_lbs), 5.0) * LBS_TO_KG
    fraction = 1.0 - FORBES_C_KG**2 / (FORBES_C_KG + fm_kg) ** 2
    return float(np.clip(fraction, 0.0, 1.0))


def transition_state(state, calories_logged, config):
    """Advance morning state t to morning state t+1 using intake on day t."""
    fm, ffm, buffer_lbs, tdee_logged = state
    fat_fraction = forbes_fat_fraction(fm)
    surplus = calories_logged * config.intake_multiplier - tdee_logged
    energy_per_lb = (
        fat_fraction * E_FAT_KCAL_PER_LB
        + (1.0 - fat_fraction) * E_FFM_KCAL_PER_LB
    )
    delta_weight = surplus / energy_per_lb
    delta_fm_raw = fat_fraction * delta_weight
    delta_ffm_raw = (1.0 - fat_fraction) * delta_weight

    delta_fm = np.clip(
        delta_fm_raw,
        -max(fm, 5.0) * FM_LOSS_LIMIT_PER_LB_FM,
        FM_GAIN_LIMIT_LBS_PER_DAY,
    )
    delta_ffm = np.clip(
        delta_ffm_raw,
        -FFM_LOSS_LIMIT_LBS_PER_DAY,
        FFM_GAIN_LIMIT_LBS_PER_DAY,
    )
    buffer_decay = np.exp(-np.log(2.0) / config.buffer_half_life_days)
    return np.array(
        [
            fm + delta_fm,
            ffm + delta_ffm,
            buffer_decay * buffer_lbs,
            tdee_logged,
        ],
        dtype=float,
    )


def transition_jacobian(state, calories_logged, config):
    """Deterministic central-difference Jacobian of the actual transition."""
    steps = np.array([1e-4, 1e-4, 1e-4, 1e-2], dtype=float)
    jacobian = np.empty((STATE_DIM, STATE_DIM), dtype=float)
    for j, step in enumerate(steps):
        plus = state.copy()
        minus = state.copy()
        plus[j] += step
        minus[j] -= step
        jacobian[:, j] = (
            transition_state(plus, calories_logged, config)
            - transition_state(minus, calories_logged, config)
        ) / (2.0 * step)
    return jacobian


def process_covariance(config):
    return np.diag(
        [
            config.process_fm_sd_lbs**2,
            config.process_ffm_sd_lbs**2,
            config.process_buffer_sd_lbs**2,
            config.process_tdee_sd_kcal**2,
        ]
    )


def stabilize_covariance(covariance):
    """Symmetrize and repair only numerical loss of positive semidefiniteness."""
    symmetric = (covariance + covariance.T) / 2.0
    values, vectors = np.linalg.eigh(symmetric)
    repaired = bool(values.min() < -1e-10)
    if values.min() < 1e-12:
        values = np.maximum(values, 1e-12)
        symmetric = vectors @ np.diag(values) @ vectors.T
        symmetric = (symmetric + symmetric.T) / 2.0
    return symmetric, repaired


def joseph_update(state, covariance, observation, design, error_covariance):
    """Linear observation update with Joseph-form posterior covariance."""
    observation = np.atleast_1d(np.asarray(observation, dtype=float))
    design = np.atleast_2d(np.asarray(design, dtype=float))
    error_covariance = np.atleast_2d(
        np.asarray(error_covariance, dtype=float)
    )

    innovation = observation - design @ state
    innovation_covariance = (
        design @ covariance @ design.T + error_covariance
    )
    cross_covariance = covariance @ design.T
    gain = np.linalg.solve(
        innovation_covariance.T, cross_covariance.T
    ).T
    posterior_state = state + gain @ innovation
    identity_minus_kh = np.eye(STATE_DIM) - gain @ design
    posterior_covariance = (
        identity_minus_kh @ covariance @ identity_minus_kh.T
        + gain @ error_covariance @ gain.T
    )
    posterior_covariance, repaired = stabilize_covariance(
        posterior_covariance
    )
    nis = float(
        innovation.T
        @ np.linalg.solve(innovation_covariance, innovation)
    )
    return (
        posterior_state,
        posterior_covariance,
        innovation,
        innovation_covariance,
        nis,
        repaired,
    )


def device_error_covariance(device):
    values = COMPOSITION_ERROR_SD_LBS.get(str(device))
    fallback = values is None
    if fallback:
        values = UNKNOWN_DEVICE_ERROR_SD_LBS
    covariance = np.diag([values["weight"] ** 2, values["fm"] ** 2])
    return covariance, fallback


def _initial_state_and_covariance(inputs, first_index, config):
    device = inputs.composition_device[first_index]
    device_covariance, fallback = device_error_covariance(device)
    fm = inputs.composition_fm_lbs[first_index]

    has_scale = np.isfinite(inputs.scale_weight_lbs[first_index])
    if has_scale:
        total_weight = inputs.scale_weight_lbs[first_index]
        weight_variance = config.scale_weight_sd_lbs**2
        update_source = "initial_scale+composition_fm"
    else:
        total_weight = inputs.composition_weight_lbs[first_index]
        weight_variance = device_covariance[0, 0]
        update_source = "initial_composition_weight+fm"

    fm_variance = device_covariance[1, 1]
    buffer_variance = config.initial_buffer_sd_lbs**2
    state = np.array(
        [fm, total_weight - fm, 0.0, config.initial_tdee_kcal],
        dtype=float,
    )

    # Transform independent uncertainties in total weight, FM, and buffer
    # through FFM = total weight - FM - buffer.
    covariance = np.zeros((STATE_DIM, STATE_DIM), dtype=float)
    covariance[0, 0] = fm_variance
    covariance[1, 1] = weight_variance + fm_variance + buffer_variance
    covariance[2, 2] = buffer_variance
    covariance[3, 3] = config.initial_tdee_sd_kcal**2
    covariance[0, 1] = covariance[1, 0] = -fm_variance
    covariance[1, 2] = covariance[2, 1] = -buffer_variance
    covariance, repaired = stabilize_covariance(covariance)
    return state, covariance, has_scale, update_source, fallback, repaired


def run_forward_filter(inputs, config=ModelConfig()):
    n = len(inputs.dates)
    x_predicted = np.full((n, STATE_DIM), np.nan)
    p_predicted = np.full((n, STATE_DIM, STATE_DIM), np.nan)
    x_filtered = np.full((n, STATE_DIM), np.nan)
    p_filtered = np.full((n, STATE_DIM, STATE_DIM), np.nan)
    jacobian = np.full((n, STATE_DIM, STATE_DIM), np.nan)

    weight_innovation = np.full(n, np.nan)
    weight_innovation_variance = np.full(n, np.nan)
    composition_fm_innovation = np.full(n, np.nan)
    composition_fm_innovation_variance = np.full(n, np.nan)
    composition_joint_nis = np.full(n, np.nan)

    scale_weight_used = np.zeros(n, dtype=bool)
    composition_weight_used = np.zeros(n, dtype=bool)
    composition_fm_used = np.zeros(n, dtype=bool)
    update_source = np.full(n, "none", dtype=object)
    device_error_fallback = np.zeros(n, dtype=bool)
    covariance_repaired_predicted = np.zeros(n, dtype=bool)
    covariance_repaired_filtered = np.zeros(n, dtype=bool)

    composition_indices = np.where(np.isfinite(inputs.composition_fm_lbs))[0]
    if not len(composition_indices):
        raise ValueError("at least one composition scan is required")
    first_index = int(composition_indices[0])

    (
        state,
        covariance,
        initial_scale_used,
        initial_source,
        initial_fallback,
        initial_repaired,
    ) = _initial_state_and_covariance(inputs, first_index, config)
    x_predicted[first_index] = state
    p_predicted[first_index] = covariance
    x_filtered[first_index] = state
    p_filtered[first_index] = covariance
    jacobian[first_index] = np.eye(STATE_DIM)
    scale_weight_used[first_index] = initial_scale_used
    composition_weight_used[first_index] = not initial_scale_used
    composition_fm_used[first_index] = True
    update_source[first_index] = initial_source
    device_error_fallback[first_index] = initial_fallback
    covariance_repaired_filtered[first_index] = initial_repaired

    q_covariance = process_covariance(config)
    weight_design = np.array([[1.0, 1.0, 1.0, 0.0]])
    fm_design = np.array([[1.0, 0.0, 0.0, 0.0]])
    joint_design = np.array(
        [
            [1.0, 1.0, 1.0, 0.0],
            [1.0, 0.0, 0.0, 0.0],
        ]
    )

    for i in range(first_index + 1, n):
        previous_state = x_filtered[i - 1]
        previous_covariance = p_filtered[i - 1]
        f_matrix = transition_jacobian(
            previous_state, inputs.calories_logged[i - 1], config
        )
        state = transition_state(
            previous_state, inputs.calories_logged[i - 1], config
        )
        covariance = (
            f_matrix @ previous_covariance @ f_matrix.T + q_covariance
        )
        covariance, repaired = stabilize_covariance(covariance)

        jacobian[i] = f_matrix
        x_predicted[i] = state
        p_predicted[i] = covariance
        covariance_repaired_predicted[i] = repaired

        has_scale = np.isfinite(inputs.scale_weight_lbs[i])
        has_composition = np.isfinite(inputs.composition_fm_lbs[i])

        if has_scale:
            (
                state,
                covariance,
                innovation,
                innovation_covariance,
                _,
                repaired,
            ) = joseph_update(
                state,
                covariance,
                [inputs.scale_weight_lbs[i]],
                weight_design,
                [[config.scale_weight_sd_lbs**2]],
            )
            weight_innovation[i] = innovation[0]
            weight_innovation_variance[i] = innovation_covariance[0, 0]
            scale_weight_used[i] = True
            update_source[i] = "scale_weight"
            covariance_repaired_filtered[i] |= repaired

        if has_composition:
            device_covariance, fallback = device_error_covariance(
                inputs.composition_device[i]
            )
            device_error_fallback[i] = fallback
            composition_fm_used[i] = True

            if has_scale:
                # The scan's total weight is not assimilated because the same-day
                # scale weight has already provided the total-mass observation.
                (
                    state,
                    covariance,
                    innovation,
                    innovation_covariance,
                    _,
                    repaired,
                ) = joseph_update(
                    state,
                    covariance,
                    [inputs.composition_fm_lbs[i]],
                    fm_design,
                    [[device_covariance[1, 1]]],
                )
                composition_fm_innovation[i] = innovation[0]
                composition_fm_innovation_variance[i] = (
                    innovation_covariance[0, 0]
                )
                update_source[i] = "scale_weight+composition_fm"
                covariance_repaired_filtered[i] |= repaired
            else:
                (
                    state,
                    covariance,
                    innovation,
                    innovation_covariance,
                    joint_nis,
                    repaired,
                ) = joseph_update(
                    state,
                    covariance,
                    [
                        inputs.composition_weight_lbs[i],
                        inputs.composition_fm_lbs[i],
                    ],
                    joint_design,
                    device_covariance,
                )
                weight_innovation[i] = innovation[0]
                weight_innovation_variance[i] = innovation_covariance[0, 0]
                composition_fm_innovation[i] = innovation[1]
                composition_fm_innovation_variance[i] = (
                    innovation_covariance[1, 1]
                )
                composition_joint_nis[i] = joint_nis
                composition_weight_used[i] = True
                update_source[i] = "composition_weight+fm"
                covariance_repaired_filtered[i] |= repaired

        x_filtered[i] = state
        p_filtered[i] = covariance

    return FilterResult(
        x_predicted=x_predicted,
        p_predicted=p_predicted,
        x_filtered=x_filtered,
        p_filtered=p_filtered,
        transition_jacobian=jacobian,
        weight_innovation=weight_innovation,
        weight_innovation_variance=weight_innovation_variance,
        composition_fm_innovation=composition_fm_innovation,
        composition_fm_innovation_variance=composition_fm_innovation_variance,
        composition_joint_nis=composition_joint_nis,
        scale_weight_used=scale_weight_used,
        composition_weight_used=composition_weight_used,
        composition_fm_used=composition_fm_used,
        update_source=update_source,
        device_error_fallback=device_error_fallback,
        covariance_repaired_predicted=covariance_repaired_predicted,
        covariance_repaired_filtered=covariance_repaired_filtered,
        first_state_index=first_index,
    )


def run_rts_smoother(result):
    x_smoothed = result.x_filtered.copy()
    p_smoothed = result.p_filtered.copy()
    repaired_flags = np.zeros(len(x_smoothed), dtype=bool)
    first_index = result.first_state_index
    valid_indices = np.where(np.isfinite(result.x_filtered[:, 0]))[0]
    last_index = int(valid_indices[-1])

    for i in range(last_index - 1, first_index - 1, -1):
        f_next = result.transition_jacobian[i + 1]
        cross_covariance = result.p_filtered[i] @ f_next.T
        smoothing_gain = np.linalg.solve(
            result.p_predicted[i + 1].T,
            cross_covariance.T,
        ).T
        x_smoothed[i] = (
            result.x_filtered[i]
            + smoothing_gain
            @ (x_smoothed[i + 1] - result.x_predicted[i + 1])
        )
        covariance = (
            result.p_filtered[i]
            + smoothing_gain
            @ (p_smoothed[i + 1] - result.p_predicted[i + 1])
            @ smoothing_gain.T
        )
        covariance, repaired = stabilize_covariance(covariance)
        p_smoothed[i] = covariance
        repaired_flags[i] = repaired

    return SmoothResult(
        x_smoothed=x_smoothed,
        p_smoothed=p_smoothed,
        covariance_repaired_smoothed=repaired_flags,
    )


def _state_standard_deviations(covariances):
    diagonal = np.diagonal(covariances, axis1=1, axis2=2)
    return np.sqrt(np.maximum(diagonal, 0.0))


def build_state_output(inputs, result, smooth, config):
    frame = pd.DataFrame(
        {
            "date": pd.to_datetime(inputs.dates),
            "segment": inputs.segment,
            "calories_logged": inputs.calories_logged,
            "intake_multiplier_assumed": config.intake_multiplier,
            "scale_weight_lbs_observed": inputs.scale_weight_lbs,
            "composition_weight_lbs_observed": inputs.composition_weight_lbs,
            "composition_fm_lbs_observed": inputs.composition_fm_lbs,
            "composition_device": inputs.composition_device,
            "days_since_prior_scale_weight": inputs.days_since_prior_scale_weight,
            "scale_weight_used": result.scale_weight_used,
            "composition_weight_used": result.composition_weight_used,
            "composition_fm_used": result.composition_fm_used,
            "update_source": result.update_source,
            "device_error_fallback": result.device_error_fallback,
            "weight_innovation_lbs": result.weight_innovation,
            "weight_innovation_variance": result.weight_innovation_variance,
            "weight_standardized_innovation": (
                result.weight_innovation
                / np.sqrt(result.weight_innovation_variance)
            ),
            "weight_nis_component": (
                result.weight_innovation**2
                / result.weight_innovation_variance
            ),
            "composition_fm_innovation_lbs": result.composition_fm_innovation,
            "composition_fm_innovation_variance": (
                result.composition_fm_innovation_variance
            ),
            "composition_fm_standardized_innovation": (
                result.composition_fm_innovation
                / np.sqrt(result.composition_fm_innovation_variance)
            ),
            "composition_fm_nis_component": (
                result.composition_fm_innovation**2
                / result.composition_fm_innovation_variance
            ),
            "composition_joint_nis": result.composition_joint_nis,
            "covariance_repaired_predicted": (
                result.covariance_repaired_predicted
            ),
            "covariance_repaired_filtered": (
                result.covariance_repaired_filtered
            ),
            "covariance_repaired_smoothed": (
                smooth.covariance_repaired_smoothed
            ),
        }
    )

    estimates = {
        "predicted": (result.x_predicted, result.p_predicted),
        "filtered": (result.x_filtered, result.p_filtered),
        "smoothed": (smooth.x_smoothed, smooth.p_smoothed),
    }
    for estimate_name, (states, covariances) in estimates.items():
        standard_deviations = _state_standard_deviations(covariances)
        for j, state_name in enumerate(STATE_NAMES):
            frame[f"{state_name}_{estimate_name}"] = states[:, j]
            frame[f"{state_name}_sd_{estimate_name}"] = (
                standard_deviations[:, j]
            )

    frame["predicted_weight_lbs"] = (
        frame["fm_lbs_predicted"]
        + frame["ffm_lbs_predicted"]
        + frame["buffer_lbs_predicted"]
    )
    weight_design = np.array([1.0, 1.0, 1.0, 0.0])
    frame["predicted_weight_sd_lbs"] = np.sqrt(
        np.einsum(
            "i,nij,j->n",
            weight_design,
            result.p_predicted,
            weight_design,
        )
    )
    frame["date"] = frame["date"].dt.strftime("%Y-%m-%d")
    return frame


def build_covariance_output(inputs, result, smooth):
    frames = []
    estimates = {
        "predicted": result.p_predicted,
        "filtered": result.p_filtered,
        "smoothed": smooth.p_smoothed,
    }
    for estimate_name, covariances in estimates.items():
        data = {
            "date": pd.to_datetime(inputs.dates).strftime("%Y-%m-%d"),
            "estimate": estimate_name,
        }
        for i, j in UPPER_TRIANGLE:
            data[f"cov_{STATE_NAMES[i]}__{STATE_NAMES[j]}"] = (
                covariances[:, i, j]
            )
        frames.append(pd.DataFrame(data))
    return pd.concat(frames, ignore_index=True)


def _innovation_metrics(innovation, variance):
    valid = np.isfinite(innovation) & np.isfinite(variance) & (variance > 0)
    residual = innovation[valid]
    predictive_variance = variance[valid]
    if not len(residual):
        return {
            "n": 0,
            "mae": np.nan,
            "rmse": np.nan,
            "mean_nlpd": np.nan,
            "standardized_mean": np.nan,
            "standardized_rms": np.nan,
            "standardized_acf1": np.nan,
            "coverage_80": np.nan,
            "coverage_95": np.nan,
        }

    standardized = residual / np.sqrt(predictive_variance)
    if len(standardized) >= 2:
        standardized_acf1 = np.corrcoef(
            standardized[:-1], standardized[1:]
        )[0, 1]
    else:
        standardized_acf1 = np.nan
    return {
        "n": int(len(residual)),
        "mae": float(np.mean(np.abs(residual))),
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "mean_nlpd": float(
            np.mean(
                0.5
                * (
                    np.log(2.0 * np.pi * predictive_variance)
                    + residual**2 / predictive_variance
                )
            )
        ),
        "standardized_mean": float(np.mean(standardized)),
        "standardized_rms": float(np.sqrt(np.mean(standardized**2))),
        "standardized_acf1": float(standardized_acf1),
        "coverage_80": float(np.mean(np.abs(standardized) <= 1.2815515655)),
        "coverage_95": float(np.mean(np.abs(standardized) <= 1.9599639845)),
    }


def build_temporal_validation(inputs, result):
    rows = []
    segments = ("development", "calibration", "locked_test", "all")
    channels = {
        "total_weight": (
            result.weight_innovation,
            result.weight_innovation_variance,
        ),
        "composition_fm": (
            result.composition_fm_innovation,
            result.composition_fm_innovation_variance,
        ),
    }
    for segment in segments:
        segment_mask = (
            np.ones(len(inputs.dates), dtype=bool)
            if segment == "all"
            else inputs.segment == segment
        )
        for channel, (innovation, variance) in channels.items():
            metrics = _innovation_metrics(
                np.where(segment_mask, innovation, np.nan),
                np.where(segment_mask, variance, np.nan),
            )
            rows.append(
                {
                    "segment": segment,
                    "channel": channel,
                    "scored_pre_update": True,
                    **metrics,
                }
            )
    return pd.DataFrame(rows)


def build_config_output(config):
    row = {
        "config_id": "fixed_primary_v1",
        "selected_by_sweep": False,
        "locked_test_used_for_selection": False,
        **config.__dict__,
        "composition_error_sd_lbs_json": json.dumps(
            COMPOSITION_ERROR_SD_LBS, sort_keys=True
        ),
        "unknown_device_error_sd_lbs_json": json.dumps(
            UNKNOWN_DEVICE_ERROR_SD_LBS, sort_keys=True
        ),
        "input_paths_json": json.dumps(
            {
                key: str(path.relative_to(ROOT))
                for key, path in RAW_INPUT_PATHS.items()
            },
            sort_keys=True,
        ),
    }
    return pd.DataFrame([row])


def _round_numeric(frame, decimals):
    result = frame.copy()
    numeric = result.select_dtypes(include=[np.number]).columns
    result[numeric] = result[numeric].round(decimals)
    return result


def save_outputs(inputs, result, smooth, config):
    state_frame = build_state_output(inputs, result, smooth, config)
    covariance_frame = build_covariance_output(inputs, result, smooth)
    validation_frame = build_temporal_validation(inputs, result)
    config_frame = build_config_output(config)

    _round_numeric(state_frame, 8).to_csv(STATE_OUTPUT, index=False)
    _round_numeric(covariance_frame, 10).to_csv(
        COVARIANCE_OUTPUT, index=False
    )
    _round_numeric(validation_frame, 8).to_csv(
        VALIDATION_OUTPUT, index=False
    )
    config_frame.to_csv(CONFIG_OUTPUT, index=False)
    return state_frame, covariance_frame, validation_frame, config_frame


def main():
    config = ModelConfig()
    print("Loading raw intake, scale-weight, and composition inputs...")
    inputs = load_raw_inputs()
    print(
        f"  {len(inputs.dates):,} days; "
        f"{np.isfinite(inputs.scale_weight_lbs).sum():,} scale weights; "
        f"{np.isfinite(inputs.composition_fm_lbs).sum():,} composition scans"
    )

    print("Running fixed-config forward filter...")
    result = run_forward_filter(inputs, config)
    print("Running separate retrospective RTS smoother...")
    smooth = run_rts_smoother(result)
    _, covariance, validation, _ = save_outputs(
        inputs, result, smooth, config
    )

    initialized = np.isfinite(result.x_filtered[:, 0])
    filtered_smoothed_difference = np.nanmax(
        np.abs(result.x_filtered - smooth.x_smoothed)
    )
    print(
        f"  initialized {initialized.sum():,} days from "
        f"{pd.Timestamp(inputs.dates[result.first_state_index]).date()}"
    )
    print(
        "  covariance repair flags: "
        f"predicted={result.covariance_repaired_predicted.sum()}, "
        f"filtered={result.covariance_repaired_filtered.sum()}, "
        f"smoothed={smooth.covariance_repaired_smoothed.sum()}"
    )
    print(
        f"  max filtered-vs-smoothed state difference: "
        f"{filtered_smoothed_difference:.3f}"
    )
    print(f"Wrote {STATE_OUTPUT}")
    print(f"Wrote {COVARIANCE_OUTPUT} ({len(covariance):,} rows)")
    print(f"Wrote {VALIDATION_OUTPUT}")
    print(f"Wrote {CONFIG_OUTPUT}")

    print("\nPre-update temporal scores:")
    print(
        validation[
            [
                "segment",
                "channel",
                "n",
                "mae",
                "mean_nlpd",
                "standardized_mean",
                "standardized_rms",
                "coverage_80",
                "coverage_95",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
