"""Deterministic semantic and isolation tests for P4b.

Run:
    uv run --with numpy --with pandas python analysis/test_P4b_forward_state_space.py
"""

from dataclasses import replace
import hashlib
from pathlib import Path
import unittest

import numpy as np

import P4b_forward_state_space as p4b


class P4bTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = p4b.ModelConfig()
        cls.inputs = p4b.load_raw_inputs()
        cls.result = p4b.run_forward_filter(cls.inputs, cls.config)
        cls.smooth = p4b.run_rts_smoother(cls.result)

    def test_raw_input_isolation(self):
        relative = {
            key: str(path.relative_to(p4b.ROOT))
            for key, path in p4b.RAW_INPUT_PATHS.items()
        }
        self.assertEqual(
            relative,
            {
                "intake": "intake/intake_daily.csv",
                "weight": "weight/weight.csv",
                "composition": "composition/composition.csv",
            },
        )
        source = Path(p4b.__file__).read_text()
        forbidden = (
            '"analysis/P1',
            '"analysis/P2',
            '"analysis/P3',
            '"analysis/P4_',
            '"RMR/',
            '"temperature/',
        )
        for token in forbidden:
            self.assertNotIn(token, source)

    def test_primary_multiplier_is_one(self):
        self.assertEqual(self.config.intake_multiplier, 1.0)

    def test_energy_indexing(self):
        state = np.array([70.0, 140.0, 0.0, 2200.0])
        unchanged = p4b.transition_state(state, 2200.0, self.config)
        np.testing.assert_allclose(unchanged, state, atol=1e-12)
        surplus = p4b.transition_state(state, 2550.0, self.config)
        self.assertGreater(surplus[0] + surplus[1], state[0] + state[1])

    def test_no_future_leakage_in_filtered_states(self):
        cutoff = np.datetime64("2024-01-01")
        future = self.inputs.dates > cutoff
        calories = self.inputs.calories_logged.copy()
        scale = self.inputs.scale_weight_lbs.copy()
        composition_weight = self.inputs.composition_weight_lbs.copy()
        composition_fm = self.inputs.composition_fm_lbs.copy()
        calories[future] += 777.0
        scale[future & np.isfinite(scale)] += 11.0
        composition_weight[future & np.isfinite(composition_weight)] += 9.0
        composition_fm[future & np.isfinite(composition_fm)] += 7.0
        perturbed = replace(
            self.inputs,
            calories_logged=calories,
            scale_weight_lbs=scale,
            composition_weight_lbs=composition_weight,
            composition_fm_lbs=composition_fm,
        )
        rerun = p4b.run_forward_filter(perturbed, self.config)
        through_cutoff = self.inputs.dates <= cutoff
        np.testing.assert_allclose(
            self.result.x_filtered[through_cutoff],
            rerun.x_filtered[through_cutoff],
            atol=1e-12,
            rtol=0.0,
            equal_nan=True,
        )
        np.testing.assert_allclose(
            self.result.p_filtered[through_cutoff],
            rerun.p_filtered[through_cutoff],
            atol=1e-12,
            rtol=0.0,
            equal_nan=True,
        )

    def test_repeated_runs_are_deterministic(self):
        rerun = p4b.run_forward_filter(self.inputs, self.config)
        first = hashlib.sha256(
            np.nan_to_num(self.result.x_filtered, nan=-9999.0).tobytes()
            + np.nan_to_num(self.result.p_filtered, nan=-9999.0).tobytes()
        ).hexdigest()
        second = hashlib.sha256(
            np.nan_to_num(rerun.x_filtered, nan=-9999.0).tobytes()
            + np.nan_to_num(rerun.p_filtered, nan=-9999.0).tobytes()
        ).hexdigest()
        self.assertEqual(first, second)

    def test_covariances_are_symmetric_psd(self):
        covariance_sets = (
            self.result.p_predicted,
            self.result.p_filtered,
            self.smooth.p_smoothed,
        )
        for covariances in covariance_sets:
            valid = np.isfinite(covariances[:, 0, 0])
            for covariance in covariances[valid]:
                np.testing.assert_allclose(
                    covariance, covariance.T, atol=1e-9, rtol=0.0
                )
                self.assertGreaterEqual(
                    np.linalg.eigvalsh(covariance).min(), -1e-8
                )

    def test_filtered_and_smoothed_semantics(self):
        first = self.result.first_state_index
        last = np.where(np.isfinite(self.result.x_filtered[:, 0]))[0][-1]
        self.assertGreater(
            np.nanmax(
                np.abs(
                    self.result.x_filtered[first:last]
                    - self.smooth.x_smoothed[first:last]
                )
            ),
            1e-6,
        )
        np.testing.assert_allclose(
            self.result.x_filtered[last],
            self.smooth.x_smoothed[last],
            atol=1e-12,
            rtol=0.0,
        )
        np.testing.assert_allclose(
            self.result.p_filtered[last],
            self.smooth.p_smoothed[last],
            atol=1e-12,
            rtol=0.0,
        )

    def test_no_duplicate_total_weight_or_lm_assimilation(self):
        both = np.isfinite(self.inputs.scale_weight_lbs) & np.isfinite(
            self.inputs.composition_fm_lbs
        )
        self.assertTrue(both.any())
        self.assertTrue(self.result.scale_weight_used[both].all())
        self.assertFalse(self.result.composition_weight_used[both].any())
        self.assertTrue(self.result.composition_fm_used[both].all())
        self.assertNotIn("lm", " ".join(self.result.update_source.astype(str)))

    def test_every_known_device_has_error_mapping(self):
        devices = set(self.inputs.composition_device)
        devices.discard("")
        self.assertEqual(
            devices - set(p4b.COMPOSITION_ERROR_SD_LBS),
            set(),
        )
        self.assertFalse(self.result.device_error_fallback.any())

    def test_exported_standard_deviations_are_complete(self):
        frame = p4b.build_state_output(
            self.inputs, self.result, self.smooth, self.config
        )
        initialized = np.arange(len(frame)) >= self.result.first_state_index
        sd_columns = [
            column
            for column in frame.columns
            if "_sd_predicted" in column
            or "_sd_filtered" in column
            or "_sd_smoothed" in column
        ]
        self.assertEqual(len(sd_columns), 12)
        self.assertFalse(frame.loc[initialized, sd_columns].isna().any().any())

    def test_temporal_scores_are_pre_update(self):
        validation = p4b.build_temporal_validation(
            self.inputs, self.result
        )
        self.assertTrue(validation["scored_pre_update"].all())
        locked = validation[validation["segment"] == "locked_test"]
        self.assertGreater(locked["n"].sum(), 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
