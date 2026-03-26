"""Causality and correctness tests for binge_analysis.

Verifies that set point distance and deficit features are strictly trailing
(no future data leakage), and that the hand-rolled logistic regression / AUC
produce correct results on synthetic data.
"""
import importlib
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "analysis"))
binge_analysis = importlib.import_module("C_binge_analysis")


class BingeFeatureTests(unittest.TestCase):
    def test_set_point_is_trailing_not_centered(self):
        fat_mass = list(range(10, 110))
        fat_mass[-1] = 500
        daily = pd.DataFrame({
            "date": pd.date_range("2024-01-01", periods=100, freq="D"),
            "calories": [2200] * 100,
            "protein_g": [100] * 100,
            "fat_g": [70] * 100,
            "carbs_g": [200] * 100,
            "fat_mass_lbs_filtered": fat_mass,
            "fat_mass_lbs": fat_mass,
            "tdee_filtered": [2200] * 100,
            "tdee": [2200] * 100,
            "effective_level": [0] * 100,
        })
        features = binge_analysis.compute_features(daily)
        # Day 98 should not be affected by the extreme day 100 value.
        self.assertLess(features.loc[97, "dist_90d"], 50)

    def test_deficit_uses_prior_days_only(self):
        daily = pd.DataFrame({
            "date": pd.date_range("2024-01-01", periods=10, freq="D"),
            "calories": [2000] * 9 + [5000],
            "protein_g": [100] * 10,
            "fat_g": [70] * 10,
            "carbs_g": [200] * 10,
            "fat_mass_lbs_filtered": np.linspace(20, 21, 10),
            "fat_mass_lbs": np.linspace(20, 21, 10),
            "tdee_filtered": [2200] * 10,
            "tdee": [2200] * 10,
            "effective_level": [0] * 10,
        })
        features = binge_analysis.compute_features(daily)
        self.assertLess(features.loc[9, "cum_deficit_7d"], 0)


class LogisticRegressionTests(unittest.TestCase):
    def test_auc_perfect_separation(self):
        y = np.array([0, 0, 1, 1], dtype=float)
        scores = np.array([0.1, 0.2, 0.8, 0.9])
        self.assertAlmostEqual(binge_analysis.roc_auc_score_np(y, scores), 1.0)

    def test_logistic_regression_fits_signal(self):
        X = np.array([[-2.0], [-1.0], [1.0], [2.0]])
        y = np.array([0.0, 0.0, 1.0, 1.0])
        beta = binge_analysis.fit_logistic_regression(X, y, l2_penalty=0.1)
        probs = binge_analysis.predict_logistic_regression(X, beta)
        self.assertLess(probs[0], 0.5)
        self.assertGreater(probs[-1], 0.5)


if __name__ == "__main__":
    unittest.main()
