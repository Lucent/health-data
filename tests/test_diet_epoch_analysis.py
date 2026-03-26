import importlib
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "analysis"))
diet_epoch_analysis = importlib.import_module("O_diet_epoch_analysis")


class DietEpochAnalysisTests(unittest.TestCase):
    def test_normalize_family(self):
        self.assertEqual(diet_epoch_analysis.normalize_family("potato_diet_4"), "potato_diet")
        self.assertEqual(diet_epoch_analysis.normalize_family("keto_phase_2"), "keto_phase")
        self.assertEqual(diet_epoch_analysis.normalize_family("weekend_fasting"), "weekend_fasting")

    def test_potato_window_summary_has_three_phases_per_epoch(self):
        daily = diet_epoch_analysis.load_daily()
        epochs = diet_epoch_analysis.load_epochs()
        potato = diet_epoch_analysis.potato_window_summary(daily, epochs)

        counts = potato.groupby("label")["phase"].nunique()
        self.assertTrue((counts == 3).all())
        self.assertIn("epoch", set(potato["phase"]))
        self.assertIn("pre", set(potato["phase"]))
        self.assertIn("post", set(potato["phase"]))

    def test_potato_epoch_contrast_groups(self):
        daily = diet_epoch_analysis.load_daily()
        epochs = diet_epoch_analysis.load_epochs()
        potato = diet_epoch_analysis.potato_window_summary(daily, epochs)
        contrast = diet_epoch_analysis.potato_contrast(daily, potato)

        groups = set(contrast["group"])
        self.assertIn("potato_epoch_days", groups)
        self.assertIn("matched_non_potato_days", groups)
        self.assertIn("potato_pre_window_mean", groups)
        self.assertIn("potato_epoch_window_mean", groups)
        self.assertIn("potato_post_window_mean", groups)


if __name__ == "__main__":
    unittest.main()
