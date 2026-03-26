import importlib
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "analysis"))
tirzepatide_dynamics = importlib.import_module("G_tirzepatide_dynamics")


class TirzepatideDynamicsTests(unittest.TestCase):
    def test_classify_state(self):
        self.assertEqual(tirzepatide_dynamics.classify_state(1700), "restriction")
        self.assertEqual(tirzepatide_dynamics.classify_state(2000), "typical")
        self.assertEqual(tirzepatide_dynamics.classify_state(2500), "high")
        self.assertEqual(tirzepatide_dynamics.classify_state(3000), "binge")

    def test_detect_restriction_runs(self):
        daily = pd.DataFrame(
            {
                "date": pd.date_range("2025-01-01", periods=12, freq="D"),
                "calories": [1700, 1600, 1500, 2200, 2300, 2400, 2500, 1800, 1700, 1600, 2200, 2300],
                "on_tirz": [False] * 12,
                "effective_level": [0.0] * 12,
                "dose_mg": [None] * 12,
                "days_since_injection": [None] * 12,
            }
        )
        runs = tirzepatide_dynamics.detect_restriction_runs(daily)
        self.assertEqual(len(runs), 1)
        self.assertEqual(int(runs.iloc[0]["run_days"]), 3)
        self.assertEqual(int(runs.iloc[0]["next7_binge"]), 0)


if __name__ == "__main__":
    unittest.main()
