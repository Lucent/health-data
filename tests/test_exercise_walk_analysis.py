import importlib.util
from pathlib import Path
import unittest

import pandas as pd


MODULE_PATH = Path(__file__).resolve().parent.parent / "analysis" / "V_exercise_walk_analysis.py"
SPEC = importlib.util.spec_from_file_location("V_exercise_walk_analysis", MODULE_PATH)
exercise_walk_analysis = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(exercise_walk_analysis)


class ExerciseWalkAnalysisTests(unittest.TestCase):
    def test_type_labels(self):
        self.assertEqual(
            exercise_walk_analysis.TYPE_LABELS[1001],
            ("walking", "high"),
        )
        self.assertEqual(
            exercise_walk_analysis.TYPE_LABELS[1002],
            ("running", "high"),
        )

    def test_classify_daylight_walk_regimes(self):
        exercises = pd.DataFrame(
            [
                {
                    "exercise_type": 1001,
                    "start_time": pd.Timestamp("2025-06-10 14:00:00"),
                    "end_time": pd.Timestamp("2025-06-10 14:32:00"),
                    "duration_min": 32.0,
                    "count": 3200,
                    "distance": 2500,
                    "source_type": 8.0,
                },
                {
                    "exercise_type": 1001,
                    "start_time": pd.Timestamp("2025-06-10 14:50:00"),
                    "end_time": pd.Timestamp("2025-06-10 15:20:00"),
                    "duration_min": 30.0,
                    "count": 3000,
                    "distance": 2300,
                    "source_type": 8.0,
                },
                {
                    "exercise_type": 1001,
                    "start_time": pd.Timestamp("2025-06-11 10:00:00"),
                    "end_time": pd.Timestamp("2025-06-11 10:20:00"),
                    "duration_min": 20.0,
                    "count": 2000,
                    "distance": 1500,
                    "source_type": 8.0,
                },
            ]
        )
        out = exercise_walk_analysis.classify_daylight_walk_regimes(exercises)
        self.assertEqual(len(out), 1)
        self.assertEqual(out.iloc[0]["date"], "2025-06-10")
        self.assertEqual(out.iloc[0]["walk_regime"], "paired_daylight_walk")
        self.assertAlmostEqual(out.iloc[0]["regime_steps"], 6200.0)


if __name__ == "__main__":
    unittest.main()
