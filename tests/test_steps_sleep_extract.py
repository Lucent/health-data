import importlib.util
from pathlib import Path
import unittest


MODULE_PATH = Path(__file__).resolve().parent.parent / "steps-sleep" / "extract.py"
SPEC = importlib.util.spec_from_file_location("steps_sleep_extract", MODULE_PATH)
steps_extract = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(steps_extract)


class StepsSleepExtractTests(unittest.TestCase):
    def test_parse_offset(self):
        td = steps_extract.parse_offset("UTC-0500")
        self.assertEqual(td.total_seconds(), -5 * 3600)

    def test_infer_exercise_label(self):
        self.assertEqual(
            steps_extract.infer_exercise_label("1001"),
            ("walking", "high"),
        )
        self.assertEqual(
            steps_extract.infer_exercise_label("11007"),
            ("bike", "high"),
        )
        self.assertEqual(
            steps_extract.infer_exercise_label("15003"),
            ("indoor_bike", "high"),
        )
        self.assertEqual(
            steps_extract.infer_exercise_label("99999"),
            ("unknown", "low"),
        )

    def test_summarize_exercises_daily(self):
        exercises = [
            {
                "date": "2024-01-01",
                "duration_min": 30,
                "count": 3000,
                "distance": 2000,
                "calorie": 150,
            },
            {
                "date": "2024-01-01",
                "duration_min": 20,
                "count": 2000,
                "distance": 1200,
                "calorie": 100,
            },
        ]
        steps = [{"date": "2024-01-01", "steps": "8000"}]
        out = steps_extract.summarize_exercises_daily(exercises, steps)
        self.assertEqual(len(out), 1)
        self.assertEqual(out[0]["exercise_sessions"], 2)
        self.assertAlmostEqual(out[0]["exercise_step_fraction"], 0.625)


if __name__ == "__main__":
    unittest.main()
