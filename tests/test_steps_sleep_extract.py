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
            "walking",
        )
        self.assertEqual(
            steps_extract.infer_exercise_label("11007"),
            "bike",
        )
        self.assertEqual(
            steps_extract.infer_exercise_label("15003"),
            "indoor_bike",
        )
        self.assertEqual(
            steps_extract.infer_exercise_label("99999"),
            "unknown",
        )

    def test_exercise_type_map_includes_other(self):
        self.assertEqual(
            steps_extract.infer_exercise_label("0"),
            "other",
        )


if __name__ == "__main__":
    unittest.main()
