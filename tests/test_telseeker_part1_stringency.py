import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import TelSeekerPart1


def require_part1_api(name):
    if not hasattr(TelSeekerPart1, name):
        raise AssertionError(f"TelSeekerPart1.{name} is required")
    return getattr(TelSeekerPart1, name)


class TeloReadWindowCountTests(unittest.TestCase):
    def test_window_count_detects_spaced_telomere_repeats_without_double_marker(self):
        classify = require_part1_api("classify_read_by_window_counts")

        sequence = ("TTAGGGAA" * 6) + ("C" * 200)
        self.assertNotIn("TTAGGGTTAGGG", sequence)

        result = classify(
            sequence,
            motif="TTAGGG",
            window_size=1000,
            min_window_repeats=5,
            min_window_density=5.0,
        )

        self.assertTrue(result["passed"])
        self.assertEqual(result["side"], "right")
        self.assertEqual(result["forward_count"], 6)
        self.assertEqual(result["reverse_count"], 0)

    def test_window_count_keeps_reverse_complement_direction(self):
        classify = require_part1_api("classify_read_by_window_counts")

        sequence = ("CCCTAAAA" * 6) + ("G" * 200)

        result = classify(
            sequence,
            motif="TTAGGG",
            window_size=1000,
            min_window_repeats=5,
            min_window_density=5.0,
        )

        self.assertTrue(result["passed"])
        self.assertEqual(result["side"], "left")
        self.assertEqual(result["forward_count"], 0)
        self.assertEqual(result["reverse_count"], 6)

    def test_window_count_rejects_low_density_noise(self):
        classify = require_part1_api("classify_read_by_window_counts")

        result = classify(
            "A" * 500 + "TTAGGG" + "C" * 500,
            motif="TTAGGG",
            window_size=1000,
            min_window_repeats=5,
            min_window_density=5.0,
        )

        self.assertFalse(result["passed"])
        self.assertIsNone(result["side"])


if __name__ == "__main__":
    unittest.main()
