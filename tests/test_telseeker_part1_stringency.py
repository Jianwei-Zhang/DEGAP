import sys
import subprocess
import tempfile
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


class TeloReadWorkerStringencyTests(unittest.TestCase):
    def test_normal_worker_uses_window_count_without_double_marker_requirement(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            input_file = tmp_path / "reads.fa"
            input_file.write_text(
                ">spaced_forward\n"
                + ("TTAGGGAA" * 6)
                + ("C" * 200)
                + "\n"
            )

            extractor = TelSeekerPart1.TeloReadsExtractor(
                output_dir=str(tmp_path),
                motif="TTAGGG",
                telo_read_stringency="normal",
            )
            script_path = extractor._create_processing_script(str(tmp_path))

            result = subprocess.run(
                [sys.executable, script_path, str(input_file)],
                capture_output=True,
                text=True,
            )

            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            left_output = tmp_path / f"left_{input_file.name}"
            self.assertNotIn(">spaced_forward", left_output.read_text())
            right_output = tmp_path / f"right_{input_file.name}"
            self.assertIn(">spaced_forward", right_output.read_text())


class TeloReadCliStringencyTests(unittest.TestCase):
    def test_part1_help_exposes_single_stringency_parameter(self):
        result = subprocess.run(
            [sys.executable, str(BIN_DIR / "TelSeekerPart1.py"), "--help"],
            capture_output=True,
            text=True,
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("--telo-read-stringency", result.stdout)
        self.assertIn("strict", result.stdout)
        self.assertIn("relaxed", result.stdout)


if __name__ == "__main__":
    unittest.main()
