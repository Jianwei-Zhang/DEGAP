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


class TeloReadTerminalRatioTests(unittest.TestCase):
    def test_terminal_ratio_detects_forward_telomere_end(self):
        classify = require_part1_api("classify_read_by_terminal_telomere")

        result = classify(
            ("TTAGGG" * 60) + ("C" * 300),
            motif="TTAGGG",
            tel_n=100,
            tel_r=0.6,
            tel_mm=0,
        )

        self.assertTrue(result["passed"])
        self.assertEqual(result["side"], "right")
        self.assertEqual(result["forward_hits"], 60)
        self.assertEqual(result["reverse_hits"], 0)

    def test_terminal_ratio_detects_reverse_complement_telomere_end(self):
        classify = require_part1_api("classify_read_by_terminal_telomere")

        result = classify(
            ("CCCTAA" * 60) + ("G" * 300),
            motif="TTAGGG",
            tel_n=100,
            tel_r=0.6,
            tel_mm=0,
        )

        self.assertTrue(result["passed"])
        self.assertEqual(result["side"], "left")
        self.assertEqual(result["forward_hits"], 0)
        self.assertEqual(result["reverse_hits"], 60)

    def test_terminal_ratio_rejects_internal_telomere_signal(self):
        classify = require_part1_api("classify_read_by_terminal_telomere")

        result = classify(
            ("A" * 600) + ("TTAGGG" * 80) + ("C" * 600),
            motif="TTAGGG",
            tel_n=100,
            tel_r=0.6,
            tel_mm=0,
        )

        self.assertFalse(result["passed"])
        self.assertIsNone(result["side"])

    def test_terminal_scan_slides_one_base_after_a_missed_window(self):
        rotations = set(TelSeekerPart1.motif_rotations("TTAGGG"))
        count_hits = require_part1_api("count_terminal_telomere_hits")

        self.assertEqual(
            count_hits("TTAGGGATAGGGT", rotations, motif_length=6, mismatch=0),
            2,
        )

    def test_terminal_ratio_can_allow_one_mismatch_per_motif_unit(self):
        classify = require_part1_api("classify_read_by_terminal_telomere")

        strict_result = classify("TTAGGA", "TTAGGG", tel_n=1, tel_r=1.0, tel_mm=0)
        tolerant_result = classify("TTAGGA", "TTAGGG", tel_n=1, tel_r=1.0, tel_mm=1)

        self.assertFalse(strict_result["passed"])
        self.assertTrue(tolerant_result["passed"])

    def test_terminal_ratio_reports_both_direction_passes(self):
        classify = require_part1_api("classify_read_by_terminal_telomere")

        result = classify(
            ("CCCTAA" * 3) + ("A" * 20) + ("TTAGGG" * 3),
            "TTAGGG",
            tel_n=3,
            tel_r=1.0,
            tel_mm=0,
        )

        self.assertTrue(result["left_pass"])
        self.assertTrue(result["right_pass"])


class TeloReadWorkerTerminalRatioTests(unittest.TestCase):
    def test_worker_uses_terminal_ratio_without_double_marker_requirement(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            input_file = tmp_path / "reads.fa"
            input_file.write_text(
                ">spaced_forward\n"
                + ("TTAGGG" * 3)
                + ("C" * 200)
                + "\n"
            )

            extractor = TelSeekerPart1.TeloReadsExtractor(
                output_dir=str(tmp_path),
                motif="TTAGGG",
                tel_n=4,
                tel_r=0.5,
                tel_mm=0,
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


class TeloReadCliTerminalRatioTests(unittest.TestCase):
    def test_part1_help_exposes_terminal_ratio_parameters(self):
        result = subprocess.run(
            [sys.executable, str(BIN_DIR / "TelSeekerPart1.py"), "--help"],
            capture_output=True,
            text=True,
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("--tel-n", result.stdout)
        self.assertIn("--tel-r", result.stdout)
        self.assertIn("--tel-mm", result.stdout)


if __name__ == "__main__":
    unittest.main()
