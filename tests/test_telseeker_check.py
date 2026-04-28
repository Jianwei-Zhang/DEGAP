import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import TelSeekerCheck


def require_checker_api(name):
    if not hasattr(TelSeekerCheck, name):
        raise AssertionError(f"TelSeekerCheck.{name} is required")
    return getattr(TelSeekerCheck, name)


class TelSeekerCheckWindowTests(unittest.TestCase):
    def test_counts_forward_and_reverse_motif_repeats_per_window(self):
        count_windows = require_checker_api("count_telomere_repeat_windows")

        sequence = "TTAGGG" * 4 + "CCCTAA" * 2
        windows = count_windows(sequence, "TTAGGG", window_size=24)

        self.assertEqual(len(windows), 2)
        self.assertEqual(windows[0]["start"], 0)
        self.assertEqual(windows[0]["end"], 24)
        self.assertEqual(windows[0]["forward_count"], 4)
        self.assertEqual(windows[0]["reverse_count"], 0)
        self.assertEqual(windows[0]["total_count"], 4)
        self.assertAlmostEqual(windows[0]["density_per_kb"], 4 / 0.024)
        self.assertEqual(windows[1]["forward_count"], 0)
        self.assertEqual(windows[1]["reverse_count"], 2)
        self.assertEqual(windows[1]["total_count"], 2)

    def test_classification_requires_terminal_enrichment_over_internal_repeats(self):
        count_windows = require_checker_api("count_telomere_repeat_windows")
        classify = require_checker_api("classify_terminal_telomere_signal")

        internal_only = "A" * 60 + "TTAGGG" * 10 + "A" * 60
        internal_windows = count_windows(internal_only, "TTAGGG", window_size=60)

        internal_result = classify(
            internal_windows,
            end="Left",
            terminal_window_count=1,
            min_terminal_repeats=4,
            min_terminal_density=50.0,
            min_terminal_to_internal_ratio=2.0,
        )

        self.assertEqual(internal_result["status"], "untelomeric")
        self.assertEqual(internal_result["terminal_repeat_count"], 0)
        self.assertGreaterEqual(internal_result["internal_max_density"], 50.0)
        self.assertIn("terminal", internal_result["reason"].lower())

        terminal = "TTAGGG" * 10 + "A" * 120
        terminal_windows = count_windows(terminal, "TTAGGG", window_size=60)

        terminal_result = classify(
            terminal_windows,
            end="Left",
            terminal_window_count=1,
            min_terminal_repeats=4,
            min_terminal_density=50.0,
            min_terminal_to_internal_ratio=2.0,
        )

        self.assertEqual(terminal_result["status"], "telomeric")
        self.assertEqual(terminal_result["terminal_repeat_count"], 10)
        self.assertEqual(terminal_result["internal_max_density"], 0.0)
        self.assertEqual(terminal_result["confidence"], "high")


class TelomereCheckerIntegrationTests(unittest.TestCase):
    def test_checker_uses_terminal_window_enrichment_not_internal_repeats(self):
        checker = TelSeekerCheck.TelomereChecker(
            genome_file="unused.fa",
            motif="TTAGGG",
            output_dir="unused",
            check_length=240,
            window_size=60,
            terminal_windows=1,
            min_terminal_repeats=4,
            min_terminal_density=50.0,
            min_terminal_to_internal_ratio=2.0,
        )

        internal_only = "A" * 60 + "TTAGGG" * 10 + "A" * 120
        result = checker._check_end("Chr01", internal_only, "Left")

        self.assertEqual(result["status"], "untelomeric")
        self.assertEqual(result["terminal_repeat_count"], 0)
        self.assertGreaterEqual(result["internal_max_density"], 50.0)
        self.assertEqual(result["method"], "window_count")

    def test_checker_report_keeps_trc_and_adds_review_fields(self):
        with self.subTest("full checker run writes richer review output"):
            import tempfile

            with tempfile.TemporaryDirectory() as tmp:
                tmp_path = Path(tmp)
                genome = tmp_path / "genome.fa"
                genome.write_text(">Chr01\n" + "TTAGGG" * 10 + "A" * 120 + "\n")

                checker = TelSeekerCheck.TelomereChecker(
                    genome_file=str(genome),
                    motif="TTAGGG",
                    output_dir=str(tmp_path / "out"),
                    check_length=180,
                    extract_length=30,
                    window_size=60,
                    terminal_windows=1,
                    min_terminal_repeats=4,
                    min_terminal_density=50.0,
                    min_terminal_to_internal_ratio=2.0,
                )

                checker.run()

                check_dir = tmp_path / "out" / "genome.telomere.check"
                csv_text = (check_dir / "genome.telomere.check.csv").read_text()
                self.assertIn("Chr01,Left,telomeric,TRC=", csv_text)
                self.assertIn("Method=window_count", csv_text)
                self.assertIn("Motif=TTAGGG", csv_text)
                self.assertIn("ForwardCount=10", csv_text)
                self.assertIn("ReverseCount=0", csv_text)
                self.assertIn("TerminalRepeatCount=10", csv_text)
                self.assertIn("InternalMaxDensity=0.00", csv_text)
                self.assertIn("TerminalToInternalRatio=Inf", csv_text)
                self.assertIn("Confidence=high", csv_text)

                need_extension = (
                    check_dir / "need_extension_chr_end.txt"
                ).read_text().splitlines()
                self.assertEqual(need_extension, ["Chr01.R"])


if __name__ == "__main__":
    unittest.main()
