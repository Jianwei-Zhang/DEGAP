import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from TelSeekerTargets import parse_target_ends


class TelSeekerTargetsTests(unittest.TestCase):
    def write_genome(self, tmp_path: Path) -> Path:
        genome = tmp_path / "genome.fa"
        genome.write_text(">Chr01\nACGT\n>Chr02\nACGT\n")
        return genome

    def test_parse_multiple_cli_target_ends(self):
        with tempfile.TemporaryDirectory() as tmp:
            genome = self.write_genome(Path(tmp))

            result = parse_target_ends(["Chr01.L", "Chr01.R", "Chr02.L"], genome)

            self.assertEqual(result, ["Chr01.L", "Chr01.R", "Chr02.L"])

    def test_parse_target_ends_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            genome = self.write_genome(tmp_path)
            targets = tmp_path / "targets.txt"
            targets.write_text("# target ends\nChr01.L\n\nChr02.R\n")

            result = parse_target_ends([str(targets)], genome)

            self.assertEqual(result, ["Chr01.L", "Chr02.R"])

    def test_parse_target_ends_deduplicates_in_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            genome = self.write_genome(Path(tmp))

            result = parse_target_ends(["Chr01.L", "Chr01.L", "Chr02.R"], genome)

            self.assertEqual(result, ["Chr01.L", "Chr02.R"])

    def test_empty_target_file_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            genome = self.write_genome(tmp_path)
            targets = tmp_path / "targets.txt"
            targets.write_text("\n# none\n")

            with self.assertRaisesRegex(ValueError, "No target ends"):
                parse_target_ends([str(targets)], genome)

    def test_invalid_target_format_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            genome = self.write_genome(Path(tmp))

            with self.assertRaisesRegex(ValueError, "Invalid target end"):
                parse_target_ends(["Chr01.left"], genome)

    def test_unknown_chromosome_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            genome = self.write_genome(Path(tmp))

            with self.assertRaisesRegex(ValueError, "not found in genome"):
                parse_target_ends(["Chr99.L"], genome)


if __name__ == "__main__":
    unittest.main()
