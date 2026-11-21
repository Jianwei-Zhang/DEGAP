from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from ExtensionFamilies import AttemptResult
from FamilyTelomere import run_telomere
from TelSeekerPart2 import TelomereExtender, ConnectionChecker


class TelomereFamilyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.seed = self.root / "seed.fa"
        self.seed.write_text(">seed\nAAAA\n")
        (self.root / "telo.reads.fa").write_text(">telo\nCCCC\n")
        self.context = patch("FamilyTelomere.prepare_context", return_value={"fixture": 1})
        self.context.start()
        self.addCleanup(self.context.stop)

    def extender(self, end="R", rounds=3, length=20):
        item = TelomereExtender.__new__(TelomereExtender)
        item.chr_end = "chr." + end
        item.end_type = end
        item.work_path = self.root
        item.extension_dir = self.root / "extension"
        item.seed_file = self.seed
        item.telo_reads_file = self.root / "telo.reads.fa"
        item.base = SimpleNamespace(out=str(item.extension_dir), terminalSeq=str(self.seed))
        item.max_rounds, item.max_extension_length = rounds, length
        item.edge, item.kmer_num, item.kmer_size = 500, 10, 41
        item._check_connection = lambda work, sequence: None
        return item

    def growth(self, family, directory):
        path = directory / "extended.fa"
        current = Path(family.current.path).read_text().splitlines()[1]
        path.write_text(">extended\n" + current + "CC\n")
        return AttemptResult(str(path), ("read-" + str(family.attempts),))

    def connect(self, work, sequence):
        path = work / "linker.fa"
        path.write_text(Path(sequence).read_text())
        return {"success": True, "read_id": "telo", "linker_file": path,
                "alignment": {"identity": 99.0, "len1": 4}}

    def test_single_family_no_seed_self_closure_and_correct_direction(self):
        for end, direction in (("L", "leftward"), ("R", "rightward")):
            with self.subTest(end=end):
                item = self.extender(end, rounds=1)
                item.extension_dir = self.root / end
                seen = []
                def extend(family, directory):
                    seen.append(family.direction)
                    return self.growth(family, directory)
                result = run_telomere(item, extend=extend)
                self.assertEqual(seen, [direction])
                self.assertFalse(result["success"])
                self.assertEqual(result["rounds"], 1)
                self.assertEqual(result["total_extension_length"], 2)
                self.assertEqual(result["stop_reason"], "max_rounds_reached")

    def test_terminal_round_growth_is_counted_and_can_cross_length_limit(self):
        item = self.extender(length=1)
        item._check_connection = self.connect
        result = run_telomere(item, extend=self.growth)
        self.assertTrue(result["success"])
        self.assertEqual(result["rounds"], 1)
        self.assertEqual(result["total_extension_length"], 2)
        self.assertEqual(result["stop_reason"], "telomere_connected")
        self.assertEqual(Path(result["final_seq_file"]).read_text(), ">extended\nAAAACC\n")

    def test_failed_attempt_counted_without_growth_or_check(self):
        item = self.extender()
        item._check_connection = lambda *a: self.fail("No sequence to check")
        result = run_telomere(item, extend=lambda *a: AttemptResult(stop_reason="no_reads"))
        self.assertFalse(result["success"])
        self.assertEqual((result["rounds"], result["total_extension_length"]), (1, 0))
        self.assertEqual(result["stop_reason"], "no_reads")

    def test_zero_budget_does_not_attempt_or_self_close(self):
        item = self.extender(rounds=0)
        result = run_telomere(item, extend=lambda *a: self.fail("Budget exhausted"))
        self.assertEqual(result["rounds"], 0)
        self.assertFalse(result["success"])

    def test_resume_reuses_terminal_receipt_and_rejects_corrupt_linker(self):
        item = self.extender()
        item._check_connection = self.connect
        result = run_telomere(item, extend=self.growth)
        again = self.extender()
        again._check_connection = lambda *a: self.fail("Detection already completed")
        restored = run_telomere(again, extend=lambda *a: self.fail("Assembly already completed"))
        self.assertEqual(result, restored)
        Path(result["connection_info"]["linker_file"]).write_text("broken")
        with self.assertRaisesRegex(ValueError, "artifact changed"):
            run_telomere(again, extend=lambda *a: self.fail("Cannot recompute corrupt evidence"))

    def test_read_history_is_retained_between_attempts(self):
        item = self.extender(rounds=2)
        seen = []
        def extend(family, directory):
            seen.append((list(family.used_reads), list(family.last_reads)))
            return self.growth(family, directory)
        run_telomere(item, extend=extend)
        self.assertEqual(seen, [([], []), (["read-1"], ["read-1"])])

    def test_tool_error_remains_pending_not_biological_stop(self):
        item = self.extender()
        def fail(*args):
            raise RuntimeError("tool failed")
        with self.assertRaisesRegex(RuntimeError, "tool failed"):
            run_telomere(item, extend=fail)
        resumed = self.extender()
        resumed._check_connection = self.connect
        result = run_telomere(resumed, extend=self.growth)
        self.assertEqual(result["rounds"], 1)
        self.assertTrue(result["success"])

    def test_report_reads_family_manifest_not_legacy_round_folders(self):
        from TelSeekerVisualizer import TelSeekerPart2LogParser
        item = self.extender(rounds=2)
        item._check_connection = self.connect
        run_telomere(item, extend=self.growth)
        data = TelSeekerPart2LogParser(self.root).parse()
        self.assertTrue(data["final_info"]["success"])
        self.assertEqual(data["final_info"]["method"], "extension")
        self.assertEqual(data["rounds"][0]["extension_length"], 2)
        self.assertTrue(data["rounds"][0]["connection_found"])
        self.assertEqual(data["config"]["max_rounds"], 2)

    def test_old_directory_rejected_unchanged(self):
        item = self.extender()
        item.extension_dir.mkdir()
        log = item.extension_dir / "extension.log"
        log.write_text("legacy")
        with self.assertRaisesRegex(ValueError, "Legacy"):
            run_telomere(item, extend=self.growth)
        self.assertEqual(list(item.extension_dir.iterdir()), [log])
        self.assertEqual(log.read_text(), "legacy")

    def test_changed_telomeric_input_rejects_resume(self):
        item = self.extender(rounds=0)
        run_telomere(item)
        item.telo_reads_file.write_text(">changed\nGGGG\n")
        with self.assertRaises(ValueError):
            run_telomere(self.extender(rounds=0))


class TelomereAlignmentTests(unittest.TestCase):
    def test_coordinates_trim_overhangs_without_right_off_by_one(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ref, qry, out = [root / name for name in ("ref.fa", "qry.fa", "out.fa")]
            checker = ConnectionChecker()
            ref.write_text(">ref\nAAAACCCC\n")
            qry.write_text(">qry\nTCCCCGGG\n")
            aln = {"qry_id": "qry", "s1": 5, "e1": 8, "s2": 2, "e2": 5,
                   "len1": 4, "identity": 100}
            self.assertTrue(checker.create_linked_sequence(ref, qry, aln, "R", out))
            self.assertEqual(out.read_text().splitlines()[1], "AAAACCCCGGG")
            ref.write_text(">ref\nCCCCGGGG\n")
            qry.write_text(">qry\nAAAACCCCT\n")
            aln.update(s1=1, e1=4, s2=5, e2=8)
            self.assertTrue(checker.create_linked_sequence(ref, qry, aln, "L", out))
            self.assertEqual(out.read_text().splitlines()[1], "AAAACCCCGGGG")

    def test_checked_missing_input_is_not_no_connection(self):
        checker = ConnectionChecker(strict_tools=True)
        with self.assertRaises(FileNotFoundError):
            checker.run_mummer_alignment(Path("/nonexistent/tel-ref"), Path("/missing/tel-qry"), Path("unused"))


if __name__ == "__main__":
    unittest.main()
