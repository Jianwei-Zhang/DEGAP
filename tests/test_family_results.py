import json
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))
from ExtensionFamilies import AttemptResult, Connection, FamilySearch
from FamilyResult import load_result, publish_search
from AutoGapfillerPart3 import GapFillingAnalyzer, JobExecutor, JobManager
from AutoGapfillerVisualizer import GapFillingVisualizer


def fasta_records(path):
    records = {}
    key = None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            key = line[1:].split()[0]
            records[key] = ""
        elif line:
            records[key] += line
    return records


class FamilyResultTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.left, self.right = self.root / "left.fa", self.root / "right.fa"
        self.left.write_text(">left\nAAAACCCC\n")
        self.right.write_text(">right\nCCCCGGGG\n")
        self.job = self.root / "03.gap_filling" / "chr1.gap1.bidirectional"

    def search(self, closed=True):
        search = FamilySearch(self.job, [("L", "rightward", self.left), ("R", "leftward", self.right)],
                              primary="L", context={"fixture": 1}, max_rounds=0)
        def compare(a, b):
            return [Connection("L", 0, "R", 0, 100, 4, 0, 4, 8, 0, 4)] if closed else []
        search.run(lambda f, d: self.fail("zero budget cannot extend"), compare)
        return search

    def assert_agp_reconstructs_fasta(self, directory, result):
        components = fasta_records(directory / "family-components.fa")
        reconstructed = {}
        for row in (directory / "process.agp").read_text().splitlines():
            obj, start, end, part, kind, component, first, last, orientation = row.split("\t")
            self.assertEqual((kind, orientation), ("W", "+"))
            sequence = components[component][int(first)-1:int(last)]
            self.assertEqual(len(sequence), int(end)-int(start)+1)
            existing = reconstructed.get(obj, "")
            self.assertEqual(len(existing)+1, int(start))
            reconstructed[obj] = existing + sequence
        self.assertEqual(reconstructed, fasta_records(directory / result["final_sequence"]))

    def test_closed_publication_and_agp_reconstruct_exact_sequence(self):
        result = publish_search(self.search())
        self.assertEqual(list(fasta_records(self.job / result["final_sequence"]).values()), ["AAAACCCCGGGG"])
        self.assertEqual((result["status"], result["gap_length"]), ("closed", -4))
        self.assert_agp_reconstructs_fasta(self.job, result)
        self.assertEqual(load_result(self.job), result)
        self.assertIn("L/0", (self.job / "family-report.html").read_text())

    def test_unclosed_publication_keeps_two_separate_sequences(self):
        result = publish_search(self.search(False))
        self.assertEqual(len(fasta_records(self.job / result["final_sequence"])), 2)
        self.assertIsNone(result["gap_length"])
        self.assertIsNone(result["winner"])
        self.assertNotIn("GAP can be closed!", (self.job / "process.log").read_text())
        self.assert_agp_reconstructs_fasta(self.job, result)

    def test_result_resume_reuses_verified_outputs(self):
        search = self.search()
        result = publish_search(search)
        path = self.job / result["final_sequence"]
        stamp = path.stat().st_mtime_ns
        self.assertEqual(publish_search(search), result)
        self.assertEqual(path.stat().st_mtime_ns, stamp)

    def test_corruption_is_not_success_even_if_log_claims_closed(self):
        search = self.search()
        result = publish_search(search)
        (self.job / result["final_sequence"]).write_text(">wrong\nNNNN\n")
        with self.assertRaisesRegex(ValueError, "artifact changed"):
            load_result(self.job)
        with self.assertRaises(ValueError):
            publish_search(search)
        executor = JobExecutor.__new__(JobExecutor)
        self.assertFalse(executor._check_job_success(self.job))

    def test_auto_analyzer_uses_manifest_not_single_side_log_rounds(self):
        result = publish_search(self.search())
        gap = {"chromosome": "chr1", "gap_id": "gap1", "left": None, "right": None,
               "bidirectional": {"job_dir": self.job, "job_name": self.job.name},
               "gap_info": {"start": 9, "end": 12, "left_len": 8, "right_len": 8, "gap_len": 4}}
        analyzer = GapFillingAnalyzer(self.root)
        analyzed = analyzer._analyze_single_gap(gap)
        self.assertTrue(analyzed["filled"])
        self.assertEqual((analyzed["fill_direction"], analyzed["report_direction"]), ("bidirectional", "bidirectional"))
        self.assertEqual((analyzed["rounds"], analyzed["fill_length"]), (0, -4))
        self.assertEqual(analyzed["final_sequence_file"], str(self.job / result["final_sequence"]))
        self.assertEqual((analyzed["replace_start"], analyzed["replace_end"]), (0, 19))
        from AutoGapfillerPart4 import GenomeIntegrator
        integrator = GenomeIntegrator.__new__(GenomeIntegrator)
        original = "AAAACCCC" + "N" * 4 + "CCCCGGGG"
        filled = (self.job / result["final_sequence"]).read_text().splitlines()[1]
        self.assertEqual(integrator._build_new_sequence(original, [analyzed], "chr1"), filled)
        internal = dict(analyzed, replace_start=2, replace_end=21)
        self.assertEqual(integrator._build_new_sequence("TT" + original + "AA", [internal], "chr1"),
                         "TT" + filled + "AA")

    def test_auto_job_generation_creates_one_job_without_flag(self):
        manager = JobManager.__new__(JobManager)
        manager.output_dir = self.root / "jobs"
        manager.output_dir.mkdir()
        manager.degap_script = BIN / "DEGAP.py"
        manager.degap_params = {"hifi": str(self.left), "MaximumExtensionRound": 10}
        manager.job_dirs = []
        with patch.object(manager, "_create_symlinks"), patch.object(manager, "_copy_flanking_sequences"):
            manager._create_working_directories({"chr1": {"count": 1, "gaps": [{"id": "gap1"}]}})
        self.assertEqual(len(manager.job_dirs), 1)
        job = manager.job_dirs[0]
        self.assertEqual(job["direction"], "bidirectional")
        script = (Path(job["job_dir"]) / "run.sh").read_text()
        self.assertNotIn("--flag", script)
        self.assertIn("--MaximumExtensionRound 10", script)
        from AutoGapfillerPart3 import JobExecutor
        executor = JobExecutor.__new__(JobExecutor)
        executor.job_config = {"jobs": manager.job_dirs}
        executor.job_list = []
        executor.workers = 1
        with patch.object(executor, "_is_job_completed", return_value=False):
            self.assertEqual(executor._collect_jobs(), 1)
        self.assertEqual(executor.job_list[0]["job"]["direction"], "bidirectional")
        self.assertFalse(executor.job_list[0]["has_right"])

    def test_auto_report_rebases_family_artifact_links(self):
        publish_search(self.search(False))
        visualizer = GapFillingVisualizer.__new__(GapFillingVisualizer)
        visualizer.main_output_dir = self.root
        visualizer.output_dir = self.root / "reports"
        visualizer.output_dir.mkdir()
        visualizer.gap_filling_results = {"gaps": [{"chromosome": "chr1", "gap_id": "gap1", "filled": False,
                                                   "fill_direction": "none", "report_direction": "bidirectional"}]}
        visualizer._copy_subtask_reports()
        copied = visualizer.output_dir / "chr1.gap1.bidirectional.html"
        for link in re.findall(r'href="([^"]+)"', copied.read_text()):
            self.assertTrue((copied.parent / link).exists(), link)

    def test_manifest_metadata_corruption_rejected(self):
        publish_search(self.search())
        path = self.job / "family-result.json"
        result = json.loads(path.read_text())
        result["total_rounds"] = 999
        path.write_text(json.dumps(result))
        with self.assertRaisesRegex(ValueError, "manifest"):
            load_result(self.job)


@unittest.skipUnless(all(shutil.which(t) for t in ("nucmer", "delta-filter", "show-coords")), "MUMmer required")
class PublicGapFillerTests(unittest.TestCase):
    def test_gapfiller_routes_to_families_and_resumes_zero_round_direct_closure(self):
        from Bio import SeqIO
        from GapFiller import GapFiller
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rng = random.Random(193)
            truth = "".join(rng.choices("ACGT", k=22000))
            left, right, reads = root / "left.fa", root / "right.fa", root / "reads.fa"
            left.write_text(">left\n" + truth[:17000] + "\n")
            right.write_text(">right\n" + truth[5000:] + "\n")
            reads.write_text(">read\n" + truth + "\n")
            index = SeqIO.index(str(reads), "fasta")
            self.addCleanup(index.close)
            out = root / "gap"
            params = ["gapfiller", 2, 1, str(reads), str(out), str(left), str(right), "right", 500,
                      None, None, 0, 0, "hifi", None, index, 22000, 15000, None, None, None]
            gap = GapFiller(params, [41, 4, False])
            self.assertEqual(gap.Elongation.search.primary, "L")
            self.assertEqual(gap.Elongation.status, "closed")
            self.assertEqual(gap.Elongation.roundNum, 0)
            self.assertEqual(list(fasta_records(gap.Elongation.finalSeq).values()), [truth])
            again = GapFiller(params, [41, 4, False])
            self.assertEqual(again.Elongation.result, gap.Elongation.result)
            self.assertFalse((out / "process").exists())

    def test_cli_flag_is_ignored_with_notice(self):
        result = subprocess.run([sys.executable, str(BIN / "DEGAP.py"), "--flag", "right", "--help"],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--flag is deprecated and ignored", result.stdout)

    def test_legacy_cli_output_rejected_before_preparation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            old = root / "checkpoint.info"
            old.write_text("round:7")
            result = subprocess.run([sys.executable, str(BIN / "DEGAP.py"), "--mode", "gapfiller", "--out", str(root)],
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Legacy single-sided output", result.stderr)
            self.assertEqual(list(root.iterdir()), [old])
            self.assertEqual(old.read_text(), "round:7")


if __name__ == "__main__":
    unittest.main()
