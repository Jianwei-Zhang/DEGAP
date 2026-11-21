from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from ExtensionFamilies import Family, Version
from FamilyRound import SingleRound
from FamilyTools import shell_command, rg_extract_script
from GapFillerClass import mummer
from FindExtensionReads import FindExtensionReads


class CheckedToolTests(unittest.TestCase):
    def test_failed_upstream_pipeline_is_not_empty_success(self):
        with self.assertRaises(subprocess.CalledProcessError):
            shell_command("false | cat >/dev/null", strict=True)
        self.assertEqual(shell_command("true | cat >/dev/null", strict=True), 0)

    @unittest.skipUnless(shutil.which("minimap2") and shutil.which("samtools"), "mapping tools required")
    def test_mapping_receipt_reuses_complete_bam_and_rejects_corruption(self):
        import pysam
        import FamilyTools
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "reads.fa"
            source.write_text(">read\nACGT\n")
            bam = root / "mapped.bam"
            def execute(command, strict):
                with pysam.AlignmentFile(str(bam), "wb", header={"HD": {"VN": "1.0"}, "SQ": [{"SN": "seed", "LN": 4}]}):
                    pass
                return 0
            with patch.object(FamilyTools, "shell_command", side_effect=execute) as run:
                FamilyTools.mapped_bam("fixture mapping command", bam, [source])
                FamilyTools.mapped_bam("fixture mapping command", bam, [source])
                self.assertEqual(run.call_count, 1)
                bam.write_text("corrupt")
                with self.assertRaisesRegex(ValueError, "mapping changed"):
                    FamilyTools.mapped_bam("fixture mapping command", bam, [source])
                self.assertEqual(run.call_count, 1)

    def test_required_mummer_input_failure_propagates(self):
        with self.assertRaises(FileNotFoundError):
            mummer("/nonexistent-degap-family-input", "also-missing", "unused", strict=True)

    @unittest.skipUnless(shutil.which("rg"), "rg required")
    def test_rg_no_match_is_normal_but_missing_input_is_error(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "kmers").write_text("CCCC\n")
            (root / "reads").write_text(">read\nAAAA\n")
            script = (f'set -e\nKMER_FILE="{root}/kmers"\nPART_FILE="{root}/reads"\n'
                      f'TEMP_OUT="{root}/output"\n' + rg_extract_script("$TEMP_OUT"))
            shell_command(script, strict=True)
            self.assertEqual((root / "output").read_text(), "")
            with self.assertRaises(subprocess.CalledProcessError):
                shell_command(script.replace(f'{root}/reads', f'{root}/missing'), strict=True)


class SingleRoundTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.seed = self.root / "seed.fa"
        self.seed.write_text(">seed\n" + "ACGT" * 10 + "\n")
        self.shared_index = {"read1": "fixture"}
        self.base = SimpleNamespace(out=str(self.root), data_type="hifi", seedLen=20,
                                    hifiSeedLen=20, ontSeedLen=20, thread=1, kmer_size=5,
                                    kmer_num=4, kmer_filter=False, readsDict=self.shared_index,
                                    flag="unchanged", edge=500, reads="fixture.fa")
        self.runner = SingleRound(self.base, context={"fixture": 1})
        self.family = Family("L", "rightward", [Version.capture("L", 0, self.seed)], attempts=1)
        self.calls = []
        self.attempt = self.root / "attempt"
        self.attempt.mkdir()

    def reads(self, seed, last_reads, used_reads, kmer_size, kmer_num, out):
        self.calls.append("reads")
        self.assertIs(seed.elongation.base.readsDict, self.shared_index)
        self.assertIsNone(seed.elongation.base.terminalSeq)
        self.assertTrue(seed.elongation.base.family_mode)
        path = Path(seed.elongation.roundDir) / "reads.fa"
        path.write_text(">read1\n" + "ACGT" * 10 + "AAAAA\n")
        return SimpleNamespace(note="", extensionReads=str(path), extensionReadsNum=1,
                               roundInput=seed)

    def contigs(self, reads, data_type, ont_reads, original_reads_info):
        self.calls.append("assembly")
        seed = reads.roundInput
        from Bio import SeqIO
        sequence = str(SeqIO.read(seed.inputSeq, "fasta").seq)
        sequence = sequence + "AAAAA" if seed.elongation.base.flag == "left" else "AAAAA" + sequence
        path = Path(seed.elongation.roundDir) / "extension.fa"
        path.write_text(">ctg\n" + sequence + "\n")
        return SimpleNamespace(extensionReads=reads, selectContigNote="", extensionLength=5,
                               extensionSequence=str(path), extensionContigID=["ctg"],
                               hifiasmResultDict={"ctg": ["read1"]})

    def run_adapter(self, family=None):
        with patch("FindExtensionReads.FindExtensionReads", side_effect=self.reads), \
             patch("FindExtensionContigs.FindExtensionContigs", side_effect=self.contigs):
            return self.runner(family or self.family, self.attempt)

    def test_left_round_uses_existing_seed_and_splice_without_fixed_terminal(self):
        with patch("GapFillerClass.OutputSequence.linkgap", side_effect=AssertionError("fixed target used")):
            result = self.run_adapter()
        self.assertEqual(result.used_reads, ("read1",))
        self.assertEqual(Path(result.sequence_path).read_text().splitlines()[1], "ACGT" * 10 + "AAAAA")
        self.assertEqual(self.base.flag, "unchanged")
        self.assertFalse(hasattr(self.base, "family_mode"))

    def test_right_round_uses_left_end_and_prepends(self):
        family = Family("R", "leftward", [Version.capture("R", 0, self.seed)], attempts=1)
        result = self.run_adapter(family)
        self.assertEqual(Path(result.sequence_path).read_text().splitlines()[1], "AAAAA" + "ACGT" * 10)

    def test_mixed_different_seed_lengths_merge_against_actual_contig_seed(self):
        self.base.data_type = "mixed"
        self.base.hifiSeedLen = 20
        self.base.ontSeedLen = 30
        result = self.run_adapter()
        self.assertEqual(Path(result.sequence_path).read_text().splitlines()[1], "ACGT" * 10 + "AAAAA")

    def test_short_initial_sequence_merges_without_duplication(self):
        self.base.hifiSeedLen = self.base.seedLen = 100
        result = self.run_adapter()
        self.assertEqual(Path(result.sequence_path).read_text().splitlines()[1], "ACGT" * 10 + "AAAAA")

    def test_identical_extension_is_shared_across_endpoint_searches(self):
        self.runner = SingleRound(self.base, context={"fixture": 1}, shared_cache=self.root / "shared")
        first = self.run_adapter()
        copied_seed = self.root / "copied.fa"
        copied_seed.write_bytes(self.seed.read_bytes())
        other = Family("other-candidate", "rightward", [Version.capture("other-candidate", 0, copied_seed)], attempts=1)
        self.attempt = self.root / "second-search-attempt"
        self.attempt.mkdir()
        second = self.run_adapter(other)
        self.assertEqual(first, second)
        self.assertEqual(self.calls, ["reads", "assembly"])
        self.assertTrue((self.attempt / "shared-stage.json").exists())

    def test_completed_reads_and_assembly_stages_reused(self):
        result = self.run_adapter()
        again = self.run_adapter()
        self.assertEqual(result, again)
        self.assertEqual(self.calls, ["reads", "assembly"])

    def test_assembly_failure_preserves_completed_reads_and_archives_partial_assembly(self):
        def fail(*args):
            (self.attempt / "assembly" / "partial.gfa").write_text("partial")
            raise RuntimeError("assembly tool failure")
        with patch("FindExtensionReads.FindExtensionReads", side_effect=self.reads), \
             patch("FindExtensionContigs.FindExtensionContigs", side_effect=fail):
            with self.assertRaisesRegex(RuntimeError, "tool failure"):
                self.runner(self.family, self.attempt)
        result = self.run_adapter()
        self.assertIsNotNone(result.sequence_path)
        self.assertEqual(self.calls, ["reads", "assembly"])
        partial, = self.attempt.glob("assembly.incomplete.*/partial.gfa")
        self.assertEqual(partial.read_text(), "partial")

    def test_successful_assembler_survives_later_contig_selection_failure(self):
        calls = self.calls
        fail_once = [True]
        class StagedContigs:
            def __init__(self, reads, *args):
                self.extensionReads = reads
                reads.roundInput.elongation.base.family_assembler(self)
                if fail_once[0]:
                    fail_once[0] = False
                    raise RuntimeError("selection interrupted after assembly")
                self.selectContigNote = ""
                self.extensionLength = 5
                self.extensionContigID = ["ctg"]
                self.hifiasmResultDict = {"ctg": ["read1"]}
                path = Path(reads.roundInput.elongation.roundDir) / "extension.fa"
                path.write_text(">ctg\n" + str(reads.roundInput.inputSeedSequence.seq) + "AAAAA\n")
                self.extensionSequence = str(path)
            def hifiasm(self):
                calls.append("hifiasm")
                path = Path(self.extensionReads.roundInput.elongation.roundDir) / "assembled.fa"
                path.write_text(">ctg\nAAAAA\n")
                self.hifiasmResult = str(path)
        with patch("FindExtensionReads.FindExtensionReads", side_effect=self.reads), \
             patch("FindExtensionContigs.FindExtensionContigs", StagedContigs):
            with self.assertRaisesRegex(RuntimeError, "selection interrupted"):
                self.runner(self.family, self.attempt)
            result = self.runner(self.family, self.attempt)
        self.assertEqual(calls, ["reads", "hifiasm"])
        self.assertIsNotNone(result.sequence_path)

    def test_corrupt_completed_reads_are_not_recomputed(self):
        self.run_adapter()
        (self.attempt / "reads" / "reads.fa").write_text("corrupt")
        with self.assertRaisesRegex(ValueError, "artifact changed"):
            self.run_adapter()
        self.assertEqual(self.calls, ["reads", "assembly"])

    def test_read_tool_error_note_is_not_family_exhaustion(self):
        with patch("FindExtensionReads.FindExtensionReads", return_value=SimpleNamespace(note="initializationError")):
            with self.assertRaisesRegex(RuntimeError, "not a normal"):
                self.runner(self.family, self.attempt)

    def test_normal_read_stop_skips_assembly(self):
        with patch("FindExtensionReads.FindExtensionReads", return_value=SimpleNamespace(note="kmerFilterNoReads")), \
             patch("FindExtensionContigs.FindExtensionContigs", side_effect=AssertionError("should not assemble")):
            result = self.runner(self.family, self.attempt)
        self.assertEqual(result.stop_reason, "kmerFilterNoReads")
        self.assertIsNone(result.sequence_path)

    def test_empty_contig_selection_is_a_normal_family_stop(self):
        from FindExtensionContigs import FindExtensionContigs
        def empty_assembly(obj):
            path = Path(obj.extensionReads.roundInput.elongation.roundDir) / "empty.fa"
            path.write_text("")
            obj.hifiasmResult = str(path)
            obj.hifiasmResultDict = {}
        def empty_read_common(obj):
            obj.readCommonResult = obj.hifiasmResult
            obj.extensionSeqNote = "no accurate assembly or read"
        with patch("FindExtensionReads.FindExtensionReads", side_effect=self.reads), \
             patch.object(FindExtensionContigs, "hifiasm", empty_assembly), \
             patch.object(FindExtensionContigs, "readCommon", empty_read_common):
            result = self.runner(self.family, self.attempt)
        self.assertEqual(result.stop_reason, "no_extension_contigs")
        self.assertIsNone(result.sequence_path)

    def test_loop_prevention_is_local_to_family(self):
        self.family.used_reads = ["read1"]
        self.family.last_reads = []
        self.assertEqual(self.run_adapter().stop_reason, "no_new_extension_reads")

    def test_hifiasm_failure_does_not_fall_back_to_apparent_exhaustion(self):
        from FindExtensionContigs import FindExtensionContigs
        obj = FindExtensionContigs.__new__(FindExtensionContigs)
        obj.strict_tools = True
        obj.data_type = "hifi"
        obj.extensionReads = SimpleNamespace(roundInput=SimpleNamespace(
            elongation=SimpleNamespace(roundDir=str(self.attempt), base=self.base)))
        obj.prepareHifiasmReads = lambda: str(self.seed)
        def command(cmd, **kwargs):
            if cmd.startswith("hifiasm -i"):
                return subprocess.CompletedProcess(cmd, 1, "partial", "assembler failure")
            return subprocess.CompletedProcess(cmd, 0, "/fixture/hifiasm", "")
        with patch("FindExtensionContigs.subprocess.run", side_effect=command):
            with self.assertRaises(subprocess.CalledProcessError):
                obj.hifiasm()

    def test_zero_exit_hifiasm_error_is_not_normal_exhaustion(self):
        from FindExtensionContigs import FindExtensionContigs, HifiasmExecutionError
        for strict in (True, False):
            for channel in ("stdout", "stderr"):
                with self.subTest(strict=strict, channel=channel):
                    obj = FindExtensionContigs.__new__(FindExtensionContigs)
                    obj.strict_tools, obj.data_type = strict, "hifi"
                    obj.extensionReads = SimpleNamespace(roundInput=SimpleNamespace(
                        elongation=SimpleNamespace(roundDir=str(self.attempt), base=self.base)))
                    obj.prepareHifiasmReads = lambda: str(self.seed)
                    def command(cmd, **kwargs):
                        if cmd.startswith("hifiasm -i"):
                            error = "[ERROR] input is in fasta format rather than fastq format"
                            return subprocess.CompletedProcess(cmd, 0, error if channel=="stdout" else "",
                                                               error if channel=="stderr" else "")
                        return subprocess.CompletedProcess(cmd, 0, "/fixture/hifiasm", "")
                    with patch("FindExtensionContigs.subprocess.run", side_effect=command):
                        with self.assertRaisesRegex(HifiasmExecutionError, "exit 0"):
                            obj.hifiasm()
                    self.assertFalse(Path(obj.hifiasmResult).exists())

    def test_zero_exit_without_assembly_artifact_is_tool_failure(self):
        from FindExtensionContigs import FindExtensionContigs, HifiasmExecutionError
        obj = FindExtensionContigs.__new__(FindExtensionContigs)
        obj.strict_tools, obj.data_type = True, "hifi"
        obj.extensionReads = SimpleNamespace(roundInput=SimpleNamespace(
            elongation=SimpleNamespace(roundDir=str(self.attempt), base=self.base)))
        obj.prepareHifiasmReads = lambda: str(self.seed)
        with patch("FindExtensionContigs.subprocess.run", return_value=subprocess.CompletedProcess("tool",0,"ok","")):
            with self.assertRaisesRegex(HifiasmExecutionError,"no primary-contig GFA"):
                obj.hifiasm()
        self.assertFalse(Path(obj.hifiasmResult).exists())

    def test_semantic_assembler_failure_does_not_publish_stage_receipt(self):
        from FindExtensionContigs import FindExtensionContigs, HifiasmExecutionError
        with patch("FindExtensionReads.FindExtensionReads",side_effect=self.reads), \
             patch.object(FindExtensionContigs,"hifiasm",side_effect=HifiasmExecutionError("[ERROR] invalid input")):
            with self.assertRaises(HifiasmExecutionError):
                self.runner(self.family,self.attempt)
        self.assertTrue((self.attempt / "reads/stage.json").exists())
        self.assertFalse((self.attempt / "assembler/stage.json").exists())
        self.assertFalse((self.attempt / "assembly/stage.json").exists())

    def test_missing_hifiasm_is_not_a_normal_stop(self):
        from FindExtensionContigs import FindExtensionContigs
        obj = FindExtensionContigs.__new__(FindExtensionContigs)
        obj.strict_tools = True
        obj.data_type = "hifi"
        obj.extensionReads = SimpleNamespace(roundInput=SimpleNamespace(
            elongation=SimpleNamespace(roundDir=str(self.attempt), base=self.base)))
        obj.prepareHifiasmReads = lambda: str(self.seed)
        with patch("FindExtensionContigs.subprocess.run", return_value=subprocess.CompletedProcess("which", 1, "", "")):
            with self.assertRaisesRegex(FileNotFoundError, "hifiasm"):
                obj.hifiasm()

    def test_reads_constructor_propagates_mapping_error_in_family_mode(self):
        elongation = self.runner._elongation(self.family, self.attempt)
        seed = SimpleNamespace(elongation=elongation)
        with patch.object(FindExtensionReads, "minimap2", side_effect=RuntimeError("mapping failure")):
            with self.assertRaisesRegex(RuntimeError, "mapping failure"):
                FindExtensionReads(seed, [], [], 5, 4, out=self.base.out)


@unittest.skipUnless(all(shutil.which(tool) for tool in ("minimap2", "samtools", "hifiasm", "nucmer", "delta-filter", "show-coords")),
                     "Full extension toolchain required")
class RealSingleRoundTests(unittest.TestCase):
    def test_actual_reads_mapping_assembly_and_historical_closure(self):
        from Bio import SeqIO
        from ExtensionFamilies import FamilySearch
        from FamilyOverlap import MummerOverlap, joined_sequence
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rng = random.Random(3913)
            truth = "".join(rng.choices("ACGT", k=60000))
            left, right, reads = root / "left.fa", root / "right.fa", root / "reads.fa"
            left.write_text(">left\n" + truth[:20000] + "\n")
            right.write_text(">right\n" + truth[40000:] + "\n")
            reads.write_text("".join(f">read{i}\n{truth[start:start+20000]}\n"
                                     for i, start in enumerate(range(0, 40001, 1000))))
            index = SeqIO.index(str(reads), "fasta")
            try:
                from GapFiller import GapFiller
                from FamilyResult import load_result
                params = ["gapfiller", 2, 2, str(reads), str(root / "search"), str(left), str(right),
                          "left", 500, None, None, None, 6, "hifi", None,
                          index, 20000, 15000, None, None, None]
                gap = GapFiller(params, [41, 4, False])
                search = gap.Elongation.search
                self.assertEqual(load_result(root / "search")["status"], "closed")
                self.assertEqual(search.status, "closed")
                winner = search.winner
                a = search.families[winner.left_family].versions[winner.left_version]
                b = search.families[winner.right_family].versions[winner.right_version]
                self.assertEqual(joined_sequence(winner, a, b), truth)
                self.assertGreater(search.total_rounds, 0)
                self.assertGreaterEqual(search.families["L"].current.number, 1)
                self.assertGreaterEqual(search.families["R"].current.number, 1)
            finally:
                index.close()


if __name__ == "__main__":
    unittest.main()
