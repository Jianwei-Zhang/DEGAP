from dataclasses import replace
from pathlib import Path
import random
import re
import shutil
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from Bio import SeqIO
from Bio.Seq import Seq
from ExtensionFamilies import AttemptResult, Connection
from FamilyCtgLinker import CtgFamilyPipeline, clip_segments, load_ctg_result, reverse_segments
from FamilyOverlap import read_sequence


class ExactOverlap:
    contract = {"fixture": "exact-suffix-prefix-v1"}

    def __init__(self, minimum=10):
        self.minimum = minimum
        self.calls = []

    def __call__(self, left, right):
        self.calls.append((left.family, left.number, right.family, right.number))
        a, b = read_sequence(left.path), read_sequence(right.path)
        for size in range(min(len(a), len(b)), self.minimum - 1, -1):
            if a[-size:] == b[:size]:
                return [Connection(left.family, left.number, right.family, right.number, 100, size, 0,
                                   len(a)-size, len(a), 0, size)]
        return []


class CtgFamilyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        rng = random.Random(9913)
        self.truth = "".join(rng.choices("ACGT", k=100))
        self.genome = self.root / "genome.fa"
        self.reads = self.root / "reads.fa"
        self.reads.write_text(">fixture-read\n" + self.truth + "\n")
        self.out = self.root / "out"
        self.base = SimpleNamespace(out=str(self.out), genomeSeq=str(self.genome), data_type="hifi",
                                    reads=str(self.reads), ont_reads=None, readsDict={}, ont_readsDict=None,
                                    thread=1, seedLen=1, hifiSeedLen=1, ontSeedLen=None,
                                    edge=5, kmer_size=5, kmer_num=4, kmer_filter=False,
                                    MaximumExtensionRound=0, MaximunExtensionLength=None)

    def inputs(self, sequences):
        self.genome.write_text("".join(f">chr{i}\n{seq}\n" for i, seq in enumerate(sequences)))

    def no_extend(self, family, directory):
        self.fail("Zero-round searches may not run assembly")

    def pipeline(self, overlap=None, extend=None):
        return CtgFamilyPipeline(self.base, overlap=overlap or ExactOverlap(), extend=extend or self.no_extend)

    def sequences(self):
        with (self.out / "DG.scaffold.fa").open() as stream:
            return [str(record.seq) for record in SeqIO.parse(stream, "fasta")]

    def test_two_contigs_join_once_and_agp_reconstructs(self):
        self.inputs([self.truth[:30], self.truth[18:48]])
        result = self.pipeline().run()
        self.assertEqual(self.sequences(), [self.truth[:48]])
        self.assertEqual([m["id"] for m in result["scaffolds"][0]["members"]], ["Contig-0", "Contig-1"])
        self.assertEqual(load_ctg_result(self.out), result)
        self.assertEqual(sum(e["status"] == "closed" for e in result["events"]), 1)
        self.assertTrue(all(e["total_rounds"] == 0 for e in result["events"]))

    def test_reverse_complement_candidate_retains_source_strand(self):
        self.inputs([self.truth[:30], str(Seq(self.truth[18:48]).reverse_complement())])
        result = self.pipeline().run()
        self.assertEqual(self.sequences(), [self.truth[:48]])
        self.assertEqual(result["scaffolds"][0]["members"], [{"id": "Contig-0", "strand": "+"}, {"id": "Contig-1", "strand": "-"}])
        self.assertTrue(any(s["strand"] == "-" for s in result["scaffolds"][0]["segments"]))
        load_ctg_result(self.out)

    def test_three_contigs_do_not_duplicate_seed_or_absorbed_members(self):
        self.inputs([self.truth[:30], self.truth[18:48], self.truth[36:66]])
        result = self.pipeline().run()
        self.assertEqual(self.sequences(), [self.truth[:66]])
        members = result["scaffolds"][0]["members"]
        self.assertEqual([m["id"] for m in members], ["Contig-0", "Contig-1", "Contig-2"])
        self.assertEqual(len(members), len(set(m["id"] for m in members)))
        self.assertEqual(sum(e["status"] == "closed" for e in result["events"]), 2)

    def test_stopped_primary_can_be_closed_by_candidate_extension(self):
        self.inputs([self.truth[40:70], self.truth[:20]])
        self.base.MaximumExtensionRound = 2
        calls = []
        def extend(family, directory):
            calls.append((family.key, family.direction))
            if family.key == "primary":
                return AttemptResult(stop_reason="no_reads")
            path = directory / "extended.fa"
            path.write_text(">candidate\n" + self.truth[:50] + "\n")
            return AttemptResult(str(path), ("fixture-read",))
        result = self.pipeline(extend=extend).run()
        self.assertEqual(self.sequences(), [self.truth[:70]])
        closed = next(e for e in result["events"] if e["status"] == "closed")
        self.assertEqual(closed["total_rounds"], 2)
        self.assertEqual(closed["winner"]["left_version"], 1)
        self.assertEqual(closed["winner"]["right_version"], 0)
        self.assertEqual(calls, [("primary", "leftward"), ("candidate0", "rightward")])
        self.assertEqual(result["scaffolds"][0]["members"], [{"id": "Contig-1", "strand": "+"}, {"id": "Contig-0", "strand": "+"}])
        self.assertTrue(any(v["kind"] == "extension" for v in result["sources"].values()))

    def test_completed_scaffold_can_be_absorbed_without_duplicating_members(self):
        self.inputs([self.truth[:30], self.truth[18:48], self.truth[70:100]])
        self.base.MaximumExtensionRound = 1
        def extend(family, directory):
            if read_sequence(family.current.path) == self.truth[70:100] and family.direction == "leftward":
                path = directory / "extension.fa"
                path.write_text(">extended\n" + self.truth[35:100] + "\n")
                return AttemptResult(str(path))
            return AttemptResult(stop_reason="no_reads")
        result = self.pipeline(extend=extend).run()
        self.assertEqual(self.sequences(), [self.truth])
        self.assertTrue(any(e["selected"] and e["selected"]["assembly"] == "DG1" for e in result["events"]))
        self.assertEqual([m["id"] for m in result["scaffolds"][0]["members"]], ["Contig-0", "Contig-1", "Contig-2"])
        load_ctg_result(self.out)

    def test_reverse_completed_scaffold_composes_original_intervals(self):
        last = str(Seq(self.truth[70:100]).reverse_complement())
        self.inputs([self.truth[:30], self.truth[18:48], last])
        self.base.MaximumExtensionRound = 1
        def extend(family, directory):
            if read_sequence(family.current.path) == last and family.direction == "rightward":
                path = directory / "extension.fa"
                path.write_text(">extended\n" + str(Seq(self.truth[35:100]).reverse_complement()) + "\n")
                return AttemptResult(str(path))
            return AttemptResult(stop_reason="no_reads")
        result = self.pipeline(extend=extend).run()
        self.assertEqual(self.sequences(), [str(Seq(self.truth).reverse_complement())])
        self.assertTrue(any(e["selected"] and e["selected"]["assembly"] == "DG1" and e["selected"]["strand"] == "-"
                            for e in result["events"]))
        self.assertEqual(result["scaffolds"][0]["members"],
                         [{"id": "Contig-2", "strand": "+"}, {"id": "Contig-1", "strand": "-"}, {"id": "Contig-0", "strand": "-"}])
        load_ctg_result(self.out)

    def test_budget_is_shared_per_search_not_given_to_every_candidate(self):
        self.inputs([self.truth[:30], self.truth[40:60], self.truth[70:90]])
        self.base.MaximumExtensionRound = 3
        def extend(family, directory):
            seq = read_sequence(family.current.path)
            seq = seq + "A" if family.direction == "rightward" else "A" + seq
            path = directory / "extended.fa"
            path.write_text(">extended\n" + seq + "\n")
            return AttemptResult(str(path))
        class NoOverlap(ExactOverlap):
            def __call__(self, a, b):
                return []
        result = self.pipeline(overlap=NoOverlap(), extend=extend).run()
        self.assertEqual(len(result["scaffolds"]), 3)
        self.assertTrue(all(e["total_rounds"] == 3 for e in result["events"]))
        self.assertEqual(result["total_rounds"], 18)
        self.assertEqual(self.sequences(), [self.truth[:30], self.truth[40:60], self.truth[70:90]])
        # Unconfirmed histories are retained, not inserted or discarded.
        self.assertTrue(list((self.out / "project").glob("*/attempts/*/*/extended.fa")))

    def test_endpoint_resume_after_closure_before_project_commit_reuses_search(self):
        self.inputs([self.truth[:30], self.truth[18:48]])
        overlap = ExactOverlap()
        pipeline = self.pipeline(overlap=overlap)
        with patch.object(pipeline, "_apply_connection", side_effect=RuntimeError("interrupted before project commit")):
            with self.assertRaisesRegex(RuntimeError, "interrupted"):
                pipeline.run()
        calls = len(overlap.calls)
        result = self.pipeline(overlap=overlap).run()
        self.assertEqual(len(overlap.calls), calls)
        self.assertEqual(self.sequences(), [self.truth[:48]])
        self.assertEqual(len(result["scaffolds"]), 1)

    def test_finished_resume_does_not_recompare_or_rewrite_outputs(self):
        self.inputs([self.truth[:30], self.truth[18:48]])
        overlap = ExactOverlap()
        result = self.pipeline(overlap=overlap).run()
        stamp = (self.out / "DG.scaffold.fa").stat().st_mtime_ns
        calls = len(overlap.calls)
        self.assertEqual(self.pipeline(overlap=overlap).run(), result)
        self.assertEqual(len(overlap.calls), calls)
        self.assertEqual((self.out / "DG.scaffold.fa").stat().st_mtime_ns, stamp)

    def test_unlinked_inputs_remain_separate_and_membership_complete(self):
        self.inputs([self.truth[:20], self.truth[50:70]])
        result = self.pipeline().run()
        self.assertEqual(self.sequences(), [self.truth[:20], self.truth[50:70]])
        self.assertEqual(len(result["scaffolds"]), 2)
        self.assertTrue(all(e["status"] == "round_limit" for e in result["events"]))

    def test_current_19_field_cli_contract_preserves_shared_ont_index(self):
        from CtgLinker import CtgLinker
        index = {"ont-read": object()}
        params = ["ctglinker", 2, 1, str(self.reads), str(self.out), str(self.genome), 500,
                  None, None, None, 2, "mixed", str(self.reads), {}, 20000, 10000, 20000, None, index]
        with patch("FamilyCtgLinker.CtgFamilyPipeline") as driver:
            driver.return_value.run.return_value = {}
            CtgLinker(params, [41, 4, False])
            self.assertIs(driver.call_args.args[0].ont_readsDict, index)

    def test_input_n_splitting_and_seed_filter_follow_existing_contract(self):
        self.inputs([self.truth[:8] + "NNNN" + self.truth[30:60]])
        self.base.hifiSeedLen = self.base.seedLen = 20
        result = self.pipeline().run()
        self.assertEqual(self.sequences(), [self.truth[30:60]])
        self.assertEqual(result["input_contigs"], ["Contig-0"])
        self.assertEqual(result["sources"]["Contig-0"]["origin"], "Contig-0\tchr0")

    def test_seed_filter_falls_back_to_all_contigs_when_none_pass(self):
        self.inputs([self.truth[:20], self.truth[50:70]])
        self.base.hifiSeedLen = self.base.seedLen = 100
        result = self.pipeline().run()
        self.assertEqual(len(result["scaffolds"]), 2)
        self.assertEqual(self.sequences(), [self.truth[:20], self.truth[50:70]])

    def test_no_candidates_do_not_trigger_single_ended_extension(self):
        self.inputs([self.truth[:30]])
        self.base.MaximumExtensionRound = None
        result = self.pipeline().run()
        self.assertEqual(self.sequences(), [self.truth[:30]])
        self.assertEqual([e["status"] for e in result["events"]], ["no_candidates", "no_candidates"])

    def test_input_change_and_result_corruption_rejected(self):
        self.inputs([self.truth[:30]])
        self.pipeline().run()
        self.genome.write_text(">changed\n" + self.truth[1:31] + "\n")
        with self.assertRaisesRegex(ValueError, "changed"):
            self.pipeline().run()
        (self.out / "DG.scaffold.agp").write_text("corrupt")
        with self.assertRaisesRegex(ValueError, "artifact changed"):
            load_ctg_result(self.out)

    def test_legacy_project_is_not_overwritten(self):
        self.inputs([self.truth[:30]])
        self.out.mkdir()
        old = self.out / "used.txt"
        old.write_text("legacy state")
        with self.assertRaisesRegex(ValueError, "Legacy"):
            self.pipeline()
        self.assertEqual(old.read_text(), "legacy state")
        self.assertFalse((self.out / "ctg-families.json").exists())

    def test_checker_and_visualizer_use_verified_family_results(self):
        from CtgLinkerCheck import CtgLinkerCheck
        from CtgLinkerVisualizer import CtgLinkerVisualizer
        self.inputs([self.truth[:30], self.truth[18:48]])
        self.pipeline().run()
        check = CtgLinkerCheck(self.out).run()
        self.assertEqual(len(check["valid_scaffolds"]), 1)
        self.assertEqual(check["valid_scaffolds"][0]["contigs"], ["Contig-0", "Contig-1"])
        report = Path(CtgLinkerVisualizer(self.out).generate_report())
        self.assertEqual(report.name, "ctglinker-family-report.html")
        for link in re.findall(r'href="([^"]+)"', report.read_text()):
            self.assertTrue((report.parent / link).exists(), link)

    def test_interval_clipping_handles_reverse_complement_coordinates(self):
        segments = [{"component": "a", "start": 5, "end": 15, "strand": "+"},
                    {"component": "b", "start": 10, "end": 30, "strand": "-"}]
        self.assertEqual(clip_segments(segments, 12, 18), [{"component": "b", "start": 22, "end": 28, "strand": "-"}])
        self.assertEqual(reverse_segments(reverse_segments(segments)), segments)
        with self.assertRaises(ValueError):
            clip_segments(segments, 0, 31)


@unittest.skipUnless(all(shutil.which(t) for t in ("minimap2", "samtools", "hifiasm", "nucmer", "delta-filter", "show-coords")),
                     "Full extension tools required")
class RealCtgLinkerTests(unittest.TestCase):
    def test_public_ctglinker_extends_candidates_and_reconstructs_gap(self):
        from CtgLinker import CtgLinker
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rng = random.Random(3913)
            truth = "".join(rng.choices("ACGT", k=60000))
            genome, reads = root / "genome.fa", root / "reads.fa"
            genome.write_text(">left\n" + truth[:20000] + "\n>right\n" + truth[40000:] + "\n")
            reads.write_text("".join(f">read{i}\n{truth[start:start+20000]}\n"
                                     for i, start in enumerate(range(0, 40001, 1000))))
            index = SeqIO.index(str(reads), "fasta")
            try:
                out = root / "out"
                params = ["ctglinker", 2, 2, str(reads), str(out), str(genome), 500,
                          None, None, None, 2, "hifi", None, index, 20000, 15000, None, None, None]
                linked = CtgLinker(params, [41, 4, False])
                result = load_ctg_result(out)
                self.assertEqual(len(result["scaffolds"]), 1)
                with open(linked.scaffoldSeq) as stream:
                    self.assertEqual(str(SeqIO.read(stream, "fasta").seq), truth)
                closed = next(e for e in result["events"] if e["status"] == "closed")
                winner = closed["winner"]
                self.assertGreater(winner["left_version"], 0)
                self.assertGreater(winner["right_version"], 0)
                self.assertLessEqual(closed["total_rounds"], 2)
            finally:
                index.close()


if __name__ == "__main__":
    unittest.main()
