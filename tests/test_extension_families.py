import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import ExtensionFamilies as ef
from ExtensionFamilies import AttemptResult, Connection, FamilySearch


class FamilySearchTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.calls = []
        self.pairs = []
        self.seeds = [("L", "rightward", self.fasta("L", 40)),
                      ("R", "leftward", self.fasta("R", 40))]

    def fasta(self, name, length):
        path = self.root / (name + ".fa")
        path.write_text(">" + name + "\n" + "A" * length + "\n")
        return str(path)

    def search(self, **kwargs):
        seeds = kwargs.pop("seeds", self.seeds)
        primary = kwargs.pop("primary", "L")
        context = kwargs.pop("context", {"reads_sha256": "fixture", "adapter": 1})
        directory = kwargs.pop("directory", self.root / "search")
        return FamilySearch(directory, seeds, primary=primary, context=context, **kwargs)

    def extend(self, family, directory):
        self.calls.append((family.key, family.attempts, tuple(family.used_reads)))
        path = directory / "output.fa"
        path.write_text(">extended\n" + "A" * (family.current.length + 10) + "\n")
        return AttemptResult(str(path), ("shared-read", f"{family.key}-{family.attempts}"))

    def compare(self, left, right):
        pair = (left.family, left.number, right.family, right.number)
        self.assertNotIn(pair, self.pairs, "A version pair must not be recomputed")
        self.pairs.append(pair)
        return []

    def hit(self, left, right, identity=99.5, overlap=20, overhang=0):
        return Connection(left.family, left.number, right.family, right.number,
                          identity, overlap, overhang,
                          left.length - overlap, left.length, 0, overlap)

    def test_alternating_shared_rounds_and_independent_reads(self):
        search = self.search(max_rounds=4).run(self.extend, self.compare)
        self.assertEqual([c[0] for c in self.calls], ["L", "R", "L", "R"])
        self.assertEqual(self.calls[1][2], ())
        self.assertNotIn("L-1", self.calls[3][2])
        self.assertIn("shared-read", self.calls[3][2])
        self.assertEqual((search.total_rounds, search.total_length, search.status), (4, 40, "round_limit"))
        self.assertEqual(len(self.pairs), 9)

    def test_stopped_side_keeps_history_and_peer_closes(self):
        def extend(family, directory):
            if family.key == "L" and family.attempts == 2:
                self.calls.append((family.key, family.attempts))
                return AttemptResult(stop_reason="no_extension_reads")
            return self.extend(family, directory)
        def compare(left, right):
            self.compare(left, right)
            return [self.hit(left, right)] if (left.number, right.number) == (1, 3) else []
        search = self.search().run(extend, compare)
        self.assertEqual([c[0] for c in self.calls], ["L", "R", "L", "R", "R"])
        self.assertEqual(search.status, "closed")
        self.assertEqual(search.families["L"].stop_reason, "no_extension_reads")
        self.assertEqual((search.total_rounds, search.total_length), (5, 40))
        self.assertEqual(search.winner.left_version, 1)

    def test_l3_connects_r2_not_already_existing_r3(self):
        def compare(left, right):
            self.compare(left, right)
            return [self.hit(left, right)] if (left.number, right.number) == (3, 2) else []
        search = self.search(primary="R").run(self.extend, compare)
        self.assertEqual(search.families["R"].current.number, 3)
        self.assertIn(("L", 3, "R", 3), self.pairs)
        self.assertEqual((search.winner.left_version, search.winner.right_version), (3, 2))
        self.assertEqual((search.total_rounds, search.total_length), (6, 60))

    def test_initial_overlap_needs_no_rounds_even_with_zero_budgets(self):
        search = self.search(max_rounds=0, max_length=0).run(
            self.extend, lambda l, r: [self.hit(l, r)])
        self.assertEqual(search.status, "closed")
        self.assertEqual(self.calls, [])
        self.assertEqual((search.total_rounds, search.total_length), (0, 0))

    def test_zero_budget_without_overlap(self):
        search = self.search(max_rounds=0).run(self.extend, self.compare)
        self.assertEqual(search.status, "round_limit")
        self.assertEqual(len(self.pairs), 1)
        self.assertEqual(self.calls, [])

    def test_length_limit_permits_last_whole_round_then_checks_closure(self):
        def compare(left, right):
            self.compare(left, right)
            return [self.hit(left, right)] if right.number == 1 else []
        search = self.search(max_length=15).run(self.extend, compare)
        self.assertEqual((search.status, search.total_length, search.total_rounds), ("closed", 20, 2))

    def test_length_limit_without_closure_stops_before_next_round(self):
        search = self.search(max_length=15).run(self.extend, self.compare)
        self.assertEqual((search.status, search.total_length, search.total_rounds), ("length_limit", 20, 2))

    def test_both_stopped_is_not_closed(self):
        search = self.search().run(lambda f, d: AttemptResult(stop_reason="no_contig"), self.compare)
        self.assertEqual((search.status, search.total_rounds), ("exhausted", 2))
        self.assertIsNone(search.winner)
        self.assertTrue(all(len(f.versions) == 1 for f in search.families.values()))

    def test_zero_growth_stops_only_that_family_without_new_version(self):
        def extend(family, directory):
            if family.key == "L":
                return AttemptResult(family.current.path)
            return self.extend(family, directory)
        search = self.search(max_rounds=3).run(extend, self.compare)
        self.assertEqual(search.families["L"].stop_reason, "no_growth")
        self.assertEqual(len(search.families["L"].versions), 1)
        self.assertEqual(search.total_length, 20)

    def test_all_candidates_checked_before_selection(self):
        seeds = self.seeds + [("C", "leftward", self.fasta("C", 40))]
        def compare(left, right):
            self.compare(left, right)
            return [self.hit(left, right, identity=99.9 if right.family == "C" else 99.0)]
        search = self.search(seeds=seeds).run(self.extend, compare)
        self.assertEqual(search.winner.right_family, "C")
        self.assertEqual(len(self.pairs), 2)
        self.assertFalse(self.calls)

    def test_ranking_identity_then_length_then_overhang_then_stable_coordinates(self):
        search = self.search()
        left, right = search.families["L"].current, search.families["R"].current
        candidates = [self.hit(left, right, 99.5, 30, 0),
                      self.hit(left, right, 99.9, 12, 0),
                      self.hit(left, right, 99.9, 20, 2),
                      self.hit(left, right, 99.9, 20, 0)]
        search.run(self.extend, lambda l, r: list(reversed(candidates)))
        self.assertEqual(search.winner, candidates[-1])

    def test_same_quality_history_prefers_fixed_version_order_not_visit_order(self):
        def compare(left, right):
            return [self.hit(left, right)] if left.number == 2 else []
        search = self.search().run(self.extend, compare)
        self.assertEqual((search.winner.left_version, search.winner.right_version), (2, 0))

    def test_multi_candidate_budget_and_no_candidate_candidate_checks(self):
        seeds = self.seeds + [("C", "leftward", self.fasta("C", 40))]
        search = self.search(seeds=seeds, max_rounds=6).run(self.extend, self.compare)
        self.assertEqual([c[0] for c in self.calls], ["L", "R", "L", "C", "L", "R"])
        self.assertTrue(all(p[0] == "L" for p in self.pairs))
        self.assertEqual(search.total_length, 60)

    def test_primary_exhaustion_leaves_candidates_round_robin(self):
        seeds = self.seeds + [("C", "leftward", self.fasta("C", 40))]
        def extend(family, directory):
            if family.key == "L":
                self.calls.append((family.key, family.attempts))
                return AttemptResult(stop_reason="no_reads")
            return self.extend(family, directory)
        self.search(seeds=seeds, max_rounds=5).run(extend, self.compare)
        self.assertEqual([c[0] for c in self.calls], ["L", "R", "C", "R", "C"])

    def test_single_ended_terminal_has_no_fake_self_comparison(self):
        def extend(family, directory):
            result = self.extend(family, directory)
            return AttemptResult(result.sequence_path, result.used_reads, terminal_reason="telomere")
        search = self.search(seeds=self.seeds[:1]).run(extend, self.compare)
        self.assertEqual((search.status, search.total_rounds), ("terminal", 1))
        self.assertFalse(self.pairs)
        self.assertIsNone(search.winner)

    def test_tool_error_is_pending_not_normal_exhaustion(self):
        search = self.search(max_rounds=1)
        with self.assertRaisesRegex(RuntimeError, "assembler failed"):
            search.run(lambda f, d: (_ for _ in ()).throw(RuntimeError("assembler failed")), self.compare)
        resumed = self.search(max_rounds=1)
        self.assertEqual(resumed.total_rounds, 1)
        self.assertIsNotNone(resumed.pending)
        self.assertIsNone(resumed.families["L"].stop_reason)
        resumed.run(self.extend, self.compare)
        self.assertEqual((resumed.total_rounds, resumed.total_length), (1, 10))

    def test_resume_after_receipt_before_checkpoint_does_not_repeat_assembly(self):
        original = ef.atomic_json
        def interrupted(path, value):
            if Path(path).name == "families.json" and value["completed_attempts"] == 1:
                raise OSError("simulated checkpoint interruption")
            return original(path, value)
        search = self.search(max_rounds=2)
        with patch.object(ef, "atomic_json", side_effect=interrupted):
            with self.assertRaisesRegex(OSError, "simulated"):
                search.run(self.extend, self.compare)
        resumed = self.search(max_rounds=2).run(self.extend, self.compare)
        self.assertEqual([c[0] for c in self.calls], ["L", "R"])
        self.assertEqual((resumed.total_rounds, resumed.total_length), (2, 20))

    def test_partial_family_comparison_resumes_without_early_wrong_winner(self):
        seeds = self.seeds + [("C", "leftward", self.fasta("C", 40))]
        def interrupted(left, right):
            if right.family == "C":
                raise RuntimeError("comparison interrupted")
            return [self.hit(left, right, 99)]
        with self.assertRaises(RuntimeError):
            self.search(seeds=seeds).run(self.extend, interrupted)
        def resume_compare(left, right):
            self.assertEqual(right.family, "C")
            return [self.hit(left, right, 100)]
        search = self.search(seeds=seeds).run(self.extend, resume_compare)
        self.assertEqual(search.winner.right_family, "C")
        self.assertEqual(search.total_rounds, 0)

    def test_completed_search_resume_does_no_work(self):
        self.search(max_rounds=2).run(self.extend, self.compare)
        calls = len(self.calls)
        self.search(max_rounds=2).run(self.extend, self.compare)
        self.assertEqual(len(self.calls), calls)

    def test_context_or_budget_change_rejected_without_overwrite(self):
        search = self.search(max_rounds=1)
        before = search.checkpoint.read_bytes()
        for kwargs in ({"max_rounds": 2}, {"max_rounds": 1, "context": {"reads_sha256": "changed"}}):
            with self.assertRaisesRegex(ValueError, "changed"):
                self.search(**kwargs)
            self.assertEqual(search.checkpoint.read_bytes(), before)

    def test_history_content_change_rejected(self):
        search = self.search(max_rounds=1).run(self.extend, self.compare)
        Path(search.families["L"].current.path).write_text(">changed\n" + "C" * 50 + "\n")
        with self.assertRaisesRegex(ValueError, "changed"):
            self.search(max_rounds=1)

    def test_invalid_checkpoint_not_silently_restarted(self):
        search = self.search()
        search.checkpoint.write_text("{broken")
        with self.assertRaisesRegex(ValueError, "receipt"):
            self.search()
        self.assertEqual(search.checkpoint.read_text(), "{broken")

    def test_corrupt_budget_rejected(self):
        search = self.search(max_rounds=1).run(self.extend, self.compare)
        state = json.loads(search.checkpoint.read_text())
        state["families"][0]["added_bases"] = 999
        search.checkpoint.write_text(json.dumps(state))
        with self.assertRaisesRegex(ValueError, "budget"):
            self.search(max_rounds=1)

    def test_legacy_task_rejected_without_mutation(self):
        directory = self.root / "legacy"
        directory.mkdir()
        old = directory / "checkpoint.info"
        old.write_text("round:7")
        with self.assertRaisesRegex(ValueError, "Legacy"):
            self.search(directory=directory)
        self.assertEqual(old.read_text(), "round:7")
        self.assertFalse((directory / "families.json").exists())

    def test_adapter_cannot_mutate_other_state(self):
        def extend(family, directory):
            family.used_reads.append("bad-side-effect")
            family.attempts = 999
            return AttemptResult(stop_reason="normal_stop")
        search = self.search(max_rounds=1).run(extend, self.compare)
        self.assertEqual(search.total_rounds, 1)
        self.assertEqual(search.families["L"].used_reads, [])

    def test_live_search_lock_rejects_second_constructor(self):
        import fcntl
        search = self.search()
        with (search.directory / "families.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            with self.assertRaisesRegex(RuntimeError, "already running"):
                self.search()

    def test_corrupted_pending_result_is_not_accepted(self):
        original = ef.atomic_json
        def interrupted(path, value):
            if Path(path).name == "families.json" and value["completed_attempts"] == 1:
                raise OSError("interrupt")
            return original(path, value)
        search = self.search(max_rounds=1)
        with patch.object(ef, "atomic_json", side_effect=interrupted):
            with self.assertRaises(OSError):
                search.run(self.extend, self.compare)
        output = search.directory / "attempts" / "0" / "1" / "output.fa"
        output.write_text(">corrupt\nT\n")
        resumed = self.search(max_rounds=1)
        with self.assertRaisesRegex(ValueError, "changed"):
            resumed.run(self.extend, self.compare)
        self.assertEqual(len(self.calls), 1)

    def test_growing_and_stopping_result_still_checks_closure(self):
        def extend(family, directory):
            result = self.extend(family, directory)
            return AttemptResult(result.sequence_path, result.used_reads, stop_reason="no_new_reads")
        def compare(left, right):
            return [self.hit(left, right)] if left.number == 1 else []
        search = self.search().run(extend, compare)
        self.assertEqual((search.status, search.total_rounds), ("closed", 1))
        self.assertEqual(search.families["L"].stop_reason, "no_new_reads")

    def test_invalid_limits_and_directions_rejected(self):
        for value in (-1, 0.5, True):
            with self.assertRaises(ValueError):
                self.search(max_rounds=value)
        with self.assertRaisesRegex(ValueError, "face"):
            self.search(seeds=[self.seeds[0], ("R", "rightward", self.seeds[1][2])])

    def test_shortened_result_rejected(self):
        with self.assertRaisesRegex(ValueError, "shortened"):
            self.search().run(lambda f, d: AttemptResult(self.fasta("short", 10)), self.compare)

    def test_invalid_candidate_and_fake_terminal_rejected(self):
        with self.assertRaisesRegex(ValueError, "identity"):
            self.search().run(self.extend, lambda l, r: [self.hit(l, r, float("nan"))])
        with self.assertRaisesRegex(ValueError, "overlap evidence"):
            self.search().run(lambda f, d: AttemptResult(terminal_reason="fake"), self.compare)


if __name__ == "__main__":
    unittest.main()
