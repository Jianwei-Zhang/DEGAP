from dataclasses import replace
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from ExtensionFamilies import AttemptResult, FamilySearch, Version
from FamilyOverlap import MummerOverlap, joined_sequence, parse_coords
import FamilyOverlap

TOOLS_AVAILABLE = all(shutil.which(name) for name in ("nucmer", "delta-filter", "show-coords"))


class CoordsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.left = Version("L", 1, "unused-left", "fixture", 50000)
        self.right = Version("R", 1, "unused-right", "fixture", 40000)
        self.path = self.root / "coords"

    def row(self, *, a=1, b=12000, c=8501, d=20500, identity=99.9, names=("right", "left")):
        return "\t".join(map(str, (a, b, c, d, abs(b-a)+1, abs(d-c)+1,
                                   identity, 20500, 20500, 0, 0, *names))) + "\n"

    def parse(self, text, **kwargs):
        self.path.write_text(text)
        return parse_coords(self.path, kwargs.pop("left", self.left), kwargs.pop("right", self.right),
                            left_offset=29500, left_tail_length=20500, right_tail_length=20500,
                            edge=500, **kwargs)

    def test_coordinates_restore_full_sequence_offsets(self):
        hit, = self.parse(self.row())
        self.assertEqual((hit.left_start, hit.left_end, hit.right_start, hit.right_end),
                         (38000, 50000, 0, 12000))
        self.assertEqual((hit.overlap, hit.overhang), (12000, 0))

    def test_reverse_and_nonterminal_alignments_are_not_connections(self):
        self.assertEqual(self.parse(self.row(c=20500, d=8501)), [])
        self.assertEqual(self.parse(self.row(a=1001, b=13000)), [])
        self.assertEqual(self.parse(self.row(c=7500, d=19499)), [])

    def test_short_overlap_is_not_accepted(self):
        self.assertEqual(self.parse(self.row(b=9999, c=10502)), [])

    def test_initial_identity_gate_does_not_silently_change_extension_gate(self):
        text = self.row(identity=79.9)
        self.assertEqual(self.parse(text, left=replace(self.left, number=0),
                                    right=replace(self.right, number=0)), [])
        self.assertEqual(len(self.parse(text)), 1)

    def test_malformed_coords_are_errors_not_no_overlap(self):
        for text in ("broken\n", self.row(names=("wrong", "left")), self.row(identity=float("nan")),
                     self.row(b=30000)):
            with self.subTest(text=text), self.assertRaisesRegex(ValueError, "coordinates"):
                self.parse(text)

    def test_empty_coords_mean_valid_no_overlap(self):
        self.assertEqual(self.parse(""), [])

    def test_join_uses_coordinates_not_overlap_length(self):
        a, b = self.root / "left.fa", self.root / "right.fa"
        a.write_text(">left\nAAAACCCCGG\n")
        b.write_text(">right\nTCCCCGGGG\n")
        left, right = Version.capture("L", 1, a), Version.capture("R", 2, b)
        from ExtensionFamilies import Connection
        hit = Connection("L", 1, "R", 2, 99, 4, 3, 4, 8, 1, 5)
        self.assertEqual(joined_sequence(hit, left, right), "AAAACCCCGGGG")


@unittest.skipUnless(TOOLS_AVAILABLE, "MUMmer tools are required for the real alignment fixtures")
class RealMummerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        rng = random.Random(913)
        self.truth = "".join(rng.choices("ACGT", k=52000))
        self.left = self.write("L", self.truth[:6000])
        self.right = self.write("R", self.truth[-6000:])
        self.overlap = MummerOverlap(self.root / "overlaps", 500)

    def write(self, name, sequence):
        path = self.root / (name + ".fa")
        path.write_text(f">{name}\n{sequence}\n")
        return str(path)

    def test_real_both_sides_overlap_reconstructs_known_sequence(self):
        calls = []
        def extend(family, directory):
            calls.append(family.key)
            sequence = self.truth[:31000] if family.key == "L" else self.truth[-31000:]
            path = directory / "extension.fa"
            path.write_text(">extension\n" + sequence + "\n")
            return AttemptResult(str(path), ("fixture-read",))
        search = FamilySearch(self.root / "search", [("L", "rightward", self.left),
                                                     ("R", "leftward", self.right)],
                              primary="L", context=self.overlap.contract)
        search.run(extend, self.overlap)
        self.assertEqual(calls, ["L", "R"])
        self.assertEqual((search.status, search.total_length), ("closed", 50000))
        self.assertEqual(joined_sequence(search.winner, search.families["L"].current,
                                        search.families["R"].current), self.truth)

    def test_real_l3_r2_closure_despite_bad_r3(self):
        def extend(family, directory):
            length = 6000 + family.attempts * 10000
            sequence = self.truth[:length] if family.key == "L" else self.truth[-length:]
            if family.key == "R" and family.attempts == 3:
                sequence = "A" * 10000 + self.truth[-26000:]
            path = directory / "extension.fa"
            path.write_text(">extension\n" + sequence + "\n")
            return AttemptResult(str(path))
        search = FamilySearch(self.root / "search", [("L", "rightward", self.left),
                                                     ("R", "leftward", self.right)],
                              primary="R", context=self.overlap.contract, max_rounds=6)
        search.run(extend, self.overlap)
        self.assertEqual((search.status, search.winner.left_version, search.winner.right_version),
                         ("closed", 3, 2))
        self.assertEqual(search.families["R"].current.number, 3)
        self.assertEqual(joined_sequence(search.winner, search.families["L"].versions[3],
                                        search.families["R"].versions[2]), self.truth)

    def test_completed_stages_reused_after_show_coords_failure(self):
        left = Version.capture("L", 1, self.write("extendedL", self.truth[:31000]))
        right = Version.capture("R", 1, self.write("extendedR", self.truth[-31000:]))
        real_run = subprocess.run
        calls = []
        failed = False
        def run(command, **kwargs):
            nonlocal failed
            name = Path(command[0]).name
            calls.append(name)
            if name == "show-coords" and not failed:
                failed = True
                raise subprocess.CalledProcessError(1, command)
            return real_run(command, **kwargs)
        with patch.object(FamilyOverlap.subprocess, "run", side_effect=run):
            with self.assertRaises(subprocess.CalledProcessError):
                self.overlap(left, right)
            matches = self.overlap(left, right)
        self.assertTrue(matches)
        self.assertEqual(calls.count("nucmer"), 1)
        self.assertEqual(calls.count("delta-filter"), 1)
        self.assertEqual(calls.count("show-coords"), 2)
        with patch.object(FamilyOverlap.subprocess, "run", side_effect=AssertionError("unexpected rerun")):
            self.assertEqual(self.overlap(left, right), matches)

    def test_alignment_reused_when_identical_sequences_change_family_labels(self):
        left = Version.capture("primary", 1, self.write("extendedL", self.truth[:31000]))
        right = Version.capture("candidate0", 1, self.write("extendedR", self.truth[-31000:]))
        first = self.overlap(left, right)
        with patch.object(FamilyOverlap.subprocess, "run", side_effect=AssertionError("unnecessary remapping")):
            second = self.overlap(replace(left, family="candidate2"), replace(right, family="primary"))
        self.assertEqual(len(first), len(second))
        self.assertEqual(second[0].left_family, "candidate2")
        self.assertEqual(second[0].right_family, "primary")
        self.assertEqual(second[0].identity, first[0].identity)

    def test_corrupt_completed_stage_is_not_silently_reused_or_overwritten(self):
        left = Version.capture("L", 1, self.write("extendedL", self.truth[:31000]))
        right = Version.capture("R", 1, self.write("extendedR", self.truth[-31000:]))
        self.overlap(left, right)
        coords, = (self.root / "overlaps").glob("*/alignment.coords")
        coords.write_text("corrupted")
        with self.assertRaisesRegex(ValueError, "completed overlap stage"):
            self.overlap(left, right)
        self.assertEqual(coords.read_text(), "corrupted")


if __name__ == "__main__":
    unittest.main()
