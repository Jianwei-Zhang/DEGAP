"""Original input provenance for normalized FASTA and the AutoGapfiller cache."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from AutoGapfillerPart1 import ReadsPreprocessor
from PreparedReads import contract, prepare


@unittest.skipUnless(shutil.which("seqkit"), "seqkit required for prepared-read tests")
class PreparedReadsTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.a = self.root / "a.fa"
        self.b = self.root / "b.fa"
        self.a.write_text(">readA\nACGTACGTACGT\n")
        self.b.write_text(">readB\nTGCATGCATGCA\n")
        self.prepared = self.root / "processed_reads"

    def test_ordered_sources_and_symlinked_cache(self):
        expected = contract([self.a, self.b], None)
        fasta = prepare(self.prepared, expected)["hifi"]
        self.assertLess(fasta.read_text().index(">readA"), fasta.read_text().index(">readB"))
        link = self.root / "shared"
        link.symlink_to(self.prepared, target_is_directory=True)
        stamp = fasta.stat().st_mtime_ns
        with patch("PreparedReads.subprocess.run", side_effect=AssertionError("must not convert again")):
            prepare(link, expected)
            with self.assertRaisesRegex(ValueError, "Original read inputs"):
                prepare(link, contract([self.b, self.a], None))
            with self.assertRaisesRegex(ValueError, "Original read inputs"):
                prepare(link, contract(self.a, self.b))
        self.assertEqual(fasta.stat().st_mtime_ns, stamp)
        fasta.write_text(">other\nCCCCCCCCCCCC\n")
        with self.assertRaisesRegex(ValueError, "Prepared reads changed"):
            prepare(link, expected)

    def test_interrupted_preparation_can_recover_only_with_same_sources(self):
        expected = contract(self.a, self.b)
        with patch("PreparedReads.subprocess.run", side_effect=subprocess.CalledProcessError(1, "seqkit")):
            with self.assertRaises(subprocess.CalledProcessError):
                prepare(self.prepared, expected)
        self.assertIsNone(json.loads((self.prepared / "manifest.json").read_text())["artifacts"])
        self.assertFalse((self.prepared / "hifi_reads.fa").exists())
        with self.assertRaisesRegex(ValueError, "Original read inputs"):
            prepare(self.prepared, contract(self.b, self.a))
        self.assertEqual(set(prepare(self.prepared, expected)), {"hifi", "ont"})

    def test_unbound_legacy_fasta_not_adopted(self):
        self.prepared.mkdir()
        fasta = self.prepared / "hifi_reads.fa"
        fasta.write_text(">older\nACGT\n")
        with self.assertRaisesRegex(ValueError, "lack an input manifest"):
            prepare(self.prepared, contract(self.a, None))
        self.assertFalse((self.prepared / "manifest.json").exists())
        self.assertEqual(fasta.read_text(), ">older\nACGT\n")

    def test_auto_cache_checks_inputs_before_reusing_config(self):
        pre = ReadsPreprocessor(self.root / "auto")
        with (patch.object(pre, "_build_index"), patch.object(pre, "_split_reads"),
              patch.object(pre, "_generate_stats", return_value={})):
            pre.process(self.a)
        config = pre.output_dir / "preprocessing_config.json"
        original_config = config.read_bytes()
        self.assertTrue(ReadsPreprocessor(self.root / "auto").process(self.a)["skipped"])
        self.a.write_text(">different\nCCCCCCCCCCCC\n")
        with self.assertRaisesRegex(ValueError, "preprocessing inputs or depth options changed"):
            ReadsPreprocessor(self.root / "auto").process(self.a)
        self.assertEqual(config.read_bytes(), original_config)
