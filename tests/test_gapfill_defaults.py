"""The CLI, whole-genome wrapper, and generated single-gap jobs agree."""
import argparse
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from GapFillDefaults import add_options, format_options, resolve_defaults


class GapFillDefaultsTests(unittest.TestCase):
    def parse(self, argv, *, mode=True, auto=False):
        parser = argparse.ArgumentParser()
        add_options(parser, parser, gapfill=True if auto else None)
        args = parser.parse_args(argv)
        if not auto:
            resolve_defaults(args, gapfill=mode)
        return vars(args)

    def test_same_default_for_single_gap_and_genome(self):
        expected = {'MaximumExtensionRound': 10, 'MaximumExtensionLength': 500000, 'kmer_filter': True}
        self.assertEqual(self.parse([]), expected)
        self.assertEqual(self.parse([], auto=True), expected)
        self.assertEqual(self.parse([], mode=False), {'MaximumExtensionRound': None, 'MaximumExtensionLength': None, 'kmer_filter': False})

    def test_generated_values_preserve_zero_unlimited_and_disabled(self):
        options = ['--MaximumExtensionRound', '0', '--MaximumExtensionLength', 'unlimited', '--no-kmer-filter']
        chosen = self.parse(options, auto=True)
        tokens = [arg for option in format_options(chosen) for arg in option]
        self.assertEqual(tokens, options)
        self.assertEqual(self.parse(tokens), chosen)
        self.assertEqual(chosen, {'MaximumExtensionRound': 0, 'MaximumExtensionLength': None, 'kmer_filter': False})

    def test_generated_job_passes_explicit_overrides_to_single_gap_cli(self):
        from AutoGapfillerPart3 import JobManager
        chosen = self.parse(['--MaximumExtensionRound', 'unlimited', '--MaximumExtensionLength', '100', '--no-kmer-filter'], auto=True)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            job = root / 'gap1'
            job.mkdir()
            (root / "01.reads_preprocessor").mkdir()
            detection = root / "02.gap_detection"
            detection.mkdir()
            (detection / "gap_detection_config.json").write_text(json.dumps({"gap_data_dir": str(root)}))
            manager = JobManager(root, chosen)
            manager._generate_run_script(job, "left")
            script, = job.glob('*.sh')
            command = script.read_text()
            self.assertIn('--MaximumExtensionRound unlimited', command)
            self.assertIn('--MaximumExtensionLength 100', command)
            self.assertIn('--no-kmer-filter', command)
            self.assertNotIn('--kmer_filter', command)
