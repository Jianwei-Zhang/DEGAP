import sys
import tempfile
import unittest
import importlib.util
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from TelSeeker import TelSeeker
from TelSeekerCheck import calculate_trc as check_calculate_trc
from TelSeekerCheck import patterns_to_search_topsicle as check_patterns_to_search_topsicle
from TelSeekerPart1 import calculate_trc as part1_calculate_trc
from TelSeekerPart1 import TeloReadsExtractor
from TelSeekerPart2 import finalize_chr_end_result


class TelSeekerContractTests(unittest.TestCase):
    def test_parse_linker_info_reads_extension_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            info_file = Path(tmp) / "linker.info"
            info_file.write_text(
                "\n".join(
                    [
                        "Chromosome End: Chr01.L",
                        "Connection Method: extension",
                        "Source File: source.fa",
                        "Final File: linker.fa",
                        "Status: Extension connection found",
                        "Extension Rounds: 3",
                        "Total Extension Length: 1200 bp",
                        "Connected Read: read42",
                    ]
                )
                + "\n"
            )

            parser = TelSeeker.__new__(TelSeeker)
            result = parser._parse_linker_info(info_file)

            self.assertEqual(result["connection_type"], "extension")
            self.assertEqual(result["round_num"], "3")
            self.assertEqual(result["connected_read"], "read42")

    def test_finalize_preserves_direct_connection_method(self):
        with tempfile.TemporaryDirectory() as tmp:
            work_path = Path(tmp) / "Chr01.L"
            direct_dir = work_path / "direct.check"
            direct_dir.mkdir(parents=True)
            source = direct_dir / "direct.linker.fa"
            source.write_text(">direct\nACGT\n")

            ok = finalize_chr_end_result(
                "Chr01.L",
                work_path,
                {
                    "success": True,
                    "connection_info": {
                        "method": "direct",
                        "read_id": "read-direct",
                        "linker_file": str(source),
                    },
                },
            )

            self.assertTrue(ok)
            self.assertEqual((work_path / "linker.fa").read_text(), source.read_text())
            info = (work_path / "linker.info").read_text()
            self.assertIn("Connection Method: direct\n", info)
            self.assertIn("Status: Direct connection found\n", info)
            self.assertIn("Connected Read: read-direct\n", info)

    def test_step2_complete_requires_all_expected_end_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp)
            check_dir = out_dir / "genome.telomere.check"
            check_dir.mkdir()
            (check_dir / "need_extension_chr_end.txt").write_text("Chr01.L\nChr01.R\n")

            part2_dir = out_dir / "part2.chr.end.job"
            complete_end = part2_dir / "Chr01.L"
            complete_end.mkdir(parents=True)
            (complete_end / "linker.fa").write_text(">linker\nACGT\n")
            (complete_end / "linker.info").write_text("Connection Method: direct\n")

            checker = TelSeeker.__new__(TelSeeker)
            checker.out = out_dir

            self.assertFalse(checker._check_step2_complete(part2_dir))

    def test_step2_complete_rejects_direct_check_only_without_terminal_result(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp)
            check_dir = out_dir / "genome.telomere.check"
            check_dir.mkdir()
            (check_dir / "need_extension_chr_end.txt").write_text("Chr01.L\n")

            part2_dir = out_dir / "part2.chr.end.job"
            direct_dir = part2_dir / "Chr01.L" / "direct.check"
            direct_dir.mkdir(parents=True)
            (direct_dir / "direct.check.log").write_text(
                "Status: FAILED - No direct connection found\n"
            )

            checker = TelSeeker.__new__(TelSeeker)
            checker.out = out_dir

            self.assertFalse(checker._check_step2_complete(part2_dir))

    def test_part1_trc_matches_telseeker_check_trc(self):
        sequence = "TTAGGGTTAGGGTTAGGG"
        motif = "TTAGGG"
        kmer_length = len(motif) - 2

        self.assertEqual(
            part1_calculate_trc(sequence, motif, kmer_length, len(sequence)),
            check_calculate_trc(sequence, motif, kmer_length, len(sequence)),
        )

    def test_part1_generated_worker_uses_telseeker_check_pattern_semantics(self):
        with tempfile.TemporaryDirectory() as tmp:
            extractor = TeloReadsExtractor(tmp, "TTAGGG")
            script_path = Path(extractor._create_processing_script(tmp))

            spec = importlib.util.spec_from_file_location("part1_worker", script_path)
            worker = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(worker)

            self.assertEqual(
                worker.patterns_to_search_topsicle("TTAGGG", 4),
                check_patterns_to_search_topsicle("TTAGGG", 4),
            )


if __name__ == "__main__":
    unittest.main()
