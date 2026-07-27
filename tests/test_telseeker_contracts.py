import sys
import csv
import json
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
from TelSeekerCtg import TelSeekerCtg, validate_shared_manifest_compatibility
from TelSeekerVisualizer import generate_html_visualization


class TelSeekerContractTests(unittest.TestCase):
    def make_telseeker_parameter(self, tmp_path: Path, tel_read_params=None):
        genome = tmp_path / "genome.fa"
        genome.write_text(">Chr01\nACGT\n")
        original_reads_info = {
            "original_hifi": "reads.fa",
            "original_ont": None,
            "working_hifi": "reads.fa",
            "working_ont": None,
        }
        parameter = [
            "telseeker",
            2,
            "20",
            "reads.fa",
            str(tmp_path),
            str(genome),
            "TTAGGG",
            1,
            500,
            None,
            None,
            None,
            None,
            "hifi",
            None,
            ["Chr01.L"],
        ]
        if tel_read_params is not None:
            parameter.append(tel_read_params)
        parameter.extend([
            "hifi_reads.idx",
            1000,
            100,
            None,
            original_reads_info,
            None,
        ])
        return parameter

    def test_init_reads_terminal_tel_read_parameters(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)

            runner = TelSeeker(
                self.make_telseeker_parameter(
                    tmp_path,
                    {"tel_n": 80, "tel_r": 0.7, "tel_mm": 1},
                ),
                [41, 20, False],
            )

            self.assertEqual(runner.tel_n, 80)
            self.assertEqual(runner.tel_r, 0.7)
            self.assertEqual(runner.tel_mm, 1)

    def test_step0_load_target_ends_writes_part2_input_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp)

            runner = TelSeeker.__new__(TelSeeker)
            runner.out = out_dir
            runner.target_ends = ["Chr01.L", "Chr02.R"]

            runner._step0_load_target_ends()

            target_file = out_dir / "genome.telomere.check" / "need_extension_chr_end.txt"
            self.assertEqual(target_file.read_text(), "Chr01.L\nChr02.R\n")

    def test_run_stops_after_step0_when_no_target_ends_are_available(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp)
            check_dir = out_dir / "genome.telomere.check"
            check_dir.mkdir()
            (check_dir / "need_extension_chr_end.txt").write_text("")

            runner = TelSeeker.__new__(TelSeeker)
            runner.out = out_dir
            calls = []

            runner._step0_load_target_ends = lambda: calls.append("step0")
            runner._step1_extract_telomeric_reads = lambda: calls.append("step1")
            runner._step2_extend_chromosome_ends = lambda: calls.append("step2")
            runner._step3_integrate_results = lambda: calls.append("step3")
            runner._step4_generate_visualization = lambda: calls.append("step4")

            runner.run()

            self.assertEqual(calls, ["step0"])

    def test_run_stops_after_step1_when_no_telomeric_reads_found(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp)
            check_dir = out_dir / "genome.telomere.check"
            part1_dir = out_dir / "part1.telo.reads"
            check_dir.mkdir()
            part1_dir.mkdir()
            (check_dir / "need_extension_chr_end.txt").write_text("Chr01.L\n")

            runner = TelSeeker.__new__(TelSeeker)
            runner.out = out_dir
            calls = []

            def step1():
                calls.append("step1")
                (part1_dir / "Global.left.telo.reads.fa").write_text("")
                (part1_dir / "Global.right.telo.reads.fa").write_text("")

            runner._step0_load_target_ends = lambda: calls.append("step0")
            runner._step1_extract_telomeric_reads = step1
            runner._step2_extend_chromosome_ends = lambda: calls.append("step2")
            runner._step3_integrate_results = lambda: calls.append("step3")
            runner._step4_generate_visualization = lambda: calls.append("step4")

            runner.run()

            self.assertEqual(calls, ["step0", "step1"])

    def test_step1_passes_terminal_tel_read_parameters_to_part1(self):
        import TelSeekerPart1

        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp)
            captured = {}

            class FakeExtractor:
                def __init__(self, **kwargs):
                    captured.update(kwargs)

                def run(self):
                    pass

            original_extractor = TelSeekerPart1.TeloReadsExtractor
            TelSeekerPart1.TeloReadsExtractor = FakeExtractor
            try:
                runner = TelSeeker.__new__(TelSeeker)
                runner.out = out_dir
                runner.motif = "TTAGGG"
                runner.thread = 20
                runner.tel_n = 80
                runner.tel_r = 0.7
                runner.tel_mm = 1

                runner._step1_extract_telomeric_reads()
            finally:
                TelSeekerPart1.TeloReadsExtractor = original_extractor

            self.assertEqual(captured["tel_n"], 80)
            self.assertEqual(captured["tel_r"], 0.7)
            self.assertEqual(captured["tel_mm"], 1)

    def make_telseeker_ctg_parameter(self, tmp_path: Path, ctg_inputs):
        reads = tmp_path / "reads.fa"
        reads.write_text(">read1\nTTAGGGTTAGGG\n")
        return [
            "telseeker_ctg",
            2,
            "4",
            str(reads),
            str(tmp_path),
            ctg_inputs,
            "TTAGGG",
            1,
            500,
            None,
            None,
            1000000,
            0,
            "hifi",
            None,
            {"tel_n": 100, "tel_r": 0.6, "tel_mm": 0},
            {"hifi": [str(reads)], "ont": []},
            str(tmp_path / "hifi_reads.idx"),
            1000,
            100,
            None,
            {
                "original_hifi": [str(reads)],
                "original_ont": None,
                "working_hifi": str(reads),
                "working_ont": None,
            },
            None,
        ]

    def test_telseeker_ctg_manifest_rejects_incompatible_reuse(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            old_reads = tmp_path / "old.fa"
            new_reads = tmp_path / "new.fa"
            old_reads.write_text(">r1\nACGT\n")
            new_reads.write_text(">r1\nACGT\n")
            manifest_dir = tmp_path / "telseeker_ctg.input"
            manifest_dir.mkdir()
            (manifest_dir / "cache_manifest.json").write_text(
                json.dumps({
                    "hifi": [str(old_reads.resolve())],
                    "ont": [],
                    "motif": "TTAGGG",
                    "tel_n": 100,
                    "tel_r": 0.6,
                    "tel_mm": 0,
                })
            )

            with self.assertRaisesRegex(ValueError, "hifi"):
                validate_shared_manifest_compatibility(
                    tmp_path,
                    {"hifi": [str(new_reads)], "ont": []},
                    "TTAGGG",
                    {"tel_n": 100, "tel_r": 0.6, "tel_mm": 0},
                )

    def test_telseeker_ctg_requires_single_record_ctg_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ctg = tmp_path / "ctg.fa"
            ctg.write_text(">ctg1\nACGT\n>ctg2\nACGT\n")
            runner = TelSeekerCtg(
                self.make_telseeker_ctg_parameter(tmp_path, [(str(ctg), "L")]),
                [41, 20, False],
            )

            with self.assertRaisesRegex(ValueError, "exactly one FASTA record"):
                runner._prepare_task_inputs()

    def test_telseeker_ctg_records_failure_when_no_tel_reads_exist(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ctg = tmp_path / "ctg.fa"
            ctg.write_text(">ctg1\n" + "ACGT" * 200 + "\n")
            runner = TelSeekerCtg(
                self.make_telseeker_ctg_parameter(tmp_path, [(str(ctg), "L")]),
                [41, 20, False],
            )
            runner._validate_parameters()
            runner._prepare_task_inputs()
            runner.part1_dir.mkdir()
            (runner.part1_dir / "Global.left.telo.reads.fa").write_text("")
            (runner.part1_dir / "Global.right.telo.reads.fa").write_text("")

            runnable = runner._classify_tel_read_availability()

            self.assertEqual(runnable, [])
            status_file = Path(runner.tasks[0]["job_dir"]) / "result" / "status.tsv"
            with open(status_file, newline="") as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual(rows[0]["status"], "failed")
            self.assertEqual(rows[0]["stop_reason"], "no_tel_reads")

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

    def test_final_genome_analysis_writes_manual_review_evidence_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            part3_dir = Path(tmp) / "part3.integration.results"
            part3_dir.mkdir()
            final_genome = part3_dir / "final.genome.fa"
            final_genome.write_text(">Chr01\n" + "A" * 200 + "\n")

            check_dir = part3_dir / "final.genome.telomere.check"
            check_dir.mkdir()
            for legacy_name in (
                "genome.telomere.check.csv",
                "need_extension_chr_end.txt",
                "uncertain_chr_end.txt",
            ):
                (check_dir / legacy_name).write_text("legacy\n")

            runner = TelSeeker.__new__(TelSeeker)
            runner.motif = "TTAGGG"

            runner._analyze_final_genome(final_genome, part3_dir)

            self.assertTrue(
                (check_dir / "final.genome.telomere.check.left.2kb.fa").exists()
            )
            self.assertTrue(
                (check_dir / "final.genome.telomere.check.right.2kb.fa").exists()
            )
            self.assertFalse((check_dir / "genome.telomere.check.csv").exists())
            self.assertFalse((check_dir / "need_extension_chr_end.txt").exists())
            self.assertFalse((check_dir / "uncertain_chr_end.txt").exists())

    def test_part3_summary_uses_extended_not_telomeric(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            part2_dir = tmp_path / "part2.chr.end.job"
            success_dir = part2_dir / "Chr01.L"
            failed_dir = part2_dir / "Chr01.R"
            success_dir.mkdir(parents=True)
            failed_dir.mkdir(parents=True)
            (success_dir / "linker.fa").write_text(">linker\nACGT\n")
            (success_dir / "linker.info").write_text(
                "Connection Method: direct\nConnected Read: read1\n"
            )

            runner = TelSeeker.__new__(TelSeeker)
            results = runner._collect_part2_results(part2_dir)

            self.assertEqual([row["extended"] for row in results], ["yes", "no"])
            self.assertTrue(all("telomeric" not in row for row in results))

            csv_file = tmp_path / "check_part2_jobs.csv"
            runner._write_results_csv(results, csv_file)
            self.assertEqual(
                csv_file.read_text().splitlines()[0],
                "Chr_end,Extended,Connected_read,Round_num,Connection_type",
            )

    def test_global_report_uses_manual_initial_and_final_review_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp)
            initial_dir = out_dir / "genome.telomere.check"
            final_dir = (
                out_dir
                / "part3.integration.results"
                / "final.genome.telomere.check"
            )
            visual_dir = out_dir / "visual.report"
            initial_dir.mkdir()
            final_dir.mkdir(parents=True)
            visual_dir.mkdir()

            (initial_dir / "need_extension_chr_end.txt").write_text("Chr01.L\n")
            (initial_dir / "genome.telomere.check.left.2kb.fa").write_text(
                ">Chr01_L\nACGT\n"
            )
            (initial_dir / "genome.telomere.check.right.2kb.fa").write_text(
                ">Chr01_R\nACGT\n"
            )
            (initial_dir / "all_chromosomes_combined.png").write_bytes(b"png")
            (final_dir / "final.genome.telomere.check.left.2kb.fa").write_text(
                ">Chr01_L\nACGT\n"
            )
            (final_dir / "final.genome.telomere.check.right.2kb.fa").write_text(
                ">Chr01_R\nACGT\n"
            )
            (final_dir / "all_chromosomes_combined.png").write_bytes(b"png")
            (final_dir / "Chr01_telomere_motif.png").write_bytes(b"png")

            report_file = visual_dir / "Global.report.html"
            generate_html_visualization(out_dir, report_file, "TTAGGG")
            report = report_file.read_text()

            self.assertIn("Initial Chromosome-End Manual Review", report)
            self.assertIn("Final Chromosome-End Manual Review", report)
            self.assertIn("final.genome.telomere.check.left.2kb.fa", report)
            self.assertIn("A successful extension means that a linker was produced", report)
            self.assertNotIn("TRC Score", report)
            self.assertNotIn("Chromosome End Status", report)

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
