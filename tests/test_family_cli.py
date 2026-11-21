"""Consolidated real-tool acceptance fixtures (not a biological benchmark)."""
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

from Bio import SeqIO

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))
TOOLS = ("seqkit", "minimap2", "samtools", "hifiasm", "nucmer", "delta-filter", "show-coords")


@unittest.skipUnless(all(shutil.which(tool) for tool in TOOLS), "Full CLI toolchain required")
class RealFamilyCLITests(unittest.TestCase):
    def setUp(self):
        retained = os.environ.get("DEGAP_VALIDATION_ARTIFACTS")
        if retained:
            Path(retained).mkdir(parents=True, exist_ok=True)
            self.root = Path(tempfile.mkdtemp(prefix=self._testMethodName + "-", dir=retained))
        else:
            tmp = tempfile.TemporaryDirectory()
            self.addCleanup(tmp.cleanup)
            self.root = Path(tmp.name)
        rng = random.Random(3913)
        self.truth = "".join(rng.choices("ACGT", k=60000))
        self.left = self.root / "left.fa"
        self.right = self.root / "right.fa"
        self.left.write_text(">left\n" + self.truth[:20000] + "\n")
        self.right.write_text(">right\n" + self.truth[40000:] + "\n")
        for platform in ("hifi", "ont"):
            (self.root / (platform + ".fa")).write_text("".join(
                f">{platform}{i}\n{self.truth[start:start+20000]}\n"
                for i, start in enumerate(range(0, 40001, 1000))))

    def command(self, script, args, log="command.log"):
        path = self.root / log
        with path.open("w") as stream:
            result = subprocess.run([sys.executable, str(BIN / script)] + args,
                                    stdout=stream, stderr=subprocess.STDOUT, timeout=300)
        self.assertEqual(result.returncode, 0, path.read_text()[-12000:])

    def gap(self, platforms):
        out = self.root / "gap"
        args = ["--mode", "gapfiller", "--out", str(out), "--seqleft", str(self.left),
                "--seqright", str(self.right), "--MaximumExtensionRound", "6", "-t", "2", "--remove", "2"]
        for platform in platforms:
            args += ["--" + platform, str(self.root / (platform + ".fa"))]
        self.command("DEGAP.py", args)
        from FamilyResult import load_result
        result = load_result(out)
        self.assertEqual(result["status"], "closed")
        sequences = list(SeqIO.parse(out / result["final_sequence"], "fasta"))
        self.assertEqual([str(record.seq) for record in sequences], [self.truth])
        return out, args

    def test_gap_cli_hifi_and_resume(self):
        out, args = self.gap(["hifi"])
        self.assertTrue((out / "processed_reads" / "hifi_reads.fa").is_file())
        checkpoint = (out / "families.json").read_bytes()
        stages = {str(p): p.stat().st_mtime_ns for p in out.rglob("stage.json")}
        self.command("DEGAP.py", args + ["--flag", "right"], log="resume.log")
        self.assertEqual((out / "families.json").read_bytes(), checkpoint)
        self.assertEqual({str(p): p.stat().st_mtime_ns for p in out.rglob("stage.json")}, stages)

    def test_gap_cli_ont(self):
        self.gap(["ont"])

    def test_gap_cli_mixed(self):
        self.gap(["hifi", "ont"])

    def test_ctglinker_cli_preserves_complete_scaffold_and_components(self):
        from FamilyCtgLinker import load_ctg_result
        contigs = self.root / "contigs.fa"
        contigs.write_text(self.left.read_text() + self.right.read_text())
        out = self.root / "ctg"
        self.command("DEGAP.py", ["--mode", "ctglinker", "--out", str(out), "--ctgseq", str(contigs),
                                  "--hifi", str(self.root / "hifi.fa"), "--MaximumExtensionRound", "2", "-t", "2"])
        result = load_ctg_result(out)
        self.assertEqual(len(result["scaffolds"]), 1)
        self.assertEqual(str(SeqIO.read(out / "DG.scaffold.fa", "fasta").seq), self.truth)
        self.assertTrue((out / "ctglinker-family-report.html").is_file())
        self.assertTrue((out / "ctglinker_check.json").is_file())
        self.assertTrue((out / "processed_reads" / "hifi_reads.fa").is_file())

    def test_auto_cli_single_job_to_integrated_genome(self):
        genome = self.root / "genome.fa"
        genome.write_text(">chr1\n" + self.truth[:20000] + "N" * 20000 + self.truth[40000:] + "\n")
        out = self.root / "auto"
        self.command("AutoGapfiller.py", ["--genome", str(genome), "--hifi", str(self.root / "hifi.fa"),
                                         "--output", str(out), "-w", "1", "-t", "2",
                                         "--MaximumExtensionRound", "6", "--remove", "3"])
        manifests = list(out.rglob("family-result.json"))
        self.assertEqual(len(manifests), 1)
        final = list(out.rglob("genome.filled.fasta"))
        self.assertEqual(len(final), 1)
        self.assertEqual(str(SeqIO.read(final[0], "fasta").seq), self.truth)
        self.assertTrue((out / "05.visualization" / "Global.report.html").is_file())
        scripts = [path.read_text() for path in out.rglob("*.sh") if "DEGAP.py" in path.read_text()]
        self.assertTrue(scripts)
        self.assertTrue(all("--flag" not in script for script in scripts))
        stages = {str(p): p.stat().st_mtime_ns for p in out.rglob("stage.json")}
        self.command("AutoGapfiller.py", ["--genome", str(genome), "--hifi", str(self.root / "hifi.fa"),
                                         "--output", str(out), "-w", "1", "-t", "2",
                                         "--MaximumExtensionRound", "6", "--remove", "3"], log="resume.log")
        self.assertEqual(str(SeqIO.read(final[0], "fasta").seq), self.truth)
        self.assertEqual({str(p): p.stat().st_mtime_ns for p in out.rglob("stage.json")}, stages)

    @unittest.skipUnless(shutil.which("parallel"), "GNU parallel required for telomeric read preprocessing")
    def test_telseeker_ctg_cli_direct_result_and_verified_resume(self):
        import csv
        reads = self.root / "telomeric.fa"
        expected = self.truth[:30000] + "TTAGGG" * 120
        reads.write_text(">telo\n" + expected + "\n")
        out = self.root / "ctg-tel"
        args = ["--mode", "telseeker_ctg", "--out", str(out), "--ctg", str(self.left), "R",
                "--hifi", str(reads), "--motif", "TTAGGG", "--MaximumExtensionRound", "0", "-t", "2"]
        self.command("DEGAP.py", args)
        with (out / "result" / "status.tsv").open() as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["status"], "success")
        self.assertEqual(rows[0]["connection_type"], "direct")
        outputs = list(out.rglob("extended.ctg.fa"))
        self.assertEqual(len(outputs), 1)
        self.assertEqual(str(SeqIO.read(outputs[0], "fasta").seq), expected)
        receipts = {str(p): p.stat().st_mtime_ns for p in out.rglob("receipt.json")}
        self.assertTrue(receipts)
        self.command("DEGAP.py", args + ["--resume_auto"], log="resume.log")
        self.assertEqual({str(p): p.stat().st_mtime_ns for p in out.rglob("receipt.json")}, receipts)

    def test_telomere_public_extension_uses_one_family(self):
        from TelSeekerPart2 import TelomereExtender
        from FamilyTelomere import load_telomere_result
        work = self.root / "tel"
        work.mkdir()
        shutil.copyfile(self.left, work / "seed.fa")
        shutil.copyfile(self.right, work / "telo.reads.fa")
        processed = work / "processed_reads"
        processed.mkdir()
        reads = processed / "hifi_reads.fa"
        shutil.copyfile(self.root / "hifi.fa", reads)
        index = SeqIO.index_db(str(work / "hifi_reads.idx"), [str(reads)], "fasta")
        index.close()
        (work / "HiFi.reads.stat").write_text("SeedLength 15000\nMaxLength 20000\n")
        item = TelomereExtender("chr.R", work, work, 4, 41, 2, 500,
                                max_rounds=4, max_extension_length=100000, original_hifi=str(reads))
        result = item.run_extension()
        self.assertTrue(result["success"])
        self.assertEqual(str(SeqIO.read(result["connection_info"]["linker_file"], "fasta").seq), self.truth)
        manifest = load_telomere_result(work / "extension")
        self.assertEqual(manifest["result"], result)
        checkpoint = json.loads((work / "extension" / "families.json").read_text())
        self.assertEqual(len(checkpoint["families"]), 1)


if __name__ == "__main__":
    unittest.main()
