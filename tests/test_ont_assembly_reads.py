import gzip
import json
from pathlib import Path
import sqlite3
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from Bio import SeqIO
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from ONTAssemblyReads import restore_fastq, source_contract, prepare_ont_fastq
from GapFillerClass import mummer  # Existing tool modules initialize through this entrypoint.
from FindExtensionContigs import FindExtensionContigs


class ONTQualityTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.original = self.root / "reads.fastq"
        self.sequence = "AACCGTATGGCAT"
        self.quality = "!#$%&()*+,-./"
        self.original.write_text("@r1 original header\n" + self.sequence + "\n+\n" + self.quality + "\n")
        self.selected = self.root / "selected.fa"
        self.output = self.root / "out.fastq"
        self.cache = self.root / "cache"

    def restore(self, sequence=None, read_id="r1", sources=None):
        self.selected.write_text(">" + read_id + "\n" + (sequence or self.sequence) + "\n")
        restore_fastq(self.selected, source_contract(sources or self.original), self.output, self.cache)
        return SeqIO.read(self.output, "fastq")

    def scores(self, chars):
        return [ord(c)-33 for c in chars]

    def test_forward_and_reverse_preserve_true_qualities_and_cache(self):
        record = self.restore()
        self.assertEqual(record.letter_annotations["phred_quality"], self.scores(self.quality))
        database, = self.cache.glob("*.sqlite")
        stamp = database.stat().st_mtime_ns
        with patch("ONTAssemblyReads.FastqGeneralIterator", side_effect=AssertionError("reparsed original reads")):
            record = self.restore(str(Seq(self.sequence).reverse_complement()))
        self.assertEqual(record.letter_annotations["phred_quality"], self.scores(self.quality[::-1]))
        self.assertEqual(database.stat().st_mtime_ns, stamp)

    def test_unambiguous_forward_and_reverse_trims(self):
        record = self.restore(self.sequence[2:10])
        self.assertEqual(record.letter_annotations["phred_quality"], self.scores(self.quality[2:10]))
        record = self.restore(str(Seq(self.sequence[2:10]).reverse_complement()))
        self.assertEqual(record.letter_annotations["phred_quality"], self.scores(self.quality[2:10][::-1]))

    def test_gzip_multiple_sources_and_selected_order(self):
        second = self.root / "second.fq.gz"
        with gzip.open(second, "wt") as stream:
            stream.write("@r2\nGTCAAAA\n+\n1234567\n")
        self.selected.write_text(">r2\nGTCAAAA\n>r1\n"+self.sequence+"\n")
        restore_fastq(self.selected, source_contract([self.original, second]), self.output, self.cache)
        records=list(SeqIO.parse(self.output,"fastq"))
        self.assertEqual([r.id for r in records],["r2","r1"])
        self.assertEqual(records[0].letter_annotations["phred_quality"],self.scores("1234567"))

    def test_missing_id_fails_without_partial_fastq(self):
        with self.assertRaisesRegex(ValueError,"missing"):
            self.restore(read_id="absent")
        self.assertFalse(self.output.exists())

    def test_mismatched_sequence_fails(self):
        with self.assertRaisesRegex(ValueError,"mismatches"):
            self.restore("TTTTGGGGAAAA")
        self.assertFalse(self.output.exists())

    def test_ambiguous_trimming_fails(self):
        with self.assertRaisesRegex(ValueError,"ambiguous|Ambiguous"):
            self.restore("A")

    def test_duplicate_original_ids_fail(self):
        with self.assertRaisesRegex(ValueError,"Duplicate ONT"):
            self.restore(sources=[self.original,self.original])
        self.assertFalse(list(self.cache.glob("*.json")))

    def test_duplicate_selected_ids_fail(self):
        self.selected.write_text((">r1\n"+self.sequence+"\n")*2)
        with self.assertRaisesRegex(ValueError,"Duplicate"):
            restore_fastq(self.selected,source_contract(self.original),self.output,self.cache)
        self.assertFalse(self.output.exists())

    def test_fasta_source_is_rejected_without_fake_qualities(self):
        self.original.write_text(">r1\n"+self.sequence+"\n")
        with self.assertRaisesRegex(ValueError,"FASTQ.*not FASTA"):
            self.restore()
        self.assertFalse(list(self.cache.glob("*.json")))

    def test_absent_original_source_is_actionable(self):
        base=SimpleNamespace(out=str(self.root),original_reads_info=None)
        with self.assertRaisesRegex(ValueError,"genuine qualities"):
            prepare_ont_fastq(self.selected,base,self.output)

    def test_quality_only_input_change_rejects_pinned_contract(self):
        contract=source_contract(self.original)
        self.original.write_text(self.original.read_text().replace(self.quality,"I"*len(self.sequence)))
        with self.assertRaisesRegex(ValueError,"FASTQ changed"):
            restore_fastq(self.selected,contract,self.output,self.cache)

    def test_quality_only_change_changes_contract(self):
        before=source_contract(self.original)
        self.original.write_text(self.original.read_text().replace(self.quality,"I"*len(self.sequence)))
        self.assertNotEqual(before,source_contract(self.original))

    def test_corrupt_cache_is_not_silently_rebuilt(self):
        self.restore()
        database, = self.cache.glob("*.sqlite")
        database.write_bytes(b"broken")
        with self.assertRaisesRegex(ValueError,"cache changed"):
            self.restore()

    def test_malformed_quality_is_not_indexed(self):
        self.original.write_text("@r1\nACGT\n+\n!!!\n")
        with self.assertRaisesRegex(ValueError,"valid original FASTQ"):
            self.restore()
        self.assertFalse(list(self.cache.glob("*.json")))

    def test_ont_only_and_mixed_fallback_choose_ont_preparation(self):
        for data_type,status in (("ont",None),("mixed","ont_only")):
            with self.subTest(data_type=data_type):
                obj=FindExtensionContigs.__new__(FindExtensionContigs)
                obj.data_type=data_type
                obj.extensionReads=SimpleNamespace(mixed_mode_status=status)
                obj._prepareONTOnlyReads=lambda:"original-quality.fastq"
                obj._prepareHiFiOnlyReads=lambda:self.fail("ONT was treated as HiFi")
                self.assertEqual(obj._assembly_mode(),"ont")
                self.assertEqual(obj.prepareHifiasmReads(),"original-quality.fastq")

    def test_hifi_only_mixed_fallback_does_not_require_ont_qualities(self):
        obj=FindExtensionContigs.__new__(FindExtensionContigs)
        obj.data_type="mixed"
        obj.extensionReads=SimpleNamespace(mixed_mode_status="hifi_only")
        obj._prepareHiFiOnlyReads=lambda:"hifi.fa"
        self.assertEqual(obj._assembly_mode(),"hifi")
        self.assertEqual(obj.prepareHifiasmReads(),"hifi.fa")

if __name__=="__main__":
    unittest.main()
