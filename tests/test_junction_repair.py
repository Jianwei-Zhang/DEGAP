from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
import json
import random
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"bin"))
from ExtensionFamilies import AttemptResult, Connection, Family, FamilySearch, Version
from FamilyOverlap import MummerOverlap, joined_sequence, parse_coords, read_sequence as read_sequence_for_test
from JunctionRepair import (BoundaryVerifier, POLICY, load_plan, normalize_edit, operations,
                            original_edits, query_boundary, read_support, consensus_proposals)


class JunctionLogicTests(unittest.TestCase):
    def test_cs_boundary_projection_and_inverted_edits(self):
        r = {"qs":0,"ts":10,"ops":list(operations("=AC*gt=AA+CC=GG-TT=CA"))}
        self.assertEqual(query_boundary(r, 0), 10)
        self.assertEqual(query_boundary(r, 2), 12)
        self.assertEqual(original_edits(r), [
            {"start0":2,"end0":3,"original":"T","proposed":"G"},
            {"start0":5,"end0":7,"original":"CC","proposed":""},
            {"start0":9,"end0":9,"original":"","proposed":"TT"}])
        with self.assertRaises(ValueError):list(operations("=AC~gt20ag"))

    def test_indels_left_normalize(self):
        self.assertEqual(normalize_edit("CAAAAG",4,4,"","A"),(1,1,"","A"))
        self.assertEqual(normalize_edit("CAAAAG",3,4,"A",""),(1,2,"A",""))

    def records(self, alternate=5, original=0, other=0, single_strand=False):
        values=[]
        for i in range(alternate+original+other):
            if i < alternate:cs=":5900*ac:4099"
            elif i < alternate+original:cs=":10000"
            else:cs=":5900*ag:4099"
            values.append(("ont",{"name":str(i),"primary":True,"mapq":60,"identity":.99,
                                  "ts":0,"te":10000,"strand":"+" if single_strand or i%2 else "-",
                                  "ops":list(operations(cs))}))
        return values

    def support(self, records):
        return read_support("A"*10000,[{"start0":5900,"end0":5901,"original":"A","proposed":"C"}],
                            records,side="L",left_length=6000,gap_length=1000)[0]

    def test_repair_requires_distinct_spanning_reads_and_both_strands(self):
        self.assertEqual(self.support(self.records())["status"],"accepted")
        self.assertEqual(self.support(self.records(alternate=4))["status"],"retained")
        self.assertEqual(self.support(self.records(single_strand=True))["status"],"retained")
        self.assertEqual(self.support(self.records(original=2))["status"],"retained")
        records=self.records(alternate=4);records.append(records[0])
        self.assertEqual(self.support(records)["status"],"retained")

    def test_non_proposed_alleles_count_against_support(self):
        result=self.support(self.records(other=3))
        self.assertEqual((result["status"],result["other_reads"]),("retained",3))

    def test_repeat_mapping_and_missing_interior_anchor_not_support(self):
        records=self.records()
        for _,r in records:r["mapq"]=0
        self.assertEqual(self.support(records)["spanning_reads"],0)
        records=self.records()
        for _,r in records:r["ts"]=4000
        self.assertEqual(self.support(records)["spanning_reads"],0)

    def test_small_and_recent_input_hashes_do_not_trust_stat_only_cache(self):
        import hashlib
        from FamilyTools import content_hash
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"reads.fastq";p.write_bytes(b"AAAA")
            with patch("FamilyTools._hash_at_stamp",return_value="stale"):
                self.assertEqual(content_hash(p),hashlib.sha256(b"AAAA").hexdigest())
                p.write_bytes(b"CCCC")
                self.assertEqual(content_hash(p),hashlib.sha256(b"CCCC").hexdigest())
                p.write_bytes(b"G"*(1024*1024+1))
                self.assertEqual(content_hash(p),hashlib.sha256(p.read_bytes()).hexdigest())
    def test_read_consensus_keeps_boundary_substitutions_but_not_ambiguous_indels(self):
        ref="A"*5999+"C"+"G"*4000
        records=[("hifi",{"name":str(i),"primary":True,"mapq":60,"identity":.99,"ts":0,"te":10000,
                          "strand":"+" if i%2 else "-","ops":list(operations(":5999*ca:4000"))}) for i in range(5)]
        proposed=consensus_proposals(ref,records,side="L",left_length=6000,gap_length=1000)
        self.assertEqual(len(proposed),1)
        for _,record in records:record['ops']=list(operations(":5999-c:4000"))
        self.assertEqual(consensus_proposals(ref,records,side="L",left_length=6000,gap_length=1000),[])
    def test_anchor_retry_is_budgeted_and_resumable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);a=root/"a.fa";b=root/"b.fa"
            a.write_text(">a\n"+"A"*5000+"\n");b.write_text(">b\n"+"C"*5000+"\n")
            def extend(family,directory):
                if family.key=="L" and not family.anchor_backoff:return AttemptResult(retry_reason="terminal_anchor_backoff")
                return AttemptResult(stop_reason="no_extension_contigs")
            search=FamilySearch(root/"job",[("L","rightward",a),("R","leftward",b)],primary="L",context={},max_rounds=2)
            search.run(extend,lambda a,b:[])
            self.assertEqual((search.total_rounds,search.status),(2,"round_limit"))
            self.assertEqual(search.families["L"].anchor_backoff,2000)
            resumed=FamilySearch(root/"job",[("L","rightward",a),("R","leftward",b)],primary="L",context={},max_rounds=2)
            self.assertEqual(resumed.families["L"].anchor_backoff,2000)
            from FamilyResult import publish_search
            result=publish_search(search)
            self.assertEqual(result["anchor_retries"],[])
            self.assertEqual(result["anchor_backoff_requests"][0]["executed"],False)
            self.assertIn("no extra round charged",(search.directory/"process.log").read_text())
            actual=FamilySearch(root/"job3",[("L","rightward",a),("R","leftward",b)],primary="L",context={},max_rounds=3)
            actual.run(extend,lambda a,b:[])
            result=publish_search(actual)
            self.assertEqual(result["anchor_retries"][0]["attempt"],2)
            self.assertEqual(result["anchor_backoff_requests"][0]["executed"],True)
            self.assertEqual(result["total_rounds"],3)


@unittest.skipUnless(all(shutil.which(x) for x in ("nucmer","delta-filter","show-coords","minimap2")),"alignment tools required")
class JunctionToolTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);rng=random.Random(913)
        self.truth="".join(rng.choices("ACGT",k=150000))
    def write(self,name,text):
        path=self.root/(name+".fa");path.write_text(">"+name+"\n"+text+"\n");return path
    def test_large_overlap_uses_full_fallback(self):
        a=Version.capture("L",0,self.write("L",self.truth[:60000]))
        b=Version.capture("R",1,self.write("R",self.truth[30000:140000]))
        compare=MummerOverlap(self.root/"overlaps",500)
        hits=compare(a,b)
        self.assertTrue(hits)
        self.assertEqual(compare.diagnostic(a,b)["stages"][-1]["stage"],"full")
        self.assertEqual(joined_sequence(hits[0],a,b),self.truth[:140000])
    def test_contained_left_is_recognized_without_requiring_artificial_endpoints(self):
        a=Version.capture("L",0,self.write("L",self.truth[40000:80000]))
        b=Version.capture("R",1,self.write("R",self.truth[:140000]))
        compare=MummerOverlap(self.root/"overlaps",500)
        hits=compare(a,b)
        self.assertTrue(hits)
        self.assertEqual(joined_sequence(hits[0],a,b),self.truth[40000:140000])
    def test_already_full_short_inputs_can_use_containment_without_a_redundant_fallback(self):
        left=self.write("shortL",self.truth[3000:15000]);right=self.write("shortR",self.truth[17000:20000])
        a=Version.capture("L",0,left);b=Version.capture("R",1,self.write("shortExtendedR",self.truth[:20000]))
        compare=MummerOverlap(self.root/"overlaps",500)
        compare.configure_boundaries(SimpleNamespace(thread=1,reads_info={}),[("L","rightward",left),("R","leftward",right)])
        hit,=compare(a,b)
        self.assertEqual(joined_sequence(hit,a,b),self.truth[3000:20000])
        self.assertEqual(len(compare.diagnostic(a,b)["stages"]),1)
        self.assertEqual(compare.diagnostic(a,b)["stages"][0]["stage"],"fast")
    def test_original_overlap_is_explicit_deduplication_with_original_sources(self):
        left=self.write("overlapL",self.truth[:30000]);right=self.write("overlapR",self.truth[20000:50000])
        a=Version.capture("L",0,left);b=Version.capture("R",0,right)
        compare=MummerOverlap(self.root/"overlaps",500)
        compare.configure_boundaries(SimpleNamespace(thread=1,reads_info={}),[("L","rightward",left),("R","leftward",right)])
        hit,=compare(a,b);summary=load_plan(hit.evidence)["summary"]
        self.assertEqual(joined_sequence(hit,a,b),self.truth[:50000])
        self.assertEqual((summary["mode"],summary["gap_length"],summary["overlap_trim"]),("overlap",-10000,10000))
        self.assertEqual(summary["flank_length_delta"],0)
    def test_extension_can_pass_whole_opposite_flank_without_losing_connection(self):
        left=self.write("originL",self.truth[40000:80000]);right=self.write("originR",self.truth[100000:140000])
        a=Version.capture("L",0,left);b=Version.capture("R",1,self.write("overextendedR",self.truth[:140000]))
        compare=MummerOverlap(self.root/"overlaps",500)
        compare.configure_boundaries(SimpleNamespace(thread=1,reads_info={}),[("L","rightward",left),("R","leftward",right)])
        hit,=compare(a,b)
        self.assertEqual(joined_sequence(hit,a,b),self.truth[40000:140000])
        self.assertEqual(load_plan(hit.evidence)["summary"]["gap_length"],20000)
    def clipped_candidate(self, side, *, alt_count=10, ref_count=0, both_strands=True, clip=20):
        from Bio.Seq import Seq
        truth = self.truth[:80000]
        mutate = str.maketrans("ACGT", "CGTA")
        if side == "L":
            original_l = truth[:30000-clip]+truth[30000-clip:30000].translate(mutate)
            original_r = truth[50000:]
            extended_l, extended_r = truth[:65000], original_r
            replaced = original_l+truth[30000:]
            start = 20000
        else:
            original_l = truth[:30000]
            original_r = truth[50000:50000+clip].translate(mutate)+truth[50000+clip:]
            extended_l, extended_r = original_l, truth[15000:]
            replaced = truth[:50000]+original_r
            start = 40000
        label = f"{side}-{alt_count}-{ref_count}-{both_strands}-{clip}"
        left = self.write(label+"-L", original_l)
        right = self.write(label+"-R", original_r)
        a = Version.capture("L", int(side == "L"), self.write(label+"-extendedL", extended_l) if side == "L" else left)
        b = Version.capture("R", int(side == "R"), self.write(label+"-extendedR", extended_r) if side == "R" else right)
        reads = self.root/(label+"-reads.fa")
        with reads.open("w") as stream:
            for i in range(alt_count+ref_count):
                src = truth if i < alt_count else replaced
                seq = src[start+i*10:start+i*10+20000]
                if both_strands and i % 2:
                    seq = str(Seq(seq).reverse_complement())
                stream.write(f">read{i}\n{seq}\n")
        compare = MummerOverlap(self.root/("overlaps-"+label), 500)
        compare.configure_boundaries(SimpleNamespace(thread=1, reads_info={"working_hifi": str(reads)}),
                                     [("L", "rightward", left), ("R", "leftward", right)])
        return compare, a, b, truth

    def test_soft_clipped_ends_are_repaired_with_spanning_reads_on_either_side(self):
        for side in ("L", "R"):
            with self.subTest(side=side):
                compare,a,b,truth = self.clipped_candidate(side)
                hit, = compare(a,b)
                self.assertEqual(joined_sequence(hit,a,b), truth)
                summary = load_plan(hit.evidence)["summary"]
                edit, = [e for e in summary["repairs"] if e.get("candidate_source") == "terminal_projection"]
                self.assertEqual((edit["status"], edit["side"], edit["strands"]), ("accepted", side, ["+", "-"]))
                self.assertEqual((summary["gap_length"], summary["flank_length_delta"]), (20000, 0))
                self.assertGreater(summary["left_replace_bases" if side == "L" else "right_replace_bases"], 0)

    def test_clipped_boundaries_are_not_guessed_when_support_is_weak(self):
        for kwargs in ({"alt_count": 4}, {"both_strands": False}, {"ref_count": 3}):
            with self.subTest(**kwargs):
                compare,a,b,_ = self.clipped_candidate("L", **kwargs)
                self.assertEqual(compare(a,b), [])
                self.assertEqual(compare.diagnostic(a,b)["reason"],
                                 "insufficient_read_support_for_clipped_original_L_boundary")

    def test_clipping_past_protected_window_is_not_repaired(self):
        compare,a,b,_ = self.clipped_candidate("R", clip=2010)
        self.assertEqual(compare(a,b), [])
        self.assertEqual(compare.diagnostic(a,b)["reason"], "ambiguous_or_unmapped_original_R_boundary")

    def test_collinear_blocks_combine_without_counting_unaligned_gaps(self):
        a=Version("L",1,"unused","unused",50000);b=Version("R",1,"unused","unused",40000)
        path=self.root/"coords"
        def row(a,b,c,d,identity):return "\t".join(map(str,(a,b,c,d,b-a+1,d-c+1,identity,20500,20500,0,0,"right","left")))+"\n"
        path.write_text(row(138,2975,1773,4639,95.16)+row(3089,18840,4744,20500,99.54))
        hits=parse_coords(path,a,b,left_offset=29500,left_tail_length=20500,right_tail_length=20500,edge=500)
        self.assertEqual(len(hits),1)
        self.assertEqual(hits[0].overlap,2838+15752)
        self.assertEqual(hits[0].evidence["blocks"],2)
        path.write_text(row(138,2975,1773,4639,95.16)+row(4089,18840,5744,20500,99.54))
        self.assertEqual(parse_coords(path,a,b,left_offset=29500,left_tail_length=20500,right_tail_length=20500,edge=500),[])
    def setup_candidate(self, *, error_in_original=False, support=False, edit_kind="sub"):
        truth=self.truth[:25000]
        left=truth[:8000];right=truth[13000:]
        wrong="C" if left[-100]!="C" else "G"
        if error_in_original:
            if edit_kind == "insertion":left=left[:-100]+wrong+left[-100:]
            elif edit_kind == "deletion":left=left[:-100]+left[-99:]
            else:left=left[:-100]+wrong+left[-99:]
        else:
            wrong="C" if right[500]!="C" else "G"
        a=self.write("originalL",left);b=self.write("originalR",right)
        proposed=truth[:24000]
        if not error_in_original:proposed=proposed[:13500]+wrong+proposed[13501:]
        av=Version.capture("L",1,self.write("extendedL",proposed));bv=Version.capture("R",0,b)
        reads={}
        if support:
            seq=truth[:24000]
            from Bio.Seq import Seq
            reads_path=self.root/"reads.fa"
            reads_path.write_text("".join(f">r{i}\n{seq if i%2 else str(Seq(seq).reverse_complement())}\n" for i in range(8)))
            reads={"working_hifi":str(reads_path)}
        base=SimpleNamespace(thread=1,reads_info=reads)
        compare=MummerOverlap(self.root/"comparisons",500)
        compare.configure_boundaries(base,[("L","rightward",a),("R","leftward",b)])
        return compare,av,bv,left,right,truth
    def test_no_read_evidence_preserves_original_flanks_and_gap_length(self):
        compare,a,b,left,right,truth=self.setup_candidate()
        hit,=compare(a,b)
        output=joined_sequence(hit,a,b)
        self.assertEqual(output,left+truth[8000:13000]+right)
        plan=load_plan(hit.evidence)
        self.assertEqual((plan["summary"]["gap_length"],plan["summary"]["flank_length_delta"]),(5000,0))
        self.assertTrue(any(r["status"]=="retained" for r in plan["summary"]["repairs"]))
    def test_read_backed_terminal_error_is_corrected(self):
        compare,a,b,left,right,truth=self.setup_candidate(error_in_original=True,support=True)
        hit,=compare(a,b)
        self.assertEqual(joined_sequence(hit,a,b),truth)
        plan=load_plan(hit.evidence)
        self.assertEqual(plan["summary"]["flanks"]["L"]["accepted_variants"],1)
        self.assertEqual(plan["summary"]["gap_length"],5000)
        self.assertGreaterEqual(plan["summary"]["repairs"][0]["alt_reads"],5)
    def test_inherited_seed_error_is_discovered_from_reads_not_assembly_difference(self):
        compare,a,b,left,right,truth=self.setup_candidate(error_in_original=True,support=True)
        proposed=list(read_sequence_for_test(a.path));proposed[7900]=left[7900]
        a=Version.capture("L",1,self.write("inherited-error","".join(proposed)))
        hit,=compare(a,b)
        self.assertEqual(joined_sequence(hit,a,b),truth)
        edits=[r for r in load_plan(hit.evidence)["summary"]["repairs"] if r["status"]=="accepted"]
        self.assertEqual(edits[0]["candidate_source"],"read_consensus")
    def test_right_seed_error_in_both_overlap_copies_is_also_repaired(self):
        compare,a,b,left,right,truth=self.setup_candidate(support=True)
        bad_right=list(right);bad_right[500]="C" if right[500]!="C" else "G"
        origin=self.write("wrong-originR","".join(bad_right));b=Version.capture("R",0,origin)
        compare.configure_boundaries(compare.validator.base,[("L","rightward",compare.validator.origins["L"].path),("R","leftward",origin)])
        hit,=compare(a,b)
        self.assertEqual(joined_sequence(hit,a,b),truth)
        edits=[r for r in load_plan(hit.evidence)["summary"]["repairs"] if r["status"]=="accepted"]
        self.assertEqual(edits[0]["side"],"R")
        self.assertEqual(edits[0]["candidate_source"],"read_consensus")
    def test_read_backed_indel_is_not_misreported_as_gap_length_change(self):
        for kind, delta in (("insertion",-1),("deletion",1)):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as temp:
                old=self.root;self.root=Path(temp)
                try:
                    compare,a,b,left,right,truth=self.setup_candidate(error_in_original=True,support=True,edit_kind=kind)
                    hit,=compare(a,b)
                    self.assertEqual(joined_sequence(hit,a,b),truth)
                    summary=load_plan(hit.evidence)["summary"]
                    self.assertEqual((summary["gap_length"],summary["flank_length_delta"]),(5000,delta))
                finally:self.root=old
    def test_raw_reads_can_veto_an_incorrect_assembled_flank_change(self):
        compare,a,b,left,right,truth=self.setup_candidate(support=True)
        hit,=compare(a,b)
        self.assertEqual(joined_sequence(hit,a,b),truth)
        events=load_plan(hit.evidence)["summary"]["repairs"]
        self.assertEqual(events[0]["status"],"retained")
        self.assertGreaterEqual(events[0]["ref_reads"],5)
    def test_publication_exports_precise_integration_patch_and_portable_agp(self):
        from FamilyResult import publish_search, load_result
        compare,a,b,left,right,truth=self.setup_candidate(error_in_original=True,support=True)
        original_a=compare.validator.origins["L"]
        search=FamilySearch(self.root/"job",[("L","rightward",original_a.path),("R","leftward",b.path)],
                            primary="L",context=compare.contract,max_rounds=1)
        search.run(lambda family,directory:AttemptResult(a.path),compare)
        result=publish_search(search)
        self.assertEqual(result["junction"]["flanks"]["L"]["status"],"repaired")
        from Bio import SeqIO
        patch=str(SeqIO.read(search.directory/result["integration_sequence"],"fasta").seq)
        self.assertEqual(left[:-100]+patch+right,truth)
        with (search.directory/"family-components.fa").open() as stream:
            components={r.id:str(r.seq) for r in SeqIO.parse(stream,"fasta")}
        built=""
        for line in (search.directory/"process.agp").read_text().splitlines():
            row=line.split("\t");self.assertEqual(int(row[1])-1,len(built))
            built+=components[row[5]][int(row[6])-1:int(row[7])]
        self.assertEqual(built,truth)
        self.assertEqual(load_result(search.directory),result)
        self.assertIn("Gap-only length: 5000",(search.directory/"process.log").read_text())
    def test_composition_cache_is_bound_to_complete_version_content(self):
        compare,a,b,left,right,truth=self.setup_candidate()
        hit,=compare(a,b)
        old=joined_sequence(hit,a,b)
        alternate=list(self.truth[:24000]);alternate[10000]="C" if alternate[10000]!="C" else "G"
        # Keep the same overlap/flank change and version number: only a gap base differs.
        original=list(self.truth[:24000]);original[13500]=("C" if right[500]!="C" else "G")
        original[10000]=alternate[10000]
        other=Version.capture("L",a.number,self.write("different-version-content","".join(original)))
        newer,=compare(other,b)
        new=joined_sequence(newer,other,b)
        self.assertNotEqual(hit.evidence["plan_path"],newer.evidence["plan_path"])
        self.assertNotEqual(old,new)
        self.assertEqual(new[10000],alternate[10000])
    def test_cached_composition_tampering_is_not_rehashed_as_valid(self):
        compare,a,b,*_=self.setup_candidate();hit,=compare(a,b)
        path=Path(hit.evidence["plan_path"]);path.write_text(path.read_text()+" ")
        with self.assertRaisesRegex(ValueError,"composition changed"):compare(a,b)
    def test_composition_tampering_is_rejected(self):
        compare,a,b,*_=self.setup_candidate();hit,=compare(a,b)
        path=Path(hit.evidence["plan_path"]);path.write_text(path.read_text()+" ")
        with self.assertRaisesRegex(ValueError,"receipt changed"):joined_sequence(hit,a,b)


if __name__=="__main__":unittest.main()
