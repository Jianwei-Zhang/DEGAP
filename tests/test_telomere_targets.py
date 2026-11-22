"""Target retrieval, single-ended decision and integrity; real batch fixtures."""
import csv
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"bin"))
from TelomereTargets import (POLICY, TargetIndex, TelomereTargets, seed_code,
                             query_seeds, prefix_projection, publish_check)
from Bio import SeqIO


def dna(n, seed=1):
    rng=random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(n))


def write_fasta(path, records):
    path=Path(path)
    path.write_text(''.join(f'>{name}\n{seq}\n' for name,seq in records))
    return path


class IndexTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)

    def index(self,records):
        pool=write_fasta(self.root/'pool.fa',records)
        item=TargetIndex(pool,self.root/'shared')
        self.addCleanup(item.close)
        return item

    def test_tandem_repeats_are_not_query_anchors(self):
        for pattern in ('A','AG','TTAGGG','CCCTAA','ACGTTAGCCAAT'):
            word=(pattern*30)[:21]
            self.assertIsNone(seed_code(word),pattern)
        self.assertEqual(query_seeds('TTAGGG'*200),{})
        self.assertTrue(query_seeds(dna(1000)))

    def test_complete_sampled_seed_retrieval_across_all_twenty_phases(self):
        shared=dna(40,20)
        item=self.index([(f'r{i}',dna(i,300+i)+shared+dna(31,200+i)) for i in range(20)])
        self.assertEqual(len(item.candidates(query_seeds(shared))),20)

    def test_identical_records_deduplicated_but_read_ids_not_fabricated(self):
        seq=dna(700)
        item=self.index([('one',seq),('one',seq),('two',seq),('one',dna(700,2))])
        self.assertEqual(item.manifest['raw_records'],4)
        self.assertEqual(item.manifest['indexed_records'],3)
        self.assertEqual(item.manifest['distinct_read_ids'],2)
        self.assertEqual(item.manifest['duplicate_records'],1)

    def test_index_reused_without_rebuild_and_unrelated_sources_share_by_content(self):
        pool=write_fasta(self.root/'a.fa',[('r',dna(1000))])
        other=self.root/'b.fa';shutil.copyfile(pool,other)
        first=TargetIndex(pool,self.root/'shared');first.close()
        with patch.object(TargetIndex,'_build',side_effect=AssertionError('rebuild')):
            second=TargetIndex(other,self.root/'shared')
            second.close()
        self.assertEqual(first.directory,second.directory)

    def test_corrupt_index_is_not_rebuilt_or_accepted(self):
        item=self.index([('r',dna(700))])
        path=item.directory/'targets.fa';path.write_text(path.read_text().replace('A','C',1))
        with self.assertRaisesRegex(ValueError,'Changed completed'):
            TargetIndex(item.pool,self.root/'shared')

    def test_projection_preserves_chromosome_and_handles_insertions(self):
        ref='ACGTTGCACT'
        self.assertEqual(prefix_projection(ref,ref+'GG'),len(ref))
        self.assertEqual(prefix_projection(ref,'ACGTTAGCACTGG'),len(ref)+1)
        self.assertIsNone(prefix_projection('A','CAA'))
        self.assertIsNone(prefix_projection(ref,'T'*20))

    def test_invalid_sequences_fail_before_index_publication(self):
        with self.assertRaisesRegex(ValueError,'Invalid telomere'):
            self.index([('r','ACGT?')])
        self.assertFalse(list((self.root/'shared').glob('*/index.json')))


    def test_workers_share_one_completed_index_construction(self):
        from concurrent.futures import ThreadPoolExecutor
        pool=write_fasta(self.root/'shared-pool.fa',[('a',dna(1500)),('b',dna(1500,2))])
        def build(_):
            item=TargetIndex(pool,self.root/'shared')
            try:
                return item.directory,item.manifest['indexed_records']
            finally:
                item.close()
        with ThreadPoolExecutor(max_workers=2) as workers:
            results=list(workers.map(build,range(2)))
        self.assertEqual(results[0],results[1])
        self.assertEqual(len(list((self.root/'shared').glob('*/index.json'))),1)


    def test_ambiguous_iupac_bases_are_retained_but_never_anchor_seeds(self):
        sequence=dna(700)+'RYSWKMBDHVN'+dna(700,2)
        item=self.index([('ambiguous',sequence)])
        self.assertEqual(item.record(1)[1],sequence)
        self.assertEqual(item.candidates(query_seeds(sequence)),[1])
        self.assertIsNone(seed_code('ACGTRYSWKMBDHVNACGTAC'))

    def test_bulk_postings_are_exact_definition_not_approximate_prefilter(self):
        records=[('a',dna(1400)),('b',dna(1700,3)),('c','TTAGGG'*300)]
        item=self.index(records)
        expected=[]
        for number,(_,seq) in enumerate(records,1):
            codes={seed_code(seq[pos:pos+21]) for pos in range(0,len(seq)-20,20)}
            expected.extend((code,number) for code in codes if code is not None)
        actual=[tuple(row) for row in item.db.execute('SELECT code,number FROM seeds ORDER BY code,number')]
        self.assertEqual(actual,sorted(expected))
        self.assertEqual(item.manifest['builder'],'sequential-postings-then-covering-index')


TOOLS=all(shutil.which(t) for t in ('nucmer','delta-filter','show-coords'))


@unittest.skipUnless(TOOLS,'MUMmer is required for real batched telomere fixtures')
class BatchSearchTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.chrom=dna(18000,81)
        self.tail=dna(400,82)+'TTAGGG'*60
        self.query=write_fasta(self.root/'chr.fa',[('chr',self.chrom)])

    def targets(self,records,*,name='pool',end='R'):
        pool=write_fasta(self.root/(name+'.fa'),records)
        item=TelomereTargets(pool,self.root/'shared',self.root/(name+'-alignments'))
        self.addCleanup(item.close)
        return item

    def test_best_read_not_first_hit_and_single_batched_tool_chain(self):
        bad=list(self.chrom[-12000:])
        for i in range(100,len(bad)-100,91):
            bad[i]='A' if bad[i]!='A' else 'C'
        records=[('bad-first',''.join(bad)+self.tail),('z-best',self.chrom[-12000:]+self.tail),
                 ('y-support',self.chrom[-11000:]+self.tail)]
        item=self.targets(records)
        real=subprocess.run
        with patch('FamilyOverlap.subprocess.run',wraps=real) as calls:
            root,result=item.search(self.query,'R')
        self.assertEqual(len(calls.call_args_list),3)
        self.assertEqual(result['read_id'],'z-best')
        self.assertEqual(str(SeqIO.read(result['linker_file'],'fasta').seq),self.chrom+self.tail)
        self.assertEqual(result['support_count_lower_bound'],2)
        d=json.loads((root/'connection-diagnostics.json').read_text())
        self.assertEqual(d['refined_records'],3)
        self.assertIsNone(d['candidate_cap'])
        other=self.targets(list(reversed(records)),name='reordered')
        _,again=other.search(self.query,'R')
        self.assertEqual(again['read_id'],result['read_id'])
        self.assertEqual(str(SeqIO.read(again['linker_file'],'fasta').seq),self.chrom+self.tail)

    def test_later_supported_candidate_beats_unsupported_highest_rank(self):
        other=dna(400,902)+'TTAGGG'*60
        item=self.targets([('unsupported',self.chrom[-13000:]+other),
                           ('supported',self.chrom[-12000:]+self.tail),
                           ('witness',self.chrom[-11000:]+self.tail)])
        root,result=item.search(self.query,'R')
        self.assertEqual(result['read_id'],'supported')
        support=json.loads((root/'connection-support.json').read_text())['attempted']
        self.assertEqual([r['status'] for r in support],['insufficient_support','supported'])

    def test_duplicate_read_id_is_not_independent_support(self):
        seq=self.chrom[-12000:]+self.tail
        item=self.targets([('same',seq),('same',seq)])
        root,result=item.search(self.query,'R')
        self.assertIsNone(result)
        d=json.loads((root/'connection-diagnostics.json').read_text())
        self.assertEqual(d['reason'],'insufficient_read_support')
        self.assertEqual(d['distinct_read_ids'],1)

    def test_left_end_preserves_entire_chromosome(self):
        prefix='CCCTAA'*60+dna(400,90)
        item=self.targets([('best',prefix+self.chrom[:12000]),('support',prefix+self.chrom[:11000])])
        _,result=item.search(self.query,'L')
        self.assertIsNotNone(result)
        self.assertEqual(str(SeqIO.read(result['linker_file'],'fasta').seq),prefix+self.chrom)

    def test_large_unrelated_pool_not_scanned_one_read_at_a_time(self):
        records=[(f'decoy{i}',dna(300,1000+i)+'TTAGGG'*40) for i in range(1500)]
        records += [('best',self.chrom[-12000:]+self.tail),('support',self.chrom[-11000:]+self.tail)]
        item=self.targets(records)
        root,result=item.search(self.query,'R')
        d=json.loads((root/'connection-diagnostics.json').read_text())
        self.assertEqual(d['pool_records'],1502)
        self.assertEqual(d['retrieved_records'],2)
        self.assertEqual(result['read_id'],'best')
        self.assertFalse(list(root.glob('temp_*.fa')))

    def test_repeat_only_query_never_runs_alignment_or_assigns_locus(self):
        write_fasta(self.query,[('repeat','TTAGGG'*3000)])
        item=self.targets([('a','TTAGGG'*3000),('b','TTAGGG'*3000)])
        with patch('FamilyOverlap.subprocess.run',side_effect=AssertionError('repeat-only mapping')):
            root,result=item.search(self.query,'R')
        self.assertIsNone(result)
        self.assertEqual(json.loads((root/'connection-diagnostics.json').read_text())['reason'],
                         'no_sampled_nonrepeat_anchor')

    def test_search_and_published_evidence_resume_without_tools(self):
        item=self.targets([('best',self.chrom[-12000:]+self.tail),('support',self.chrom[-11000:]+self.tail)])
        dest=self.root/'direct.check'
        result=publish_check(item,self.query,'R',dest,0)
        stamps={p.name:p.stat().st_mtime_ns for p in dest.iterdir()}
        with patch('FamilyOverlap.subprocess.run',side_effect=AssertionError('rerun')):
            again=publish_check(item,self.query,'R',dest,0)
        self.assertEqual(result,again)
        self.assertEqual(stamps,{p.name:p.stat().st_mtime_ns for p in dest.iterdir()})
        (dest/'candidates.tsv').write_text('tampered')
        with self.assertRaisesRegex(ValueError,'Changed published'):
            publish_check(item,self.query,'R',dest,0)

    def test_tools_and_completed_cache_corruption_are_errors(self):
        item=self.targets([('best',self.chrom[-12000:]+self.tail),('support',self.chrom[-11000:]+self.tail)])
        root,_=item.search(self.query,'R')
        (root/'batch-000001'/'alignment.coords').write_text('corrupt')
        with self.assertRaisesRegex(ValueError,'Changed completed'):
            item.search(self.query,'R')

    def test_pool_mutation_during_open_index_is_rejected(self):
        item=self.targets([('best',self.chrom[-12000:]+self.tail),('support',self.chrom[-11000:]+self.tail)])
        item.index.pool.write_text(item.index.pool.read_text().replace('A','C',1))
        with self.assertRaisesRegex(ValueError,'pool changed'):
            item.search(self.query,'R')

    def test_negative_search_is_cached_with_reason(self):
        item=self.targets([('one',self.chrom[-12000:]+self.tail)])
        root,result=item.search(self.query,'R');self.assertIsNone(result)
        with patch('FamilyOverlap.subprocess.run',side_effect=AssertionError('rerun')):
            self.assertIsNone(item.search(self.query,'R')[1])
        self.assertEqual(json.loads((root/'connection-diagnostics.json').read_text())['reason'],
                         'insufficient_read_support')


    def test_round_zero_public_path_and_extension_skip_are_integrity_checked(self):
        from TelSeekerPart2 import check_single_chr_end, extend_single_chr_end
        work=self.root/'chr.R';work.mkdir()
        write_fasta(work/'seed.fa',[('chr',self.chrom)])
        write_fasta(work/'telo.reads.fa',[('best',self.chrom[-12000:]+self.tail),
                                        ('support',self.chrom[-11000:]+self.tail)])
        direct=check_single_chr_end('chr.R',work)
        self.assertTrue(direct['success'],direct)
        skipped=extend_single_chr_end('chr.R',work,self.root,10,41,1,500,max_rounds=0)
        self.assertTrue(skipped['success'],skipped)
        self.assertTrue(skipped['skipped'])
        (work/'direct.check'/'connection-support.json').write_text('{}')
        broken=extend_single_chr_end('chr.R',work,self.root,10,41,1,500,max_rounds=0)
        self.assertTrue(broken.get('error'))
        self.assertFalse(broken['success'])

    def test_single_family_checkpoint_records_new_target_evidence_and_reuses_it(self):
        from TelSeekerPart2 import TelomereExtender
        from FamilyTelomere import run_telomere, load_telomere_result
        from ExtensionFamilies import AttemptResult
        work=self.root/'chr.R';work.mkdir()
        write_fasta(work/'seed.fa',[('chr',self.chrom[:9000])])
        write_fasta(work/'telo.reads.fa',[('best',self.chrom[-12000:]+self.tail),
                                        ('support',self.chrom[-11000:]+self.tail)])
        (work/'processed_reads').mkdir()
        reads=write_fasta(work/'processed_reads'/'hifi_reads.fa',[('r',self.chrom)])
        index=SeqIO.index_db(str(work/'hifi_reads.idx'),[str(reads)],'fasta');index.close()
        (work/'HiFi.reads.stat').write_text('SeedLength 9000\nMaxLength 18000\n')
        item=TelomereExtender('chr.R',work,work,10,41,1,500,max_rounds=2,
                             max_extension_length=100000,original_hifi=str(reads))
        self.addCleanup(item.extension_checker.close)
        self.addCleanup(item.base.readsDict.close)
        middle=write_fasta(self.root/'middle.fa',[('chr',self.chrom[:14000])])
        versions=iter([middle,self.query])
        result=run_telomere(item,extend=lambda family,directory:AttemptResult(str(next(versions))))
        self.assertTrue(result['success'])
        self.assertEqual(result['rounds'],2)
        self.assertEqual(result['total_extension_length'],9000)
        manifest=load_telomere_result(work/'extension')
        self.assertEqual(manifest['rounds'][0]['target_search']['reason'],'no_eligible_alignment')
        self.assertEqual(manifest['rounds'][1]['target_search']['reason'],'connected')
        self.assertEqual(len(json.loads((work/'extension'/'families.json').read_text())['families']),1)
        with patch('FamilyOverlap.subprocess.run',side_effect=AssertionError('rerun')):
            resumed=run_telomere(item,extend=lambda *args: (_ for _ in ()).throw(AssertionError('extension rerun')))
        self.assertEqual(result,resumed)


    def test_failed_tool_is_error_not_cached_biological_no_hit(self):
        item=self.targets([('best',self.chrom[-12000:]+self.tail),('support',self.chrom[-11000:]+self.tail)])
        with patch('FamilyOverlap.subprocess.run',side_effect=subprocess.CalledProcessError(2,['nucmer'])):
            with self.assertRaises(subprocess.CalledProcessError):
                item.search(self.query,'R')
        self.assertFalse(list(item.directory.glob('*/search.json')))

    def test_malformed_tool_coordinates_fail_before_publishing_decision(self):
        item=self.targets([('best',self.chrom[-12000:]+self.tail),('support',self.chrom[-11000:]+self.tail)])
        original=item.aligner._stage
        def malformed(directory,name,command,inputs,output,**kwargs):
            original(directory,name,command,inputs,output,**kwargs)
            if name=='coords':
                Path(output).write_text('malformed\n')
        with patch.object(item.aligner,'_stage',side_effect=malformed):
            with self.assertRaisesRegex(ValueError,'Malformed'):
                item.search(self.query,'R')
        self.assertFalse(list(item.directory.glob('*/search.json')))

    def test_legacy_direct_log_is_not_a_completed_receipt(self):
        item=self.targets([('best',self.chrom[-12000:]+self.tail),('support',self.chrom[-11000:]+self.tail)])
        dest=self.root/'direct.check';dest.mkdir()
        (dest/'direct.check.log').write_text('Status: SUCCESS\n')
        with self.assertRaisesRegex(ValueError,'Legacy'):
            publish_check(item,self.query,'R',dest,0)

    def test_two_orientations_of_same_id_are_not_two_supporting_reads(self):
        seq=self.chrom[-12000:]+self.tail
        item=self.targets([('same',seq),('same',self.chrom[-11000:]+self.tail)])
        root,result=item.search(self.query,'R')
        self.assertIsNone(result)
        self.assertEqual(json.loads((root/'connection-diagnostics.json').read_text())['distinct_read_ids'],1)

    def test_direct_manager_does_not_swallow_operational_failure(self):
        from TelSeekerPart2 import DirectConnectionChecker
        manager=DirectConnectionChecker(self.root,['chr.R'],{'chr.R':self.root},max_workers=1)
        with patch('TelSeekerPart2.check_single_chr_end',return_value={'chr_end':'chr.R','success':False,'error':'tool failed'}):
            with self.assertRaisesRegex(RuntimeError,'tool failed'):
                manager.run_checks()


    def test_batches_are_not_candidate_cap_or_per_batch_first_hit(self):
        with patch.dict(POLICY,{'batch_records':1}):
            item=self.targets([('short',self.chrom[-10500:]+self.tail),
                               ('support',self.chrom[-11000:]+self.tail),
                               ('best-last-batch',self.chrom[-12000:]+self.tail)])
            root,result=item.search(self.query,'R')
        d=json.loads((root/'connection-diagnostics.json').read_text())
        self.assertEqual(d['alignment_batches'],3)
        self.assertEqual(d['refined_records'],3)
        self.assertEqual(result['read_id'],'best-last-batch')
        self.assertEqual(len(list(root.glob('batch-*'))),3)


if __name__=='__main__':
    unittest.main()
