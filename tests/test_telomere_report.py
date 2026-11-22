"""Focused unittest tests for telomere target-search reporting and TelSeekerVisualizer integration."""

import ast
import html as html_mod
import json
from pathlib import Path
import re
import sys
import tempfile
import unittest

BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from ExtensionFamilies import fingerprint
from FamilyTools import content_hash
from TelomereReport import (
    load_check_evidence,
    render_target_search_section,
    safe_relative_href,
    summarize_chr_end_evidence,
)
from TelSeekerVisualizer import (
    HAS_BIOPYTHON,
    Part2HTMLReportGenerator,
    TelSeekerPart2LogParser,
    generate_all_task_reports,
    generate_html_visualization,
    parse_part2_results,
)


def _sample_diagnostics(**overrides):
    diag = {
        "schema": "telomere-target-search-v1",
        "strategy": "shared-exact-seed-index-batched-mummer",
        "policy": {
            "version": "telomere-targets-v1",
            "k": 21,
            "stride": 20,
            "anchor_span": 500,
            "anchor_seeds": 3,
            "min_support": 2,
            "support_bases": 200,
            "support_identity": 80.0,
            "projection_limit": 1000,
            "ranking": "identity-overlap-overhang-readid-sequence-v1",
        },
        "end_type": "L",
        "pool_sha256": "a" * 64,
        "index_identity": "b" * 64,
        "pool_records": 12450,
        "indexed_records": 12400,
        "distinct_read_ids": 12320,
        "query_bases": 45000,
        "query_anchor_keys": 4120,
        "retrieved_records": 38,
        "refined_records": 38,
        "alignment_blocks": 52,
        "alignment_batches": 2,
        "eligible_records": 4,
        "rejection_counts": {
            "short_alignment": 21,
            "nonterminal_alignment": 11,
            "insufficient_nonrepeat_anchor": 4,
        },
        "candidate_cap": None,
        "global_grouping": False,
        "target_extension": False,
        "selection_scope": "all_index_retrieved_records_with_qualifying_anchor_and_read_support",
        "reason": "connected",
        "winner_read_id": "tel_read_001",
        "support_checks": 1,
        "supporting_read_ids": ["tel_read_001", "tel_read_002"],
        "support_count_lower_bound": 2,
        "selected_tail_bases": 3420,
        "junction_bases_tested": 200,
        "support_covers": "common_nonrepeat_anchor_and_junction_proximal_sequence_not_full_telomere_length",
        "chromosome_sequence_retained": True,
    }
    diag.update(overrides)
    return diag


def _write_check_dir(
    check_dir: Path,
    diag: dict,
    *,
    attempted=None,
    candidates_rows=None,
    connection=None,
):
    check_dir.mkdir(parents=True, exist_ok=True)
    diag_path = check_dir / "connection-diagnostics.json"
    diag_path.write_text(json.dumps(diag, indent=2), encoding="utf-8")

    support_path = check_dir / "connection-support.json"
    support_payload = {
        "schema": "telomere-support-v1",
        "attempted": attempted if attempted is not None else [],
    }
    support_path.write_text(json.dumps(support_payload, indent=2), encoding="utf-8")

    cand_path = check_dir / "candidates.tsv"
    lines = ["read_id\trecord\tidentity\teffective_overlap\tendpoint_overhang\tstatus"]
    for row in (candidates_rows or []):
        lines.append("\t".join(str(x) for x in row))
    cand_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    receipt_path = check_dir / "telomere-check.json"
    files = {
        name: content_hash(check_dir / name)
        for name in ("connection-diagnostics.json", "connection-support.json", "candidates.tsv")
    }
    receipt = {
        "schema": "telomere-check-v1",
        "identity": "check-id-1",
        "connection": connection,
        "files": files,
    }
    receipt["payload_sha256"] = fingerprint(receipt)
    receipt_path.write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    return diag_path


class TestTelomereReport(unittest.TestCase):
    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.tmp_path = Path(self._tmpdir.name)

    def test_schema_values_and_bilingual_scientific_accuracy(self):
        tmp_path = self.tmp_path
        report_dir = tmp_path / "visual.report"
        report_dir.mkdir()
        check_dir = tmp_path / "part2.chr.end.job" / "Chr1.L" / "direct.check"
        diag = _sample_diagnostics()
        _write_check_dir(
            check_dir,
            diag,
            attempted=[
                {
                    "read_id": "tel_read_001",
                    "number": 10,
                    "status": "supported",
                    "supporting_read_ids": ["tel_read_001", "tel_read_002"],
                    "support_is_lower_bound": True,
                }
            ],
            candidates_rows=[
                ("tel_read_001", 10, "99.4", 14200, 12, "supported"),
                ("tel_read_002", 11, "99.1", 13900, 20, "ranked_below_winner"),
                ("tel_read_bad", 12, "", "", "", "alignment_rejected"),
            ],
        )
        ev = load_check_evidence(check_dir)
        data = {
            "direct_check": {
                "checked": True,
                "success": True,
                "target_search": ev["diagnostics"],
                "support_evidence": ev["support"],
                "candidates_summary": ev["candidates_summary"],
                "evidence_files": ev["evidence_files"],
            },
            "rounds": [],
            "final_info": {"total_extension_length": 12000},
        }

        html = render_target_search_section(data, report_dir)

        # Exact schema counts, batching, and separated length fields
        for token in (
            "shared-exact-seed-index-batched-mummer",
            "12,450",
            "12,400",
            "12,320",
            "4,120",
            "38",
            "52",
            "4",
            "500",
            "200 bp",
            "3,420 bp",
            "12,000 bp",
            "tel_read_001",
            "tel_read_002",
            "common_nonrepeat_anchor_and_junction_proximal_sequence_not_full_telomere_length",
        ):
            assert token in html

        # Bounded-batch wording (no promise of exactly one MUMmer pipeline for huge pools)
        assert "256 records / 8 MiB" in html
        assert "256 条记录 / 8 MiB" in html
        assert "in one batched MUMmer comparison" not in html
        assert "执行一次批量 MUMmer 验证" not in html
        assert "Batches / Blocks / Eligible" in html and "批次 / 比对块 / 合格" in html
        assert "Alignment Batches" in html and "比对批次数" in html

        # Separated length labels & final-join retention scope
        assert "Chr-Side Assembly Growth (Excl. Tel-Read Tail)" in html
        assert "染色体侧组装延伸（不含端粒 Read 尾）" in html
        assert "Selected Tel-Read Tail" in html and "入选端粒 Read 尾长" in html
        assert "Outward Junction Span Tested" in html and "实测向外接合区跨度" in html
        assert "Current Chr-Side Seq Retained (Final Join)" in html
        assert "最终拼接保留当前染色体侧序列" in html
        assert "FINAL JOIN" in html
        assert "does NOT guarantee that earlier extension/assembly rounds preserved the initial flank" in html
        assert "并不保证此前各轮延伸/组装完整保留了初始侧翼序列" in html

        # Invariants & bilingual scientific scope
        assert "None (no cap)" in html
        assert "无上限" in html
        assert "Global Target Grouping" in html and "全局目标分组" in html
        assert "Telomeric-Read Extension" in html and "端粒 Read 延伸" in html
        assert "shared exact sampled-seed index" in html
        assert "bounded-batch MUMmer" in html
        assert "not first hit" in html
        assert "NOT full-read polishing" in html
        assert "NOT telomere-repeat-length consensus" in html
        assert "NOT guaranteed full-genome specificity" in html
        assert "并非全 read 抛光" in html
        assert "并非端粒重复长度共识" in html
        assert "不保证全基因组唯一性" in html

        # Rejection counts rendered cleanly
        assert "Short alignment" in html and "21" in html
        assert "Non-terminal alignment" in html and "11" in html
        assert "Insufficient non-repeat anchor" in html and "4" in html

        # Safe relative links to all 4 evidence files
        for fname in (
            "../part2.chr.end.job/Chr1.L/direct.check/connection-diagnostics.json",
            "../part2.chr.end.job/Chr1.L/direct.check/connection-support.json",
            "../part2.chr.end.job/Chr1.L/direct.check/candidates.tsv",
            "../part2.chr.end.job/Chr1.L/direct.check/telomere-check.json",
        ):
            assert f'href="{fname}"' in html

    def test_negative_reasons_and_rejected_candidates(self):
        tmp_path = self.tmp_path
        report_dir = tmp_path / "visual.report"
        report_dir.mkdir()

        r1_dir = tmp_path / "part2.chr.end.job" / "Chr2.R" / "direct.check"
        d1 = _sample_diagnostics(
            end_type="R",
            reason="no_sampled_nonrepeat_anchor",
            retrieved_records=0,
            refined_records=0,
            alignment_blocks=0,
            alignment_batches=0,
            eligible_records=0,
            rejection_counts={},
            winner_read_id=None,
            support_checks=0,
            supporting_read_ids=[],
            support_count_lower_bound=0,
            selected_tail_bases=None,
            junction_bases_tested=None,
        )
        _write_check_dir(r1_dir, d1)
        ev1 = load_check_evidence(r1_dir)

        r2_dir = tmp_path / "part2.chr.end.job" / "Chr2.R" / "extension" / "attempts" / "0" / "1" / "telomere" / "connection_check"
        d2 = _sample_diagnostics(
            end_type="R",
            reason="no_eligible_alignment",
            retrieved_records=15,
            refined_records=15,
            alignment_blocks=19,
            alignment_batches=1,
            eligible_records=0,
            rejection_counts={
                "low_identity": 7,
                "ambiguous_boundary_projection": 5,
                "inverted_alignment": 3,
            },
            winner_read_id=None,
            support_checks=0,
            supporting_read_ids=[],
            support_count_lower_bound=0,
            selected_tail_bases=None,
            junction_bases_tested=None,
        )
        _write_check_dir(r2_dir, d2)
        ev2 = load_check_evidence(r2_dir)

        r3_dir = tmp_path / "part2.chr.end.job" / "Chr2.R" / "extension" / "attempts" / "0" / "2" / "telomere" / "connection_check"
        d3 = _sample_diagnostics(
            end_type="R",
            reason="insufficient_read_support",
            retrieved_records=12,
            refined_records=12,
            alignment_blocks=14,
            alignment_batches=1,
            eligible_records=2,
            rejection_counts={"short_alignment": 10},
            winner_read_id=None,
            support_checks=2,
            supporting_read_ids=[],
            support_count_lower_bound=0,
            selected_tail_bases=None,
            junction_bases_tested=200,
        )
        _write_check_dir(
            r3_dir,
            d3,
            attempted=[
                {
                    "read_id": "lone_candidate_A",
                    "number": 5,
                    "status": "insufficient_support",
                    "supporting_read_ids": ["lone_candidate_A"],
                    "support_is_lower_bound": True,
                },
                {
                    "read_id": "lone_candidate_B",
                    "number": 9,
                    "status": "insufficient_support",
                    "supporting_read_ids": ["lone_candidate_B"],
                    "support_is_lower_bound": True,
                },
            ],
            candidates_rows=[
                ("lone_candidate_A", 5, "98.1", 11200, 40, "insufficient_support"),
                ("lone_candidate_B", 9, "96.5", 10400, 80, "insufficient_support"),
            ],
        )
        ev3 = load_check_evidence(r3_dir)

        data = {
            "direct_check": {
                "checked": True,
                "success": False,
                "target_search": ev1["diagnostics"],
                "support_evidence": ev1["support"],
                "candidates_summary": ev1["candidates_summary"],
                "evidence_files": ev1["evidence_files"],
            },
            "rounds": [
                {
                    "round": 1,
                    "extension_status": "extended",
                    "connection_found": False,
                    "target_search": ev2["diagnostics"],
                    "support_evidence": ev2["support"],
                    "candidates_summary": ev2["candidates_summary"],
                    "evidence_files": ev2["evidence_files"],
                },
                {
                    "round": 2,
                    "extension_status": "extended",
                    "connection_found": False,
                    "target_search": ev3["diagnostics"],
                    "support_evidence": ev3["support"],
                    "candidates_summary": ev3["candidates_summary"],
                    "evidence_files": ev3["evidence_files"],
                },
                {
                    "round": 3,
                    "extension_status": "no_extension",
                    "connection_found": False,
                },
            ],
        }

        html = render_target_search_section(data, report_dir)

        # All three negative reasons + bilingual labels
        assert "no_sampled_nonrepeat_anchor" in html
        assert "No sampled non-repeat anchor" in html and "无采样非重复锚点命中" in html
        assert "no_eligible_alignment" in html
        assert "No eligible alignment" in html and "无合格比对" in html
        assert "insufficient_read_support" in html
        assert "Insufficient read support" in html and "支持 read 不足" in html

        # Failed round and rejected candidates preserved
        assert "Not checked (no extension)" in html and "未检查（无延伸）" in html
        assert "lone_candidate_A" in html
        assert "lone_candidate_B" in html
        assert "Ambiguous boundary projection" in html
        assert "Inverted alignment" in html
        assert 'class="row-selected"' not in html

    def test_lower_bound_wording_never_exact_total(self):
        tmp_path = self.tmp_path
        report_dir = tmp_path / "visual.report"
        report_dir.mkdir()
        diag = _sample_diagnostics(support_count_lower_bound=2, supporting_read_ids=["rA", "rB"])
        html = render_target_search_section(
            {"direct_check": {"target_search": diag}, "rounds": []},
            report_dir,
        )

        assert "&ge; 2" in html
        assert "distinct read IDs (lower bound)" in html
        assert "独立 read ID（下界）" in html
        assert "not exact totals" in html
        assert "不是精确总数" in html
        assert "exact total" not in html.lower().replace("not exact total", "")

    def test_absent_data_shows_unavailable_not_zero_or_imputed(self):
        tmp_path = self.tmp_path
        report_dir = tmp_path / "visual.report"
        report_dir.mkdir()

        # Case 1: Legacy run with no target_search diagnostics at all
        legacy_data = {
            "direct_check": {"checked": True, "success": False},
            "rounds": [
                {"round": 1, "extension_status": "extended", "connection_found": False},
            ],
        }
        legacy_html = render_target_search_section(legacy_data, report_dir)
        assert "Unavailable" in legacy_html
        assert "不可用" in legacy_html
        assert "legacy result or incomplete output" in legacy_html
        assert "None (no cap)" not in legacy_html
        assert "Connected" not in legacy_html
        assert "&ge; " not in legacy_html

        # Case 2: Partial diagnostics dict missing keys (must show Unavailable, NOT 0 or None (no cap))
        partial_diag = {
            "schema": "telomere-target-search-v1",
            "reason": "no_eligible_alignment",
            "retrieved_records": 5,
        }
        partial_html = render_target_search_section(
            {"direct_check": {"target_search": partial_diag}, "rounds": []},
            report_dir,
        )
        assert "Unavailable" in partial_html
        assert "None (no cap)" not in partial_html
        # Absent support_count_lower_bound / tail / junction / batches must not become 0
        assert "(lower bound; below min support)" not in partial_html
        assert "0 bp" not in partial_html

    def test_escaping_and_safe_relative_links(self):
        tmp_path = self.tmp_path
        report_dir = tmp_path / "visual.report"
        report_dir.mkdir()

        # Path with spaces, '#', '?', and colon-prefixed directory name
        weird_dir = tmp_path / "part2.chr.end.job" / "javascript:Chr 1#L?x" / "direct.check"
        xss = '<script>alert("xss")</script>'
        img_xss = '"><img src=x onerror=alert(1)>'
        diag = _sample_diagnostics(
            strategy=f"strat-{xss}",
            reason=f"custom_{img_xss}",
            winner_read_id=xss,
            supporting_read_ids=[xss, img_xss],
            rejection_counts={f"rej_{xss}": 3},
            support_covers=f"covers_{img_xss}",
        )
        diag_path = _write_check_dir(
            weird_dir,
            diag,
            attempted=[
                {
                    "read_id": img_xss,
                    "number": 1,
                    "status": f"st_{xss}",
                    "supporting_read_ids": [xss],
                    "support_is_lower_bound": True,
                }
            ],
        )
        ev = load_check_evidence(weird_dir)
        html = render_target_search_section(
            {
                "direct_check": {
                    "target_search": ev["diagnostics"],
                    "support_evidence": ev["support"],
                    "candidates_summary": ev["candidates_summary"],
                    "evidence_files": ev["evidence_files"],
                },
                "rounds": [],
            },
            report_dir,
        )

        assert "<script>" not in html
        assert "<img " not in html
        assert '"><img' not in html
        assert "&lt;script&gt;" in html
        assert "&lt;img src=x onerror=alert(1)&gt;" in html

        # Double-check data-i18n-* attributes stay inert even after HTML-attribute decoding
        for match in re.finditer(r'data-i18n-(?:en|zh)="([^"]*)"', html):
            decoded = html_mod.unescape(match.group(1))
            assert "<script>" not in decoded
            assert "<img " not in decoded

        # Relative link is properly percent-encoded and never absolute or javascript: scheme
        href = safe_relative_href(diag_path, report_dir)
        assert href is not None
        assert href.startswith("../part2.chr.end.job/")
        assert "javascript%3AChr%201%23L%3Fx" in href
        assert "file://" not in href
        assert not href.startswith("/")
        assert not href.startswith("javascript:")
        assert f'href="{href}"' in html

    def test_no_external_dependencies_offline(self):
        tmp_path = self.tmp_path
        # AST check: TelomereReport.py only imports standard library + ReportStyle
        source = (BIN_DIR / "TelomereReport.py").read_text(encoding="utf-8")
        tree = ast.parse(source)
        allowed_roots = {
            "collections",
            "csv",
            "json",
            "os",
            "pathlib",
            "typing",
            "urllib",
            "ReportStyle",
        }
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    assert alias.name.split(".")[0] in allowed_roots
            elif isinstance(node, ast.ImportFrom):
                if node.module:
                    assert node.module.split(".")[0] in allowed_roots

        # Generated HTML contains no external URLs, CDNs, or external stylesheets/scripts
        report_dir = tmp_path / "visual.report"
        report_dir.mkdir()
        data = {
            "chr_end": "Chr1.L",
            "initial_info": {"chr_name": "Chr1", "end_type": "L", "seed_length": 15000, "telo_reads_count": 20},
            "config": {"edge": 500, "kmer_num": 10, "kmer_size": 41, "max_rounds": None, "max_extension_length": None},
            "direct_check": {"checked": True, "success": True, "target_search": _sample_diagnostics()},
            "rounds": [],
            "final_info": {"success": True, "method": "direct", "total_rounds": 0, "total_extension_length": 0, "stop_reason": "direct_connection_found", "status": "success"},
            "summary_info": {},
        }
        out_file = Part2HTMLReportGenerator(data, report_dir).generate()
        html = out_file.read_text(encoding="utf-8")
        assert not re.search(r"https?://", html, re.IGNORECASE)
        assert "//cdn" not in html.lower()
        assert "<script src=" not in html.lower()
        assert '<link rel="stylesheet"' not in html.lower()

    def test_visualizer_integration_direct_and_family_and_global(self):
        tmp_path = self.tmp_path
        out_dir = tmp_path / "telseeker_out"
        (out_dir / "genome.telomere.check").mkdir(parents=True)
        (out_dir / "genome.telomere.check" / "need_extension_chr_end.txt").write_text(
            "Chr1.L\nChr2.R\nChr3.L\n", encoding="utf-8"
        )
        (out_dir / "genome.telomere.check" / "genome.telomere.check.left.2kb.fa").write_text(
            ">Chr1_L\n" + "ACGT" * 500 + "\n>Chr3_L\n" + "ACGT" * 500 + "\n", encoding="utf-8"
        )
        (out_dir / "genome.telomere.check" / "genome.telomere.check.right.2kb.fa").write_text(
            ">Chr2_R\n" + "TGCA" * 500 + "\n", encoding="utf-8"
        )
        (out_dir / "part1.telo.reads").mkdir(parents=True)
        (out_dir / "part1.telo.reads" / "left.telo.reads.fa").write_text(">r1\n" + "A" * 12000 + "\n", encoding="utf-8")
        (out_dir / "part1.telo.reads" / "right.telo.reads.fa").write_text(">r2\n" + "T" * 12000 + "\n", encoding="utf-8")
        (out_dir / "part3.integration.results" / "final.genome.telomere.check").mkdir(parents=True)

        part2 = out_dir / "part2.chr.end.job"

        # --- Task 1: Chr1.L (Direct connection success with full new evidence) ---
        chr1_dir = part2 / "Chr1.L"
        chr1_dir.mkdir(parents=True)
        (chr1_dir / "seed.fa").write_text(">Chr1.L_seed\n" + "A" * 20000 + "\n", encoding="utf-8")
        (chr1_dir / "telo.reads.fa").write_text(">direct_read_1\n" + "A" * 15000 + "\n", encoding="utf-8")
        chr1_direct = chr1_dir / "direct.check"
        chr1_direct.mkdir()
        (chr1_direct / "direct.linker.fa").write_text(">chromosome_linked\n" + "A" * 22000 + "\n", encoding="utf-8")
        (chr1_dir / "linker.fa").write_text(">chromosome_linked\n" + "A" * 22000 + "\n", encoding="utf-8")
        (chr1_dir / "linker.info").write_text(
            "Status: Direct connection found\nExtension Rounds: 0\nConnected Read: direct_read_1\n",
            encoding="utf-8",
        )
        (chr1_direct / "direct.check.log").write_text(
            "Status: SUCCESS - Direct connection found\n"
            "Connected Read: direct_read_1\n"
            "Alignment Length: 11000 bp\n"
            "Identity: 99.50%\n"
            "Reads Checked: 12\n"
            "Total Telomeric Reads: 500\n"
            "Reference Positions: 1-11000\n"
            "Query Positions: 4001-15000\n",
            encoding="utf-8",
        )
        d1 = _sample_diagnostics(
            end_type="L",
            pool_records=500,
            indexed_records=495,
            distinct_read_ids=490,
            retrieved_records=12,
            refined_records=12,
            eligible_records=3,
            winner_read_id="direct_read_1",
            supporting_read_ids=["direct_read_1", "direct_read_2"],
            support_count_lower_bound=2,
            reason="connected",
        )
        _write_check_dir(
            chr1_direct,
            d1,
            attempted=[
                {
                    "read_id": "direct_read_1",
                    "number": 1,
                    "status": "supported",
                    "supporting_read_ids": ["direct_read_1", "direct_read_2"],
                    "support_is_lower_bound": True,
                }
            ],
            candidates_rows=[("direct_read_1", 1, "99.5", 11000, 10, "supported")],
        )

        # --- Task 2: Chr2.R (Direct failed + Family extension with negative round, no_extension, and connected round) ---
        chr2_dir = part2 / "Chr2.R"
        chr2_dir.mkdir(parents=True)
        (chr2_dir / "seed.fa").write_text(">Chr2.R_seed\n" + "C" * 20000 + "\n", encoding="utf-8")
        (chr2_dir / "telo.reads.fa").write_text(">ext_read_1\n" + "C" * 14000 + "\n", encoding="utf-8")
        chr2_direct = chr2_dir / "direct.check"
        chr2_direct.mkdir()
        (chr2_direct / "direct.check.log").write_text("Status: FAILED - No direct connection found\n", encoding="utf-8")
        _write_check_dir(
            chr2_direct,
            _sample_diagnostics(
                end_type="R",
                pool_records=800,
                indexed_records=790,
                distinct_read_ids=780,
                retrieved_records=9,
                refined_records=9,
                alignment_blocks=11,
                eligible_records=0,
                rejection_counts={"nonterminal_alignment": 11},
                reason="no_eligible_alignment",
                winner_read_id=None,
                support_checks=0,
                supporting_read_ids=[],
                support_count_lower_bound=0,
            ),
        )

        ext2 = chr2_dir / "extension"
        r1_cc = ext2 / "attempts" / "0" / "1" / "telomere" / "connection_check"
        r1_diag = _sample_diagnostics(
            end_type="R",
            pool_records=800,
            indexed_records=790,
            distinct_read_ids=780,
            retrieved_records=14,
            refined_records=14,
            alignment_blocks=16,
            eligible_records=1,
            rejection_counts={"short_alignment": 15},
            reason="insufficient_read_support",
            winner_read_id=None,
            support_checks=1,
            supporting_read_ids=[],
            support_count_lower_bound=0,
        )
        r1_diag_path = _write_check_dir(
            r1_cc,
            r1_diag,
            attempted=[
                {
                    "read_id": "single_hit_r1",
                    "number": 4,
                    "status": "insufficient_support",
                    "supporting_read_ids": ["single_hit_r1"],
                    "support_is_lower_bound": True,
                }
            ],
        )

        r3_cc = ext2 / "attempts" / "0" / "3" / "telomere" / "connection_check"
        r3_seq = ext2 / "attempts" / "0" / "3" / "extended.fa"
        r3_seq.parent.mkdir(parents=True, exist_ok=True)
        r3_seq.write_text(">ext3\n" + "C" * 32000 + "\n", encoding="utf-8")
        r3_linker = r3_cc / "round3.linker.fa"
        r3_cc.mkdir(parents=True, exist_ok=True)
        r3_linker.write_text(">chromosome_linked\n" + "C" * 36000 + "\n", encoding="utf-8")
        r3_diag = _sample_diagnostics(
            end_type="R",
            pool_records=800,
            indexed_records=790,
            distinct_read_ids=780,
            retrieved_records=22,
            refined_records=22,
            alignment_blocks=28,
            eligible_records=3,
            rejection_counts={"low_identity": 4},
            reason="connected",
            winner_read_id="ext_read_1",
            support_checks=1,
            supporting_read_ids=["ext_read_1", "ext_read_2"],
            support_count_lower_bound=2,
        )
        r3_diag_path = _write_check_dir(
            r3_cc,
            r3_diag,
            attempted=[
                {
                    "read_id": "ext_read_1",
                    "number": 7,
                    "status": "supported",
                    "supporting_read_ids": ["ext_read_1", "ext_read_2"],
                    "support_is_lower_bound": True,
                }
            ],
            candidates_rows=[
                ("ext_read_1", 7, "99.2", 12000, 15, "supported"),
                ("ext_read_2", 8, "98.8", 11800, 22, "ranked_below_winner"),
            ],
        )

        # Signed telomere-result.json matching FamilyTelomere.load_telomere_result contract
        artifact_paths = [r1_diag_path, r3_diag_path, r3_seq, r3_linker]
        payload = {
            "schema": "telomere-family-result-v1",
            "result": {
                "chr_end": "Chr2.R",
                "success": True,
                "rounds": 3,
                "total_extension_length": 12000,
                "final_seq_file": str(r3_seq),
                "connection_info": {"read_id": "ext_read_1", "linker_file": str(r3_linker)},
                "stop_reason": "telomere_connected",
            },
            "rounds": [
                {
                    "round": 1,
                    "extension_length": 6000,
                    "cumulative_extension": 6000,
                    "extension_status": "extended",
                    "connection_found": False,
                    "target_search": r1_diag,
                    "target_evidence": str(r1_diag_path),
                },
                {
                    "round": 2,
                    "extension_length": 0,
                    "cumulative_extension": 6000,
                    "extension_status": "no_extension",
                    "connection_found": False,
                },
                {
                    "round": 3,
                    "extension_length": 6000,
                    "cumulative_extension": 12000,
                    "extension_status": "extended",
                    "extended_seq_file": str(r3_seq),
                    "connection_found": True,
                    "connected_read": "ext_read_1",
                    "identity": 99.2,
                    "aln_length": 12000,
                    "linker_file": str(r3_linker),
                    "checked_reads": 22,
                    "total_reads": 800,
                    "ref_start": 20001,
                    "ref_end": 32000,
                    "qry_start": 1,
                    "qry_end": 12000,
                    "target_search": r3_diag,
                    "target_evidence": str(r3_diag_path),
                    "supporting_read_ids": ["ext_read_1", "ext_read_2"],
                    "support_count_lower_bound": 2,
                },
            ],
            "config": {
                "max_rounds": None,
                "max_extension_length": None,
                "edge": 500,
                "kmer_num": 10,
                "kmer_size": 41,
            },
            "artifacts": {str(p.resolve()): content_hash(p) for p in artifact_paths},
        }
        payload["payload_sha256"] = fingerprint(payload)
        (ext2 / "telomere-result.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
        (chr2_dir / "linker.fa").write_text(">chromosome_linked\n" + "C" * 36000 + "\n", encoding="utf-8")
        (chr2_dir / "linker.info").write_text(
            "Status: Extension connection found\nExtension Rounds: 3\nConnected Read: ext_read_1\n",
            encoding="utf-8",
        )

        # --- Task 3: Chr3.L (Legacy task without any target-search diagnostics) ---
        chr3_dir = part2 / "Chr3.L"
        (chr3_dir / "direct.check").mkdir(parents=True)
        (chr3_dir / "direct.check" / "direct.check.log").write_text("Status: FAILED - No direct connection found\n", encoding="utf-8")
        (chr3_dir / "extension" / "round1" / "connection_check").mkdir(parents=True)
        (chr3_dir / "extension" / "round1" / "round.log").write_text("Status: Success\nExtension Length: 2500 bp\n", encoding="utf-8")
        (chr3_dir / "extension" / "round1" / "connection_check" / "round1.check.log").write_text(
            "Status: FAILED - No connection found\nReads Checked: 5/20\n", encoding="utf-8"
        )
        (chr3_dir / "extension" / "extension.summary").write_text(
            "Total Rounds: 1\nTotal Extension Length: 2500 bp\nStop Reason: no_extension\nConnection Found: No\n",
            encoding="utf-8",
        )

        # Shared hidden cache directory must never be treated as a chromosome task
        (part2 / ".telomere-target-cache" / "dummy_index_entry").mkdir(parents=True)

        # 1. Verify parser preserves new + legacy structures
        parsed1 = TelSeekerPart2LogParser(chr1_dir).parse()
        assert parsed1["direct_check"]["target_search"]["reason"] == "connected"
        assert parsed1["direct_check"]["support_count_lower_bound"] == 2
        assert "diagnostics" in parsed1["direct_check"]["evidence_files"]

        parsed2 = TelSeekerPart2LogParser(chr2_dir).parse()
        assert parsed2["direct_check"]["target_search"]["reason"] == "no_eligible_alignment"
        assert len(parsed2["rounds"]) == 3
        assert parsed2["rounds"][0]["target_search"]["reason"] == "insufficient_read_support"
        assert parsed2["rounds"][2]["target_search"]["reason"] == "connected"
        assert "candidates" in parsed2["rounds"][2]["evidence_files"]

        parsed3 = TelSeekerPart2LogParser(chr3_dir).parse()
        assert "target_search" not in parsed3["direct_check"]
        assert "target_search" not in parsed3["rounds"][0]

        # 2. Generate all per-task reports
        ok, total = generate_all_task_reports(out_dir)
        assert (ok, total) == (3, 3)

        visual_dir = out_dir / "visual.report"
        assert not (visual_dir / ".telomere-target-cache.report.html").exists()
        chr1_html = (visual_dir / "Chr1.L.report.html").read_text(encoding="utf-8")
        assert "Telomere Target Search &amp; Corroboration Evidence" in chr1_html
        assert "端粒目标检索与佐证证据" in chr1_html
        assert "../part2.chr.end.job/Chr1.L/direct.check/connection-diagnostics.json" in chr1_html
        assert "&ge; 2" in chr1_html
        assert "<svg" in chr1_html  # existing visualization preserved
        assert 'class="language-btn"' in chr1_html  # existing language switcher preserved

        chr2_html = (visual_dir / "Chr2.R.report.html").read_text(encoding="utf-8")
        assert "no_eligible_alignment" in chr2_html
        assert "insufficient_read_support" in chr2_html
        assert "single_hit_r1" in chr2_html
        assert "Not checked (no extension)" in chr2_html
        assert "../part2.chr.end.job/Chr2.R/extension/attempts/0/3/telomere/connection_check/candidates.tsv" in chr2_html
        assert "Round Details (1 vs N Check)" in chr2_html  # existing round table preserved
        assert "<svg" in chr2_html

        chr3_html = (visual_dir / "Chr3.L.report.html").read_text(encoding="utf-8")
        assert "legacy result or incomplete output" in chr3_html
        assert "Unavailable" in chr3_html

        # 3. Generate global report and verify per-end evidence summaries + backward compatibility
        global_path = visual_dir / "Global.report.html"
        generate_html_visualization(out_dir, global_path, motif="TTAGGG")
        global_html = global_path.read_text(encoding="utf-8")

        assert "Linked Sequence Length" in global_html and "拼接后序列长度" in global_html
        assert "<th>Extension Length</th>" not in global_html
        assert ".telomere-target-cache" not in global_html
        assert "Target Evidence" in global_html and "目标检索证据" in global_html
        assert "renderTargetEvidenceCell" in global_html
        assert "lower bound" in global_html and "下界" in global_html
        assert "../part2.chr.end.job/Chr1.L/direct.check/connection-diagnostics.json" in global_html
        assert "../part2.chr.end.job/Chr2.R/extension/attempts/0/3/telomere/connection_check/connection-diagnostics.json" in global_html
        assert 'id="extension-table"' in global_html
        assert 'id="initial-review-table"' in global_html
        assert 'id="final-review-table"' in global_html
        assert 'id="chromosome-task-tracks"' in global_html
        assert not re.search(r"https?://", global_html, re.IGNORECASE)

        # Verify parse_part2_results backward compatibility (1-arg call) and summary semantics
        rows_compat = parse_part2_results(part2)
        assert len(rows_compat) == 3
        by_end = {r["chr_end"]: r for r in rows_compat}
        assert ".telomere-target-cache" not in by_end
        if HAS_BIOPYTHON:
            assert by_end["Chr1.L"]["extension_length"] == 22000  # numeric compatibility kept
            assert by_end["Chr2.R"]["extension_length"] == 36000
        assert by_end["Chr1.L"]["target_evidence_summary"]["available"] is True
        assert by_end["Chr1.L"]["target_evidence_summary"]["support_count_lower_bound"] == 2
        assert by_end["Chr1.L"]["target_evidence_summary"]["alignment_batches"] == 2
        assert by_end["Chr1.L"]["target_evidence_summary"]["selected_tail_bases"] == 3420
        assert by_end["Chr1.L"]["target_evidence_summary"]["junction_bases_tested"] == 200
        assert by_end["Chr2.R"]["target_evidence_summary"]["available"] is True
        assert by_end["Chr2.R"]["target_evidence_summary"]["checks_with_diagnostics"] == 3
        assert by_end["Chr3.L"]["target_evidence_summary"]["available"] is False
        assert by_end["Chr3.L"]["target_evidence_summary"]["pool_records"] is None
        assert by_end["Chr3.L"]["target_evidence_summary"]["alignment_batches"] is None
        assert by_end["Chr3.L"]["target_evidence_summary"]["selected_tail_bases"] is None
        assert by_end["Chr3.L"]["target_evidence_summary"]["junction_bases_tested"] is None
        assert by_end["Chr3.L"]["target_evidence_summary"]["support_count_lower_bound"] is None

        # Direct helper check on summarize_chr_end_evidence with report_dir
        summary2 = summarize_chr_end_evidence(chr2_dir, report_dir=visual_dir)
        assert summary2["diagnostics_href"] == "../part2.chr.end.job/Chr2.R/extension/attempts/0/3/telomere/connection_check/connection-diagnostics.json"


if __name__ == "__main__":
    unittest.main()
