"""Lightweight presentation contracts: offline HTML, evidence and long labels."""
from html.parser import HTMLParser
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from FamilyResult import render_report
from FamilyCtgLinker import render_ctg_report
from TelSeekerVisualizer import Part2HTMLReportGenerator
from AutoGapfillerVisualizer import GapFillingVisualizer


class ReportHTML(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.tags = []
        self.text = []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))

    def handle_data(self, value):
        self.text.append(value)


class ReportPresentationTests(unittest.TestCase):
    def assert_portable(self, html):
        parsed = ReportHTML(html)
        self.assertTrue(any(tag == "meta" and attrs.get("name") == "viewport" for tag, attrs in parsed.tags))
        self.assertEqual(sum(tag == "h1" for tag, attrs in parsed.tags), 1)
        for tag, attrs in parsed.tags:
            if tag in ("script", "link", "img"):
                url = attrs.get("src", attrs.get("href", ""))
                self.assertFalse(url.startswith(("https://", "http://", "//")), url)
        self.assertIn("@media", html)
        return parsed

    def family(self, status="closed"):
        families = {}
        for key, direction in (("L", "rightward"), ("R", "leftward")):
            versions = [SimpleNamespace(number=n, length=20000 + n * 5000, sha256=str(n) * 64)
                        for n in range(3)]
            families[key] = SimpleNamespace(key=key, direction=direction, attempts=2, added_bases=10000,
                                            stop_reason=None, versions=versions)
        winner = {"left_family": "L", "left_version": 2, "right_family": "R", "right_version": 1,
                  "identity": 99.75, "overlap": 12345, "overhang": 27}
        alternate = dict(winner, right_version=2, identity=98.5)
        search = SimpleNamespace(status=status, total_rounds=4, total_length=20000, families=families,
                                 checked={"one": [winner, alternate]})
        result = {"winner": winner if status == "closed" else None, "final_sequence": "sample.final.fa"}
        return search, result

    def test_family_report_preserves_all_versions_candidates_and_downloads(self):
        search, result = self.family()
        html = render_report(search, result)
        parsed = self.assert_portable(html)
        text = " ".join(parsed.text)
        for value in ("99.75", "98.5", "12345", "27", "0" * 64, "1" * 64, "2" * 64):
            self.assertIn(value, text)
        self.assertIn("闭合", text)
        # Three table headers + two families + two candidates + six versions.
        self.assertEqual(sum(tag == "tr" for tag, attrs in parsed.tags), 13)
        links = {attrs["href"] for tag, attrs in parsed.tags if tag == "a" and "href" in attrs}
        self.assertTrue({"sample.final.fa", "process.agp", "family-components.fa", "family-result.json",
                         "families.json"}.issubset(links))

    def test_unclosed_report_keeps_stop_reason_and_does_not_invent_success(self):
        search, result = self.family("round_limit")
        search.checked = {}
        search.families["L"].stop_reason = '<not enough reads & evidence>'
        html = render_report(search, result)
        parsed = self.assert_portable(html)
        self.assertIn("round_limit", " ".join(parsed.text))
        self.assertIn("<not enough reads & evidence>", " ".join(parsed.text))
        self.assertNotIn("<not enough reads & evidence>", html)
        self.assertIn("未闭合", html)

    def test_ctg_report_keeps_members_strands_and_endpoint_evidence(self):
        name = 'contig<evidence>&"' + "long_" * 30
        result = {"input_contigs": [name], "scaffolds": [{"scaffold_id": "Scaffold_1", "length": 50000,
                  "members": [{"id": name, "strand": "-"}]}],
                  "events": [{"number": 1, "primary": name, "direction": "right", "status": "closed",
                              "selected": {"assembly": "candidate_2", "strand": "+"},
                              "total_rounds": 4, "result_dir": "project/endpoint-000001"}],
                  "total_rounds": 4, "total_extension_length": 18000}
        html = render_ctg_report(result)
        parsed = self.assert_portable(html)
        text = " ".join(parsed.text)
        self.assertIn(name, text)
        self.assertNotIn(name, html)
        self.assertIn("candidate_2", text)
        self.assertIn("-", text)
        links = {attrs["href"] for tag, attrs in parsed.tags if tag == "a" and "href" in attrs}
        self.assertIn("project/endpoint-000001/family-report.html", links)
        self.assertIn("DG.scaffold.components.fa", links)

    def test_auto_global_keeps_language_and_chart_hooks_without_network(self):
        visualizer = GapFillingVisualizer.__new__(GapFillingVisualizer)
        visualizer.gap_filling_results = {"summary": {"total_gaps": 0, "filled_gaps": 0, "fill_rate": "0%"}, "gaps": []}
        visualizer.gaps_by_chrom = {}
        visualizer.preprocessing_config = None
        visualizer.integration_config = None
        html = visualizer._create_full_html()
        self.assert_portable(html)
        self.assertIn("switchLanguage", html)
        self.assertIn("lang-en", html)
        self.assertIn("lang-zh", html)
        self.assertIn('class="lang-en">Filled</th>', html)
        self.assertIn('class="lang-en">Direction</th>', html)
        self.assertNotIn('class="lang-en">Extended</th>', html)
        self.assertNotIn('class="lang-en">Flag</th>', html)

    def test_telomere_detail_keeps_round_values_and_language_hooks(self):
        data = {"chr_end": "chr01.R", "initial_info": {}, "config": {},
                "final_info": {"success": False, "method": "none", "status": "max_limit_reached",
                               "total_rounds": 1, "total_extension_length": 1234, "stop_reason": "max_rounds_reached"},
                "rounds": [{"round": 1, "extension_length": 1234, "cumulative_extension": 1234,
                            "extension_status": "extended", "connection_found": False}], "direct_check": {}}
        with tempfile.TemporaryDirectory() as tmp:
            html = Part2HTMLReportGenerator(data, Path(tmp))._build_html()
        self.assert_portable(html)
        self.assertIn("data-i18n-en", html)
        self.assertIn("data-i18n-zh", html)
        self.assertIn("chr01.R", html)
        self.assertIn("max_rounds_reached", html)
        self.assertTrue("1234" in html or "1,234" in html)


    def test_telomere_limits_and_identifiers_survive_language_switch(self):
        name = 'chr<special>&".R'
        read_id = 'read<special>&"'
        data = {"chr_end": name, "initial_info": {"chr_name": name, "end_type": "R"},
                "config": {"max_rounds": None, "max_extension_length": None},
                "final_info": {"success": False, "method": "none", "total_rounds": 1},
                "rounds": [{"round": 1, "connected_read": read_id,
                            "extension_status": "extended", "connection_found": True}],
                "direct_check": {"success": True, "read_id": read_id}}
        with tempfile.TemporaryDirectory() as tmp:
            html = Part2HTMLReportGenerator(data, Path(tmp))._build_html()
        parsed = self.assert_portable(html)
        self.assertIn(name, " ".join(parsed.text))
        self.assertIn(read_id, " ".join(parsed.text))
        self.assertFalse(any(tag == "special" for tag, _ in parsed.tags))
        unlimited = [attrs for _, attrs in parsed.tags if attrs.get("data-i18n-en") == "No limit"]
        self.assertEqual(len(unlimited), 2)
        self.assertTrue(all(attrs["data-i18n-zh"] == "无限制" for attrs in unlimited))
        # Translation strings are assigned to innerHTML: parse that second layer,
        # not merely the initial HTML, to catch single-escaped attribute values.
        task = next(attrs for _, attrs in parsed.tags if attrs.get("class") == "sub")
        for lang in ("en", "zh"):
            translated = ReportHTML(task["data-i18n-" + lang])
            self.assertIn(name, " ".join(translated.text))
            self.assertEqual([tag for tag, _ in translated.tags], ["code"])
        self.assertIn("Report generated by TelSeekerVisualizer.py.", html)
        self.assertNotIn("telo_part2_visual.py", html)

    def test_telomere_zero_limits_are_not_rendered_as_unlimited(self):
        data = {"config": {"max_rounds": 0, "max_extension_length": 0},
                "final_info": {}, "rounds": [], "direct_check": {}}
        with tempfile.TemporaryDirectory() as tmp:
            html = Part2HTMLReportGenerator(data, Path(tmp))._build_html()
        self.assertNotIn('data-i18n-en="No limit"', html)
        self.assertIn('Max Rounds</td><td>0</td>', html)
        self.assertIn('Max Extension Length</td><td>0</td>', html)

    def test_minimal_layout_folds_secondary_details_and_removes_boilerplate(self):
        search, result = self.family()
        fam_html = render_report(search, result)
        self.assertIn('data-lang="en"', fam_html)
        self.assertIn('class="summary-strip"', fam_html)
        self.assertIn('<details class="fold">', fam_html)
        self.assertNotIn("Only a closed result joins both families", fam_html)
        self.assertNotIn("Read mappings, local assemblies and their receipts remain under", fam_html)

        ctg_html = render_ctg_report({
            "input_contigs": ["c1"],
            "scaffolds": [{"scaffold_id": "Scaffold_1", "length": 1000, "members": [{"id": "c1", "strand": "+"}]}],
            "events": [],
            "total_rounds": 0,
        })
        self.assertIn('data-lang="en"', ctg_html)
        self.assertIn('class="summary-strip"', ctg_html)
        self.assertNotIn("Each endpoint search shares one budget", ctg_html)

        visualizer = GapFillingVisualizer.__new__(GapFillingVisualizer)
        visualizer.gap_filling_results = {"summary": {"total_gaps": 1, "filled_gaps": 1, "fill_rate": "100%"}, "gaps": []}
        visualizer.gaps_by_chrom = {}
        visualizer.preprocessing_config = None
        visualizer.integration_config = None
        auto_html = visualizer._create_full_html()
        self.assertLess(auto_html.index('id="chromVisualization"'), auto_html.index('id="gapFillingTasksTable"'))
        self.assertNotIn("No reads information available", auto_html)
        self.assertIn('<details class="fold">', auto_html)

        # TelSeeker detail omits empty Connection Visualization placeholder when no SVG exists
        with tempfile.TemporaryDirectory() as tmp:
            telo_html = Part2HTMLReportGenerator(
                {"chr_end": "Chr1.L", "initial_info": {}, "config": {}, "final_info": {}, "rounds": [], "direct_check": {}},
                Path(tmp),
            )._build_html()
        self.assertNotIn("No successful connection available for visualization", telo_html)
        self.assertIn('<details class="fold">', telo_html)

    def test_family_report_renders_interactive_extension_process_visual(self):
        search, result = self.family()
        result["winner"].update({
            "left_start": 17655,
            "left_end": 30000,
            "right_start": 0,
            "right_end": 12345,
            "evidence": {"stage": "full"},
        })
        result["junction"] = {
            "gap_length": 450,
            "overlap_trim": 0,
            "flank_length_delta": 0,
            "repairs": [{"side": "L", "start0": 19980, "status": "accepted"}],
        }
        result["connection_diagnostics"] = [
            {"left_family": "L", "left_version": 0, "right_family": "R", "right_version": 0, "status": "rejected", "reason": "no_alignment"},
            {"left_family": "L", "left_version": 2, "right_family": "R", "right_version": 1, "status": "eligible", "reason": "winner"},
        ]
        html = render_report(search, result)
        parsed = self.assert_portable(html)
        self.assertIn("data-viz-player", html)
        self.assertIn('class="viz-svg"', html)
        self.assertIn('data-step="0"', html)
        self.assertIn('data-viz-view="full"', html)
        self.assertIn("is-hidden", html)
        self.assertIn("lbl-pause", html)
        self.assertIn("Winning pair synteny", html)
        self.assertIn("Final junction composition", html)
        self.assertIn("detail-grid", html)
        self.assertIn("viz-matrix", html)
        self.assertIn("hash-short", html)
        # Visual section itself uses SVG + CSS Grid, never adding extra <table> rows.
        viz_start = html.index("data-viz-player")
        viz_end = html.index("</section>", viz_start)
        self.assertNotIn("<tr", html[viz_start:viz_end])

        # Closed job with preserved flanks (0 repairs, 0 overlap_trim) skips redundant junction composition bar.
        plain_search, plain_result = self.family()
        plain_result["junction"] = {"mode": "gap", "gap_length": 300, "flank_length_delta": 0, "overlap_trim": 0, "repairs": []}
        plain_html = render_report(plain_search, plain_result)
        self.assertNotIn("Final junction composition", plain_html)

        u_search, u_result = self.family("exhausted")
        u_search.checked = {}
        u_search.families["L"].versions = [SimpleNamespace(number=0, length=50000, sha256="0" * 64)]
        u_search.families["L"].attempts = 2
        u_search.families["L"].added_bases = 0
        u_search.families["R"].versions = [
            SimpleNamespace(number=0, length=50000, sha256="0" * 64),
            SimpleNamespace(number=1, length=50001, sha256="1" * 64),
        ]
        u_search.families["R"].attempts = 2
        u_search.families["R"].added_bases = 1
        u_html = render_report(u_search, u_result)
        self.assert_portable(u_html)
        self.assertIn("data-viz-player", u_html)
        self.assertIn("Unclosed · span unknown", u_html)
        self.assertIn("未闭合 · 间距未知", u_html)
        self.assertNotIn("Winning pair synteny", u_html)
        self.assertNotIn("Eligible Connections", u_html)

if __name__ == "__main__":
    unittest.main()
