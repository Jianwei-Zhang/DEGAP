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

if __name__ == "__main__":
    unittest.main()
