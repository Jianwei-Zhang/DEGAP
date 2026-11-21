"""Coordinate projection and track alignment tests for AutoGapfiller HTML report."""
from html.parser import HTMLParser
import json
from pathlib import Path
import shutil
import subprocess
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from AutoGapfillerVisualizer import GapFillingVisualizer


class ScriptExtractor(HTMLParser):
    """Extract inline <script> blocks and structural tags from rendered HTML."""

    def __init__(self, text):
        super().__init__()
        self.tags = []
        self.scripts = []
        self._in_script = False
        self._chunks = []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))
        if tag == "script":
            self._in_script = True
            self._chunks = []

    def handle_endtag(self, tag):
        if tag == "script" and self._in_script:
            self.scripts.append("".join(self._chunks))
            self._in_script = False
            self._chunks = []

    def handle_data(self, data):
        if self._in_script:
            self._chunks.append(data)


def render_fixture_html(gaps, chrom_lengths=None):
    chrom_lengths = chrom_lengths or {}
    gaps_by_chrom = {}
    filled = 0
    for gap in gaps:
        chrom = gap["chromosome"]
        bucket = gaps_by_chrom.setdefault(
            chrom,
            {"gaps": [], "total": 0, "filled": 0, "unfilled": 0, "length": chrom_lengths.get(chrom, 0)},
        )
        bucket["gaps"].append(gap)
        bucket["total"] += 1
        if gap.get("filled"):
            bucket["filled"] += 1
            filled += 1
        else:
            bucket["unfilled"] += 1
    for chrom, data in gaps_by_chrom.items():
        if data["length"] == 0:
            data["length"] = max((g.get("gap_end", 0) for g in data["gaps"]), default=0)

    visualizer = GapFillingVisualizer.__new__(GapFillingVisualizer)
    total = len(gaps)
    rate = f"{(filled / total * 100):.1f}%" if total else "0%"
    visualizer.gap_filling_results = {
        "summary": {"total_gaps": total, "filled_gaps": filled, "fill_rate": rate},
        "gaps": gaps,
    }
    visualizer.gaps_by_chrom = gaps_by_chrom
    visualizer.preprocessing_config = None
    visualizer.integration_config = None
    return visualizer._create_full_html()


def run_report_js(html, expression):
    """Evaluate an expression against the rendered report's inline JS using Node vm."""
    extractor = ScriptExtractor(html)
    assert extractor.scripts, "No inline <script> found in rendered AutoGapfiller HTML"
    report_script = extractor.scripts[-1]
    harness = (
        "const vm = require('vm');\n"
        "const payload = JSON.parse(require('fs').readFileSync(0, 'utf8'));\n"
        "const ctx = vm.createContext({\n"
        "  document: { addEventListener: () => {} },\n"
        "  console,\n"
        "});\n"
        "vm.runInContext(payload.script, ctx);\n"
        "const result = vm.runInContext('(' + payload.expr + ')', ctx);\n"
        "process.stdout.write(JSON.stringify(result));\n"
    )
    proc = subprocess.run(
        ["node", "-e", harness],
        input=json.dumps({"script": report_script, "expr": expression}),
        text=True,
        capture_output=True,
        check=True,
    )
    return json.loads(proc.stdout)


class AutoGapfillerReportCoordinateTests(unittest.TestCase):
    def test_rendered_html_preserves_hooks_and_bilingual_original_coordinate_frame(self):
        gaps = [
            {
                "chromosome": "chr1",
                "gap_id": 1,
                "gap_start": 101,
                "gap_end": 200,
                "gap_len": 100,
                "filled": True,
                "fill_length": 100,
                "fill_type": "extension",
                "fill_direction": "bidirectional",
                "rounds": 2,
                "junction": {"mode": "gap", "gap_length": 100, "flank_length_delta": 5, "overlap_trim": 0},
            },
            {
                "chromosome": "chr1",
                "gap_id": 2,
                "gap_start": 501,
                "gap_end": 600,
                "gap_len": 100,
                "filled": False,
                "fill_direction": "none",
            },
        ]
        html = render_fixture_html(gaps, {"chr1": 10000})
        extracted = ScriptExtractor(html)
        ids = {attrs.get("id") for _, attrs in extracted.tags if attrs.get("id")}
        self.assertTrue(
            {
                "initialGenomeTable",
                "gapFillingTasksTable",
                "finalGenomeTable",
                "chromVisualization",
                "chromSelect",
                "trackLengthControl",
                "tooltip",
                "chromCoordinateNote",
            }.issubset(ids)
        )
        self.assertIn("original assembly coordinate frame", html)
        self.assertIn("原始组装坐标系", html)
        self.assertIn("switchLanguage", html)
        self.assertIn("lang-en", html)
        self.assertIn("lang-zh", html)
        for helper in (
            "function getOriginalGapLength",
            "function getFlankLengthDelta",
            "function getFilledSequenceLength",
            "function getFilledGapOffsetDelta",
            "function compareGapsGenomic",
            "function projectFinalUnfilledGaps",
            "function getOriginalTrackPlacement",
            "function formatFillLengthCell",
            "function buildTooltipContent",
        ):
            self.assertIn(helper, html)
        # Track rendering must not substitute final fill_length for original axis width.
        self.assertNotIn("gapDisplayLength = gap.fill_length", html)
        self.assertNotIn("const originalGapLen = gap.gap_end - gap.gap_start;", html)

    @unittest.skipUnless(shutil.which("node"), "node runtime not available")
    def test_zero_based_track_placement_and_inclusive_length_fallback(self):
        html = render_fixture_html([], {})
        result = run_report_js(
            html,
            """(() => {
                const noLen = { gap_start: 101, gap_end: 200, filled: true, fill_length: 9999 };
                const singleBase = { gap_start: 101, gap_end: 101, filled: false };
                const firstBase = { gap_start: 1, gap_end: 50, gap_len: 50, filled: true, fill_length: 500 };
                const rowEdge = { gap_start: 1000, gap_end: 1000, filled: false };
                return {
                    inclusiveFallback: getOriginalGapLength(noLen),
                    singleBaseLen: getOriginalGapLength(singleBase),
                    firstPlacement: getOriginalTrackPlacement(firstBase),
                    filledUsesOriginalSpan: getOriginalTrackPlacement(noLen, 0, 1000),
                    rowEdgeInRow0: getOriginalTrackPlacement(rowEdge, 0, 1000),
                    rowEdgeInRow1: getOriginalTrackPlacement(rowEdge, 1000, 2000),
                };
            })()""",
        )
        self.assertEqual(result["inclusiveFallback"], 100)
        self.assertEqual(result["singleBaseLen"], 1)
        self.assertEqual(result["firstPlacement"]["start0"], 0)
        self.assertEqual(result["firstPlacement"]["span"], 50)
        self.assertEqual(result["filledUsesOriginalSpan"]["start0"], 100)
        self.assertEqual(result["filledUsesOriginalSpan"]["span"], 100)
        self.assertEqual(result["filledUsesOriginalSpan"]["clippedSpan"], 100)
        self.assertEqual(result["rowEdgeInRow0"]["relativeStart"], 999)
        self.assertIsNone(result["rowEdgeInRow1"])

    @unittest.skipUnless(shutil.which("node"), "node runtime not available")
    def test_no_plus_one_drift_after_equal_length_fill_and_legacy_data(self):
        html = render_fixture_html([], {})
        result = run_report_js(
            html,
            """(() => {
                const legacyGaps = [
                    { gap_id: 1, gap_start: 101, gap_end: 200, filled: true, fill_length: 100 },
                    { gap_id: 2, gap_start: 501, gap_end: 600, filled: false },
                ];
                return {
                    projected: projectFinalUnfilledGaps(legacyGaps),
                    legacyDelta: getFlankLengthDelta(legacyGaps[0]),
                    legacyTooltip: buildTooltipContent(legacyGaps[0], 'chr1'),
                };
            })()""",
        )
        self.assertEqual(result["legacyDelta"], 0)
        self.assertEqual(len(result["projected"]), 1)
        self.assertEqual(result["projected"][0]["offset"], 0)
        self.assertEqual(result["projected"][0]["finalStart"], 501)
        self.assertEqual(result["projected"][0]["finalEnd"], 600)
        self.assertIn("Original Gap Span: 100 bp", result["legacyTooltip"])
        self.assertIn("Filled Sequence Length: 100 bp", result["legacyTooltip"])
        self.assertNotIn("undefined", result["legacyTooltip"])
        self.assertNotIn("NaN", result["legacyTooltip"])

    @unittest.skipUnless(shutil.which("node"), "node runtime not available")
    def test_positive_and_negative_flank_delta_shift_later_unfilled_gap(self):
        html = render_fixture_html([], {})
        result = run_report_js(
            html,
            """(() => {
                const posCase = projectFinalUnfilledGaps([
                    { gap_id: 1, gap_start: 101, gap_end: 200, gap_len: 100, filled: true, fill_length: 100, junction: { flank_length_delta: 12 } },
                    { gap_id: 2, gap_start: 501, gap_end: 600, gap_len: 100, filled: false },
                ]);
                const negCase = projectFinalUnfilledGaps([
                    { gap_id: 1, gap_start: 101, gap_end: 200, gap_len: 100, filled: true, fill_length: 100, junction: { flank_length_delta: -7 } },
                    { gap_id: 2, gap_start: 501, gap_end: 600, gap_len: 100, filled: false },
                ]);
                const cellPos = formatFillLengthCell({ filled: true, fill_length: 100, junction: { flank_length_delta: 12 } });
                const tipNeg = buildTooltipContent({ gap_id: 1, gap_start: 101, gap_end: 200, gap_len: 100, filled: true, fill_length: 100, junction: { flank_length_delta: -7 } }, 'chr1');
                return { pos: posCase[0], neg: negCase[0], cellPos, tipNeg };
            })()""",
        )
        self.assertEqual(result["pos"]["offset"], 12)
        self.assertEqual(result["pos"]["finalStart"], 513)
        self.assertEqual(result["pos"]["finalEnd"], 612)
        self.assertEqual(result["neg"]["offset"], -7)
        self.assertEqual(result["neg"]["finalStart"], 494)
        self.assertEqual(result["neg"]["finalEnd"], 593)
        self.assertIn("Flank delta: +12 bp", result["cellPos"]["html"])
        self.assertIn("Flank Length Delta: -7 bp", result["tipNeg"])

    @unittest.skipUnless(shutil.which("node"), "node runtime not available")
    def test_zero_and_negative_fill_length_arithmetic_and_wording(self):
        html = render_fixture_html([], {})
        result = run_report_js(
            html,
            """(() => {
                const zeroGap = { gap_id: 1, gap_start: 101, gap_end: 200, gap_len: 100, filled: true, fill_length: 0 };
                const overlapGap = {
                    gap_id: 2, gap_start: 301, gap_end: 400, gap_len: 100, filled: true, fill_length: -150,
                    junction: { mode: 'overlap', gap_length: -150, flank_length_delta: 0, overlap_trim: 150 },
                };
                const laterUnfilled = { gap_id: 3, gap_start: 1001, gap_end: 1100, gap_len: 100, filled: false };
                return {
                    afterZero: projectFinalUnfilledGaps([zeroGap, laterUnfilled])[0],
                    afterOverlap: projectFinalUnfilledGaps([overlapGap, laterUnfilled])[0],
                    zeroCell: formatFillLengthCell(zeroGap),
                    overlapCell: formatFillLengthCell(overlapGap),
                    zeroTip: buildTooltipContent(zeroGap, 'chr<1>'),
                    overlapTip: buildTooltipContent(overlapGap, 'chr1'),
                };
            })()""",
        )
        # Zero fill removes 100 bp original gap: 1001 - 100 = 901.
        self.assertEqual(result["afterZero"]["offset"], -100)
        self.assertEqual(result["afterZero"]["finalStart"], 901)
        self.assertEqual(result["afterZero"]["finalEnd"], 1000)
        # Negative overlap (-150) removes 100 bp original gap plus 150 bp overlap trim: offset = -250.
        self.assertEqual(result["afterOverlap"]["offset"], -250)
        self.assertEqual(result["afterOverlap"]["finalStart"], 751)
        self.assertEqual(result["afterOverlap"]["finalEnd"], 850)

        self.assertIn("0 bp", result["zeroCell"]["html"])
        self.assertNotEqual(result["zeroCell"]["html"], "-")
        self.assertIn("Filled Sequence Length: 0 bp", result["zeroTip"])
        self.assertNotIn("<1>", result["zeroTip"])
        self.assertIn("&lt;1&gt;", result["zeroTip"])

        self.assertIn("-150 bp", result["overlapCell"]["html"])
        self.assertIn("Original-overlap trim (150 bp)", result["overlapCell"]["html"])
        self.assertIn("Original-Overlap Trimming: -150 bp (150 bp trimmed)", result["overlapTip"])

    @unittest.skipUnless(shutil.which("node"), "node runtime not available")
    def test_genomic_sorting_uses_numeric_gap_start_then_stable_gap_id(self):
        html = render_fixture_html([], {})
        result = run_report_js(
            html,
            """(() => {
                // gap_id=10 is upstream (start 101), gap_id=2 is downstream (start 501).
                const unsorted = [
                    { gap_id: 2, gap_start: 501, gap_end: 600, gap_len: 100, filled: false },
                    { gap_id: 10, gap_start: 101, gap_end: 200, gap_len: 100, filled: true, fill_length: 150 },
                ];
                const ties = [
                    { gap_id: 10, gap_start: 300, gap_end: 350, filled: false },
                    { gap_id: 2, gap_start: 300, gap_end: 350, filled: false },
                ];
                return {
                    projected: projectFinalUnfilledGaps(unsorted),
                    tieOrder: [...ties].sort(compareGapsGenomic).map(g => g.gap_id),
                };
            })()""",
        )
        # Upstream gap_id=10 (+50 bp) must be applied BEFORE downstream gap_id=2!
        self.assertEqual(result["projected"][0]["gap"]["gap_id"], 2)
        self.assertEqual(result["projected"][0]["offset"], 50)
        self.assertEqual(result["projected"][0]["finalStart"], 551)
        self.assertEqual(result["projected"][0]["finalEnd"], 650)
        self.assertEqual(result["tieOrder"], [2, 10])


if __name__ == "__main__":
    unittest.main()
