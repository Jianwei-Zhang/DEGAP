"""Unit tests for portable offline junction and connection diagnostics HTML helpers."""
from html.parser import HTMLParser
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from JunctionReport import render_junction_report, render_project_junctions


class FragmentHTML(HTMLParser):
    """Parse HTML fragments while separating overview text from <details> text."""

    def __init__(self, text):
        super().__init__()
        self.tags = []
        self.overview_text = []
        self.details_text = []
        self._details_depth = 0
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))
        if tag == "details":
            self._details_depth += 1

    def handle_endtag(self, tag):
        if tag == "details" and self._details_depth > 0:
            self._details_depth -= 1

    def handle_data(self, value):
        if self._details_depth > 0:
            self.details_text.append(value)
        else:
            self.overview_text.append(value)

    @property
    def all_text(self):
        return " ".join(self.overview_text + self.details_text)


class JunctionReportTests(unittest.TestCase):
    def rendered_classes(self, parsed):
        return {cls for _, attrs in parsed.tags for cls in attrs.get("class", "").split()}

    def assert_offline_accessible(self, html):
        parsed = FragmentHTML(html)
        self.assertEqual(sum(tag == "h1" for tag, _ in parsed.tags), 0)
        self.assertIn("@media", html)
        self.assertIn(":focus-visible", html)
        for tag, attrs in parsed.tags:
            self.assertNotEqual(tag, "script")
            if tag in ("link", "img", "a"):
                url = attrs.get("src", attrs.get("href", ""))
                self.assertFalse(url.startswith(("https://", "http://", "//", "javascript:", "file:")), url)
        th_tags = [attrs for tag, attrs in parsed.tags if tag == "th"]
        self.assertTrue(th_tags)
        self.assertTrue(all(attrs.get("scope") == "col" for attrs in th_tags))
        return parsed

    def test_legacy_and_empty_inputs_return_empty_string(self):
        self.assertEqual(render_junction_report(None), "")
        self.assertEqual(render_junction_report({}), "")
        self.assertEqual(
            render_junction_report({
                "status": "closed",
                "gap_length": 500,
                "winner": {"left_family": "L", "right_family": "R"},
            }),
            "",
        )
        self.assertEqual(render_junction_report({"junction": {}, "connection_diagnostics": []}), "")
        self.assertEqual(
            render_junction_report({"junction": {}, "connection_diagnostics": [], "anchor_retries": []}),
            "",
        )
        self.assertEqual(render_project_junctions(None), "")
        self.assertEqual(render_project_junctions([]), "")
        self.assertEqual(render_project_junctions({}), "")
        self.assertEqual(render_project_junctions({"events": [{"number": 1, "status": "closed"}]}), "")
        self.assertEqual(render_project_junctions({"history": [{"id": "DG1"}]}), "")

    def test_escaping_truncation_and_safe_relative_artifacts_only(self):
        long_orig = '<script>alert("o")</script>&' + ("ACGT" * 20)
        long_prop = '<img src=x onerror=1>&' + ("TGCA" * 20)
        result = {
            "junction": {
                "mode": "gap",
                "gap_length": 310,
                "flank_length_delta": 2,
                "overlap_trim": 0,
                "plan_path": "/abs/secret/composition.json",
                "policy": {"repair_window": 2000, "min_support": 5, "min_fraction": 0.8, "min_mapq": 20},
                "flanks": {
                    "L": {"status": "repaired", "proposed_variants": 1, "accepted_variants": 1, "protected_variants": 0},
                    "R": {"status": "preserved", "proposed_variants": 0, "accepted_variants": 0, "protected_variants": 0},
                },
                "repairs": [
                    {
                        "side": "L",
                        "start0": 19500,
                        "end0": 19608,
                        "original": long_orig,
                        "proposed": long_prop,
                        "status": "accepted",
                        "reason": 'independent_spanning_reads<"&>',
                        "alt_reads": 9,
                        "ref_reads": 0,
                        "other_reads": 0,
                        "spanning_reads": 9,
                        "alt_fraction": 1.0,
                        "strands": ["+", "-"],
                    }
                ],
                "warnings": ['Warn <b>check</b> & "review"'],
                "evidence_artifacts": {
                    "Absolute": "/tmp/composition.json",
                    "Traversal": "../secret.json",
                    "Scheme": "javascript:alert(1)",
                    "WinPath": "C:\\temp\\bad.json",
                    "Relative receipt": "evidence/junction-summary.json",
                },
            },
            "connection_diagnostics": [
                {
                    "left_family": 'L<bad>&"',
                    "left_version": 0,
                    "right_family": "R",
                    "right_version": 1,
                    "stages": [{"stage": "fast", "blocks": 2, "candidates": 1, "reason": "eligible<tag>"}],
                    "status": "eligible",
                    "reason": "eligible<tag>",
                }
            ],
        }
        html = render_junction_report(result)
        parsed = self.assert_offline_accessible(html)

        # Raw payloads must never appear unescaped.
        self.assertNotIn('<script>alert("o")</script>', html)
        self.assertNotIn("<img src=x onerror=1>", html)
        self.assertNotIn("<b>check</b>", html)
        self.assertNotIn('L<bad>&"', html)

        # Overview truncates long alleles; details retains full safely escaped sequences.
        overview = " ".join(parsed.overview_text)
        details = " ".join(parsed.details_text)
        self.assertNotIn(long_orig, overview)
        self.assertNotIn(long_prop, overview)
        self.assertIn("\u2026", overview)
        self.assertIn(long_orig, details)
        self.assertIn(long_prop, details)
        self.assertIn('Warn <b>check</b> & "review"', overview)
        self.assertIn('L<bad>&"/v0', details)

        # Only explicit safe relative artifact links are emitted; no raw paths linked.
        hrefs = {attrs["href"] for tag, attrs in parsed.tags if tag == "a" and "href" in attrs}
        self.assertEqual(hrefs, {"evidence/junction-summary.json"})
        self.assertNotIn("/abs/secret/composition.json", html)

    def test_preserved_vs_repaired_and_conflict_scientific_honesty(self):
        preserved_html = render_junction_report({
            "junction": {
                "mode": "gap",
                "gap_length": 250,
                "flank_length_delta": 0,
                "overlap_trim": 0,
                "flanks": {
                    "L": {"status": "preserved", "proposed_variants": 0, "accepted_variants": 0, "protected_variants": 0},
                    "R": {"status": "preserved", "proposed_variants": 0, "accepted_variants": 0, "protected_variants": 0},
                },
                "repairs": [],
                "warnings": [],
            }
        })
        self.assert_offline_accessible(preserved_html)
        self.assertIn("Preserved", preserved_html)
        self.assertNotIn("Conflict", preserved_html)

        mixed_html = render_junction_report({
            "junction": {
                "mode": "gap",
                "gap_length": 420,
                "flank_length_delta": -7,
                "overlap_trim": 0,
                "flanks": {
                    "L": {"status": "repaired", "proposed_variants": 2, "accepted_variants": 1, "protected_variants": 1},
                    "R": {"status": "conflict", "proposed_variants": 1, "accepted_variants": 0, "protected_variants": 0},
                },
                "repairs": [
                    {
                        "side": "L",
                        "start0": 19900,
                        "end0": 19907,
                        "original": "ACGTACG",
                        "proposed": "",
                        "status": "accepted",
                        "reason": "independent_spanning_reads",
                        "alt_reads": 8,
                        "ref_reads": 1,
                        "other_reads": 0,
                        "spanning_reads": 9,
                        "alt_fraction": 8 / 9,
                        "strands": ["+", "-"],
                    },
                    {
                        "side": "L",
                        "start0": 500,
                        "end0": 501,
                        "original": "A",
                        "proposed": "G",
                        "status": "protected",
                        "reason": "outside_bounded_terminal_window",
                    },
                    {
                        "side": "R",
                        "start0": 120,
                        "end0": 121,
                        "original": "C",
                        "proposed": "T",
                        "status": "retained",
                        "reason": "insufficient_or_conflicting_read_support",
                        "alt_reads": 2,
                        "ref_reads": 4,
                        "other_reads": 1,
                        "spanning_reads": 7,
                        "alt_fraction": 2 / 7,
                        "strands": ["+"],
                    },
                ],
                "warnings": ["Unconfirmed end differences retained as original sequence; not proven corrections"],
            }
        })
        parsed = self.assert_offline_accessible(mixed_html)
        text = parsed.all_text

        # Gap-only length and flank delta must stay independent (never combined as 413 bp).
        self.assertIn("420 bp", text)
        self.assertIn("-7 bp", text)
        self.assertNotIn("413 bp", text)

        # Scientific honesty: retained is NOT applied; conflict is warned; edits are not proven corrections.
        self.assertIn("Retained (not applied)", text)
        self.assertIn("Conflict (original retained)", text)
        self.assertIn("Protected (not applied)", text)
        self.assertIn("Repaired (read-backed)", text)
        self.assertIn("not proven biological corrections", text)
        self.assertIn("badge-warn", self.rendered_classes(parsed))
        self.assertIn("\u2205 (0 bp)", text)

    def test_true_overlap_trim_is_labeled_as_deliberate_deduplication(self):
        html = render_junction_report({
            "junction": {
                "mode": "overlap",
                "gap_length": -185,
                "flank_length_delta": 0,
                "overlap_trim": 185,
                "flanks": {
                    "L": {"status": "preserved", "proposed_variants": 0, "accepted_variants": 0, "protected_variants": 0},
                    "R": {"status": "preserved", "proposed_variants": 0, "accepted_variants": 0, "protected_variants": 0},
                },
                "repairs": [],
                "warnings": ["Original overlap deduplicated; left original copy retained, not polished"],
            }
        })
        parsed = self.assert_offline_accessible(html)
        text = parsed.all_text
        self.assertIn("-185 bp", text)
        self.assertIn("185 bp", text)
        self.assertIn("0 bp", text)
        self.assertIn("deliberate deduplication, not silent polishing", text.lower())
        self.assertIn("Original overlap", text)

    def test_diagnostics_fast_full_and_rejections_without_inventing_success(self):
        html = render_junction_report({
            "connection_diagnostics": [
                {
                    "left_family": "L",
                    "left_version": 0,
                    "right_family": "R",
                    "right_version": 0,
                    "stages": [
                        {"stage": "fast", "blocks": 0, "candidates": 0, "reason": "no_eligible_terminal_chain"},
                    ],
                    "status": "rejected",
                    "reason": "no_eligible_terminal_chain",
                },
                {
                    "left_family": "L",
                    "left_version": 1,
                    "right_family": "R",
                    "right_version": 0,
                    "stages": [
                        {"stage": "fast", "blocks": 1, "candidates": 0, "reason": "no_eligible_terminal_chain"},
                        {"stage": "full", "blocks": 2, "candidates": 0, "reason": "ambiguous_or_unmapped_original_L_boundary"},
                    ],
                    "status": "rejected",
                    "reason": "ambiguous_or_unmapped_original_L_boundary",
                },
            ]
        })
        parsed = self.assert_offline_accessible(html)
        text = parsed.all_text
        self.assertIn("L/v0", text)
        self.assertIn("R/v0", text)
        self.assertIn("no_eligible_terminal_chain", text)
        self.assertIn("ambiguous_or_unmapped_original_L_boundary", text)
        self.assertIn("No eligible connection", text)
        # All pairs rejected -> no rendered element may carry badge-ok, and no junction metrics invented.
        classes = self.rendered_classes(parsed)
        self.assertNotIn("badge-ok", classes)
        self.assertIn("badge-warn", classes)
        self.assertNotIn("Gap-only length", text)

    def test_missing_fields_and_partial_schema_tolerance(self):
        html = render_junction_report({
            "junction": {
                "mode": "gap",
                "gap_length": 0,
                "flank_length_delta": True,  # bool must not be misformatted as 1 bp
                "repairs": [
                    "malformed-entry",
                    {
                        "side": "R",
                        "start": 10,
                        "end": 12,
                        "original": "",
                        "proposed": "AT",
                        "status": "retained",
                        "strands": [],
                    },
                ],
            },
            "connection_diagnostics": [
                "not-a-dict",
                {"left_family": "primary", "right_family": "candidate0", "status": "rejected"},
            ],
        })
        parsed = self.assert_offline_accessible(html)
        text = parsed.all_text
        self.assertIn("0 bp", text)
        self.assertNotIn("+1 bp", text)
        self.assertIn("10..12", text)
        self.assertIn("none / 无", text)

    def test_project_junctions_aggregates_events_history_and_direct_list(self):
        events = [
            {
                "number": 1,
                "primary": "ctg1",
                "direction": "right",
                "status": "closed",
                "junction": {
                    "mode": "gap",
                    "gap_length": 600,
                    "flank_length_delta": 3,
                    "overlap_trim": 0,
                    "flanks": {
                        "L": {"status": "repaired", "proposed_variants": 1, "accepted_variants": 1, "protected_variants": 0},
                        "R": {"status": "preserved", "proposed_variants": 0, "accepted_variants": 0, "protected_variants": 0},
                    },
                    "repairs": [
                        {
                            "side": "L",
                            "start0": 9900,
                            "end0": 9900,
                            "original": "",
                            "proposed": "ATG",
                            "status": "accepted",
                            "reason": "independent_spanning_reads",
                            "alt_reads": 6,
                            "ref_reads": 0,
                            "other_reads": 0,
                            "spanning_reads": 6,
                            "alt_fraction": 1.0,
                            "strands": ["+", "-"],
                        }
                    ],
                    "warnings": [],
                },
                "connection_diagnostics": [
                    {
                        "left_family": "primary",
                        "left_version": 1,
                        "right_family": "candidate0",
                        "right_version": 0,
                        "stages": [
                            {"stage": "fast", "blocks": 0, "candidates": 0, "reason": "no_eligible_terminal_chain"},
                            {"stage": "full", "blocks": 2, "candidates": 1, "reason": "eligible_connection"},
                        ],
                        "status": "eligible",
                        "reason": "eligible_connection",
                    }
                ],
            },
            {
                "number": 2,
                "primary": "DG1",
                "direction": "left",
                "status": "closed",
                "junction": {
                    "mode": "overlap",
                    "gap_length": -90,
                    "flank_length_delta": 0,
                    "overlap_trim": 90,
                    "flanks": {
                        "L": {"status": "preserved", "proposed_variants": 0, "accepted_variants": 0, "protected_variants": 0},
                        "R": {"status": "preserved", "proposed_variants": 0, "accepted_variants": 0, "protected_variants": 0},
                    },
                    "repairs": [],
                    "warnings": ["Original overlap deduplicated; left original copy retained, not polished"],
                },
            },
            {
                "number": 3,
                "primary": "DG2",
                "direction": "right",
                "status": "no_candidates",
            },
        ]

        from_events = render_project_junctions({"events": events})
        from_history = render_project_junctions({"history": events})
        from_list = render_project_junctions(events)

        for html in (from_events, from_history, from_list):
            parsed = self.assert_offline_accessible(html)
            text = parsed.all_text
            self.assertIn("2 / 3", text)
            self.assertIn("#1 ctg1 (right)", text)
            self.assertIn("#2 DG1 (left)", text)
            self.assertNotIn("#3 DG2", text)
            self.assertIn("600 bp", text)
            self.assertIn("+3 bp", text)
            self.assertIn("-90 bp", text)
            self.assertIn("90 bp", text)
            # Must not invent combined 510 bp sum of +600 gap and -90 overlap.
            self.assertNotIn("510 bp", text)

    def test_anchor_retries_shared_budget_block_and_escaping(self):
        html = render_junction_report({
            "anchor_retries": [
                "not-a-dict",
                {
                    "family": 'L<1>&"',
                    "attempt": 2,
                    "backoff_bases": 2000,
                    "reason": "terminal_anchor_backoff",
                },
                {
                    "family": "R",
                    "attempt": 0,
                    "backoff_bases": True,
                },
            ]
        })
        parsed = self.assert_offline_accessible(html)
        text = parsed.all_text
        classes = self.rendered_classes(parsed)

        self.assertNotIn('L<1>&"', html)
        self.assertIn('L<1>&"', text)
        self.assertIn("2,000 bp", text)
        self.assertIn("terminal_anchor_backoff", text)
        self.assertIn("shared round budget", text.lower())
        self.assertIn("not free hidden retries", text.lower())
        self.assertIn("not proven corrections", text.lower())
        self.assertIn("共享轮次预算", text)
        self.assertNotIn("+1 bp", text)
        self.assertNotIn("badge-ok", classes)
        self.assertIn("badge-warn", classes)
        self.assertNotIn("Gap-only length", text)
        self.assertNotIn("Checked pairs", text)

    def test_project_junctions_aggregates_anchor_retries_on_events(self):
        events = [
            {
                "number": 1,
                "primary": "ctgA",
                "direction": "right",
                "status": "closed",
                "junction": {
                    "mode": "gap",
                    "gap_length": 400,
                    "flank_length_delta": 0,
                    "overlap_trim": 0,
                    "flanks": {
                        "L": {"status": "preserved", "proposed_variants": 0, "accepted_variants": 0, "protected_variants": 0},
                        "R": {"status": "preserved", "proposed_variants": 0, "accepted_variants": 0, "protected_variants": 0},
                    },
                    "repairs": [],
                    "warnings": [],
                },
                "anchor_retries": [
                    {"family": "primary", "attempt": 2, "backoff_bases": 2000, "reason": "terminal_anchor_backoff"},
                ],
            },
            {
                "number": 2,
                "primary": "ctgB",
                "direction": "left",
                "status": "round_limit",
                "anchor_retries": [
                    {"family": "primary", "attempt": 3, "backoff_bases": 2000, "reason": "terminal_anchor_backoff"},
                    {"family": "candidate0", "attempt": 4, "backoff_bases": 2000, "reason": "terminal_anchor_backoff"},
                ],
            },
            {
                "number": 3,
                "primary": "ctgC",
                "direction": "right",
                "status": "no_candidates",
            },
        ]
        html = render_project_junctions({"events": events})
        parsed = self.assert_offline_accessible(html)
        text = parsed.all_text

        self.assertIn("2 / 3", text)
        self.assertIn("#1 ctgA (right)", text)
        self.assertIn("#2 ctgB (left)", text)
        self.assertNotIn("#3 ctgC", text)
        self.assertIn("3 (2 evt)", text)
        self.assertIn("Anchor retries / 末端回退", text)
        self.assertIn("terminal_anchor_backoff", text)
        self.assertIn("2,000 bp", text)
        self.assertIn("shared round budget", text.lower())
        self.assertIn("not free hidden retries", text.lower())
        self.assertIn("not proven corrections", text.lower())


if __name__ == "__main__":
    unittest.main()
