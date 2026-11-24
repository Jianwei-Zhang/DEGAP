"""Offline bidirectional report, using real-schema evidence without sequence I/O."""

import copy
from html.parser import HTMLParser
from hashlib import sha256
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from FamilyResult import render_family_report_from_data as render
from ExtensionFamilies import fingerprint
from FamilyResult import SCHEMA
from render_family_report import write_preview


class Page(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.tags = []
        self.script = ""
        self.in_script = False
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))
        if tag == "script":
            self.in_script = True

    def handle_endtag(self, tag):
        if tag == "script":
            self.in_script = False

    def handle_data(self, value):
        if self.in_script:
            self.script += value

    def attributes(self, tag):
        return [a for t, a in self.tags if t == tag]


def fixture(single=False, closed=True):
    counts = [1, 0] if single else [2, 1]
    families = [
        dict(
            key=s,
            direction=direction,
            attempts=n,
            added_bases=100 * n,
            stop_reason=None,
            versions=[
                dict(number=i, length=1000 + i * 100, sha256=str(i) * 64, path="/private/reads/sequence.fa")
                for i in range(n + 1)
            ],
        )
        for s, direction, n in zip("LR", ("rightward", "leftward"), counts)
    ]
    a, b = counts
    w = dict(
        left_family="L",
        right_family="R",
        left_version=a,
        right_version=b,
        left_start=900,
        left_end=1000 + 100 * a,
        right_start=2,
        right_end=130,
        overlap=127,
        identity=99.2,
        overhang=2,
        evidence={"stage": "full", "plan_path": "/private/plan"},
    )
    diagnostics = [
        dict(
            left_family="L",
            right_family="R",
            left_version=i,
            right_version=k,
            status="eligible" if closed and (i, k) == (a, b) else "rejected",
            reason="eligible_connection" if closed and (i, k) == (a, b) else "overlap_below_minimum",
            stages=[
                {
                    "stage": "fast",
                    "blocks": 2,
                    "rejections": {"blocks": 2, "short": 2},
                    "reason": "overlap_below_minimum",
                }
            ],
        )
        for i in range(a + 1)
        for k in range(b + 1)
    ]
    return dict(
        families=families,
        winner=w if closed else None,
        status="closed" if closed else "exhausted",
        gap_length=20 if closed else None,
        total_rounds=a + b,
        final_sequence="sample.final.fa",
        search_contract="fixture",
        connection_diagnostics=diagnostics,
        evidence_artifacts={"unsafe": "other.fa"},
    )


class FamilyPlayerTests(unittest.TestCase):
    def test_shared_round_chart_counts_each_gain_once(self):
        html = render(fixture())
        for value in (
            "Round 1 · L +100 bp cumulative",
            "Round 2 · R +100 bp cumulative",
            "Round 3 · L +200 bp cumulative",
            "L/v2 · 1,200 bp",
            "R/v1 · 1,100 bp",
        ):
            self.assertIn(value, html)
        self.assertNotIn("Round 2 · L +200 bp cumulative", html)

    def test_historical_winner_is_not_replaced_by_latest_peer(self):
        data = fixture()
        data["winner"]["right_version"] = 0
        data["connection_diagnostics"] = []
        html = render(data)
        self.assertIn("L/v2 × R/v0: Winner", html)
        self.assertIn("R/v1 · 1,100 bp", html)  # live round state
        self.assertIn("R/v0</text>", html)  # actual winning pair in full-sequence alignment
        self.assertIn("No diagnostic recorded", html)

    def test_signed_legacy_overlap_and_verified_right_flank_trim(self):
        data = fixture()
        data["gap_length"] = -12
        self.assertIn("Net length vs original flanks <b>-12 bp", render(data))
        data["junction"] = {"mode": "overlap", "gap_length": -12, "overlap_trim": 12}
        html = render(data)
        for value in ("Net gap-only fill <b>0 bp", "Right flank 988 bp", "original overlap trimmed: 12 bp"):
            self.assertIn(value, html)

    def test_preview_reads_only_json_and_never_overwrites_publication(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = fixture()
            raw = json.dumps(data["connection_diagnostics"]).encode()
            data.update(schema=SCHEMA, artifacts={"connection-diagnostics.json": sha256(raw).hexdigest()})
            data["payload_sha256"] = fingerprint(data)
            result = root / "family-result.json"
            result.write_text(json.dumps(data))
            original = root / "family-report.html"
            original.write_text("original")
            side = root / "connection-diagnostics.json"
            side.write_text("[]")
            with self.assertRaisesRegex(ValueError, "artifact changed"):
                write_preview(root)
            self.assertFalse((root / "family-report.preview.html").exists())
            side.write_bytes(raw)
            preview = write_preview(root)
            self.assertIn("L/v2 × R/v1: Winner", preview.read_text())
            with self.assertRaises(FileExistsError):
                write_preview(root)
            self.assertEqual(original.read_text(), "original")
            self.assertEqual(json.loads(result.read_text()), data)
            data["total_rounds"] = 42
            result.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, "Invalid family result"):
                write_preview(root)

    def test_stale_checkpoint_or_contradictory_outcome_cannot_supply_evidence(self):
        data = fixture()
        checkpoint = {"contract_id": data["search_contract"], "families": copy.deepcopy(data["families"])}
        checkpoint["families"][0]["versions"][1]["sha256"] = "changed"
        with self.assertRaisesRegex(ValueError, "history"):
            render(data, checkpoint)
        with self.assertRaisesRegex(ValueError, "outcome"):
            render(data, {"winner": None})
        with self.assertRaisesRegex(ValueError, "Closure status"):
            render(dict(data, status="exhausted"))

    def test_closed_multiround_final_view_exact_pair_intervals_and_rounds(self):
        data = fixture()
        before = copy.deepcopy(data)
        html = render(data)
        page = Page(html)
        self.assertEqual(data, before)
        self.assertNotIn("/private/", html)
        self.assertNotIn("other.fa", html)
        self.assertIn("L [900, 1,200) ↔ R [2, 130)", html)
        self.assertIn("Stage: full", html)
        self.assertIn("L→v2 · 闭合", html)
        self.assertIn("L/v2 × R/v1: Winner", html)
        self.assertIn("fast: 2 blocks · overlap_below_minimum · short: 2", html)
        self.assertEqual(
            [a["data-caption"] for a in page.attributes("div") if "data-caption" in a and "hidden" not in a],
            ["3"],
        )
        self.assertEqual(
            [a["data-step"] for a in page.attributes("button") if a.get("aria-current") == "step"], ["3"]
        )
        self.assertTrue(all("open" not in a for a in page.attributes("details")))
        self.assertEqual(html.count('class="overlap-bridge"'), 2)

    def test_one_round_and_zero_round_closures(self):
        data = fixture(single=True)
        html = render(data)
        self.assertIn("L→v1 · 闭合", html)
        self.assertNotIn("R→v1", html)
        self.assertIn('data-from="1"><path class="overlap-bridge"', html)
        data["families"][0].update(attempts=0, added_bases=0, versions=data["families"][0]["versions"][:1])
        data["winner"].update(left_version=0, left_end=1000)
        data.update(total_rounds=0, connection_diagnostics=[])
        page = Page(render(data))
        self.assertTrue(all("disabled" in a for a in page.attributes("button") if "class" in a))

    def test_exhausted_replays_backoff_without_inventing_growth_or_overlap(self):
        data = fixture(closed=False)
        data["families"][0].update(attempts=4, stop_reason="no_extension_contigs")
        data["anchor_backoff_requests"] = [dict(family="L", after_attempt=1)]
        data["total_rounds"] = 5
        html = render(data)
        self.assertIn("Round 1 · L/v0 · +0 bp", html)
        self.assertIn("Round 2 · R/v1 · +100 bp", html)
        self.assertIn("Round 5 · L/v2 · +0 bp", html)
        self.assertIn("no_extension_contigs", html)
        self.assertIn("url(#hatch)", html)
        self.assertNotIn('class="overlap-bridge"', html)
        self.assertNotIn("Final junction composition", html)
        self.assertIn("No eligible connection", html)

    def test_conflict_and_protection_do_not_change_fill_or_flank_length(self):
        data = fixture()
        summary = dict(
            mode="gap",
            gap_length=20,
            overlap_trim=0,
            flank_length_delta=0,
            flanks={
                "R": dict(status="conflict", accepted_variants=0, proposed_variants=1, protected_variants=1)
            },
            repairs=[dict(side="R", start0=4, end0=5, original="A", proposed="T", status="protected")],
            evidence_artifacts={"malicious": "other.html"},
            policy={"repair_window": 100},
        )
        html = render(
            data,
            junction_evidence={"summary": summary},
            connection_diagnostics=data["connection_diagnostics"],
        )
        self.assertIn("R [4, 5) · protected", html)
        self.assertIn("R: conflict (0 accepted / 1 proposed)", html)
        self.assertIn("Right flank 1,000 bp", html)
        self.assertIn("Net fill 20 bp", html)
        self.assertNotIn("other.html", html)
        data["junction"] = dict(summary, gap_length=100)
        with self.assertRaisesRegex(ValueError, "Conflicting junction"):
            render(data, junction_evidence={"summary": summary})
        with self.assertRaisesRegex(ValueError, "different family"):
            render(data, {"contract_id": "wrong"})

    @unittest.skipUnless(shutil.which("node"), "Node required for offline player control test")
    def test_js_syntax_and_controls(self):
        script = Page(render(fixture())).script
        subprocess.run(["node", "--check"], input=script, text=True, check=True)
        # Small DOM surface, real event handlers, deterministic timer; no browser/network needed.
        harness = """
const assert = require('node:assert/strict');
class El {constructor(){this.dataset={};this.attrs={};} setAttribute(k,v){this.attrs[k]=v;}
  toggleAttribute(k,v){this[k]=v;} addEventListener(k,f){this[k]=f;}}
const els=()=>Array.from({length:4},()=>new El());
const pills=els(),caps=els(),layers=els(),buttons={'.play':new El(),'.prev':new El(),'.next':new El()};
layers.forEach((e,i)=>e.dataset.from=i);
const root={querySelector:s=>buttons[s],querySelectorAll:s=>({'[data-step]':pills,'[data-from]':layers,'[data-caption]':caps})[s]};
const document={querySelector:()=>root}; let tick=null;
const setInterval=f=>{tick=f;return 1;}, clearInterval=()=>{tick=null;};
"""
        assertions = """
assert.equal(caps[3].hidden,false); buttons['.play'].click(); assert.equal(caps[0].hidden,false);
assert.equal(layers[3].hidden,true); tick(); assert.equal(caps[1].hidden,false);
buttons['.play'].click(); assert.equal(tick,null); pills[2].click(); assert.equal(caps[2].hidden,false);
buttons['.prev'].click(); assert.equal(caps[1].hidden,false); buttons['.next'].click();
buttons['.play'].click();tick();assert.equal(caps[3].hidden,false);assert.equal(tick,null);
"""
        subprocess.run(["node"], input=harness + script + assertions, text=True, check=True)

    @unittest.skipUnless(
        os.getenv("CGAT_NATIVE_REPORT_INSPECTOR"),
        "Set CGAT_NATIVE_REPORT_INSPECTOR to test the actual CGAT inspector",
    )
    def test_real_cgat_transform(self):
        spec = importlib.util.spec_from_file_location(
            "cgat_offline_view", os.environ["CGAT_NATIVE_REPORT_INSPECTOR"]
        )
        inspector = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = inspector
        spec.loader.exec_module(inspector)
        for closed in (True, False):
            raw = render(fixture(closed=closed)).encode()
            page, stats = inspector.transform(raw, (), (), artifact_job="sample")
            self.assertIn(b"Content-Security-Policy", page)
            self.assertEqual(stats["scripts"], 1)
            self.assertEqual(stats["detached_native_artifact_links"], 6)
