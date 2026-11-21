# Junction checks and bounded terminal repair: SSH 180 validation

Date: **2026-11-22**. Branch: `feat/bidirectional-gap-extension`.

## What changed

- A failed cropped-end comparison now falls back to full sequences. Small, compatible collinear blocks can form a connection; full containment is recognized, including inputs already shorter than the fast window. Effective overlap requirements are not lowered just to close a case.
- Original flank cores are protected. The gap is cut using mapped original boundaries, not final length minus original flank lengths. Bounded terminal candidates come from both assembly differences and **original-read consensus**, including errors inherited unchanged by the seed and cumulative copies.
- Each applied terminal edit needs distinct spanning reads, both strands, an interior anchor and the configured allele fraction. Unconfirmed proposals retain the original allele. Input snapshots are never overwritten.
- Gap-only length, terminal-edit length delta and true-overlap trimming are separate. AutoGapfiller integrates only the gap and approved end changes; CtgLinker retains actual source-component provenance and reconstructable AGP.
- Anchor backoff is one separately budgeted attempt. A request that never executes is not reported as an executed retry or charged another round.

The default policy and limitations are in [junction-validation.md](../design/junction-validation.md). These defaults are not claimed to be empirically optimal polishing thresholds.

## Validation design and provenance

The **22 fixed-coordinate workloads** reuse exactly the inputs, masks, read files, 12-attempt / 300 kb shared budgets and tool versions from [the expanded matrix](2026-11-22-ssh180-expanded-matrix.md). They cover nine ordinary HiFi gaps, nine high-depth Klebsiella ONT gaps, both AutoGapfiller datasets, CtgLinker and one low-depth ONT negative control. No gaps were relocated or thresholds changed to improve the reported outcomes.

Four representative consumer runs were repeated after fixing the input-hash guard. Final read-consensus code then rechecked **29 closed junctions from those immutable completed assemblies** and regenerated **36 family reports/logs**, without rerunning their local assembly. Candidate geometry and gap-only sequence were checked against the original evidence; newly supported end edits were evaluated separately. These regenerated workspaces are explicitly marked derivative validation reports, not migrated/resumable old searches.

The final policy was also exercised in fresh CtgLinker and AutoGapfiller runs. The AutoGapfiller input deliberately contains one wrong flank base, described below. Runtime source manifests distinguish the intermediate and final checks; original checkpoints/results remain unchanged.

## High-depth ONT: fixed 9-gap matrix

Dataset: ENA `ERR12342517`, Klebsiella MGH 78578; comparator: `CP000647.1`. The unchanged input has about 70x MAPQ-filtered aligned-span depth, with about 42–76x across the fixed gap regions. Span depth is not effective assembly depth.

Closed gaps increase from **1/9 to 7/9** under the updated search/composition policy. This is one isolate, three positions and three gap lengths—not a general ONT success-rate estimate.

| Case | Status | Shared attempts | Selected versions | Gap-only length (bp) | Gap-only identity |
| --- | --- | ---: | --- | ---: | ---: |
| `klebs-mid-20k` | `exhausted` | 4 | — | — | — |
| `klebs-mid-50k` | `closed` | 2 | L0 + R1 | 49958 | 98.80703% |
| `klebs-mid-5k` | `closed` | 2 | L0 + R1 | 4999 | 99.96000% |
| `klebs-q1-20k` | `closed` | 2 | L0 + R1 | 20000 | 99.94501% |
| `klebs-q1-50k` | `closed` | 5 | L0 + R3 | 49996 | 99.85407% |
| `klebs-q1-5k` | `closed` | 2 | L0 + R1 | 5000 | 100.00000% |
| `klebs-q3-20k` | `closed` | 2 | L0 + R1 | 20000 | 99.96000% |
| `klebs-q3-50k` | `exhausted` | 4 | — | — | — |
| `klebs-q3-5k` | `closed` | 2 | L0 + R1 | 5000 | 99.98000% |

All seven closed gaps cover the complete comparator gap. `klebs-mid-20k` and `klebs-q3-50k` remain unclosed and retain separate family outputs. Their failures are not removed from the denominator.

**Quality warning:** `klebs-mid-50k` closes, but its 49,958 bp gap has only **98.80703%** identity to the 50 kb comparator gap. Closure is not proof of correct sequence. Differences may include assembly error or sample/reference differences; this test does not distinguish them conclusively.

The final read-consensus policy additionally inserts one `G` into the left original end of that case, at original flank boundary `[59794,59794)`: **51 alternate-supporting reads out of 52 spanning reads**, both strands (98.08%). This is an explicit bounded edit, not silent consensus replacement. Gap-only sequence remains unchanged at 49,958 bp; flank length delta is **+1 bp**. Read support is evidence, not independent proof that the comparator is wrong. No other ordinary matrix output changes under the final read-consensus pass.

## HiFi and consumer modes

- Ordinary HiFi: **9/9 closed**. Histories, global single-winner ranking, shared budgets, artifacts and resume-without-recomputation checks pass.
- HiFi AutoGapfiller: **3/3 filled**, no remaining N. Comparator non-gap regions are unchanged in the ordinary run.
- Klebsiella ONT AutoGapfiller: **1/3 filled**; 70,000 N bases remain. Its fixed gaps are q1 5 kb, mid 20 kb and q3 50 kb, so the nine-case improvement does not imply this particular three-gap input also improves.
- CtgLinker: three 100 kb contigs, two 20 kb gaps and a reversed middle input yield **one 339,995 bp scaffold**, members `+ / - / +`. Both gap regions cover the comparator completely; gap-only identities are 99.58517% and 99.97500%. AGP reconstructs the sequence exactly.
- Low-depth Salmonella ONT negative control remains unclosed. Its normal stop is not mislabeled as a tool error or closure.

## Known end errors, genuine reads

At a fixed site 100 bp before the existing HiFi mid-5 kb gap, the same artificial edit is applied to **both the initial flank and the cumulative extension copy**. The original source files and original reads are unchanged. This explicitly tests an inherited seed error, not just an already-correct assembly overwriting a bad flank.

| Controlled edit | Final decision | Spanning-read evidence | Gap-only length | Flank length delta |
| --- | --- | --- | ---: | ---: |
| Substitution | Restore baseline allele | 25 alternate, 0 reference, both strands | 5,000 bp | 0 bp |
| Extra base | Remove extra base | 25 alternate, 0 reference, both strands | 5,000 bp | -1 bp |
| Missing base | Restore missing base | 25 alternate, 0 reference, both strands | 5,000 bp | +1 bp |
| Unsupported assembled substitution | Keep original allele | 0 alternate, 25 reference | 5,000 bp | 0 bp |

A fresh AutoGapfiller run additionally starts from a genome with the mid-5 kb gap and a wrong `C` at original zero-based position **1,105,971**. It closes the gap, restores `T` from **25 supporting reads**, and integrates the corrected end without changing the reported 5,000 bp gap length. Its proposal is recorded as `read_consensus`; the normal extension had inherited the wrong base. The final genome has no remaining N. This validates consumer integration, not a population-level natural-error correction rate.

## Automated checks and integrity

- Final scientific runtime: **214 tests passed, 0 skipped**, 144.481 s, Donau job `43010294`. Tests include full/cropped/contained overlaps, fragmented chains, true-overlap deduplication, both inherited end errors, indel boundary projection, false/ambiguous read support, shared budget/recovery and exact FASTA/AGP reconstruction.
- The initial remote suite exposed an existing stat-only hash-cache hole on a coarse-timestamp filesystem. Small/recent files now use a fresh hash; both quality-only FASTQ mutation tests pass. The failed diagnostic log is preserved rather than hidden.
- Composition cache identities now include the complete participating version content, not just version numbers/coordinates. Modified completed receipts/components are rejected.
- Scientific matrix audit: **22 workloads, 28 family searches, 37 HTML files and 238 local links**. Representative audit: **4 workloads, 8 searches, 13 HTML files and 84 links**. Histories, ranking, actual HTML rows, input/source hashes, offline resources and integration reconstruction pass.
- HTML evidence was implemented by the requested AGY model **`gemini-4-neon-high`**. Final tables separate gap length, flank delta, original overlap trimming, applied/unapplied variants, fast/full rejection reasons and executed anchor retries. No browser validation or new dependency installation is claimed.

Portable machine-readable results: [2026-11-22-junction-repair-summary.json](2026-11-22-junction-repair-summary.json).

## Final HTML coordinates and unified checks

AGY `gemini-4-neon-high` also corrected AutoGapfiller's JavaScript coordinate projection: original gap ends are inclusive, downstream shifts include the separate flank length delta, zero/negative gap lengths are handled, and gaps are sorted by genomic position. Chromosome tracks consistently use the original coordinate frame. An AST/source audit confirms that this last follow-up changes only `_create_full_html` and its tests; the scientific runtime is byte-identical to the final read-consensus snapshot.

- Exact final branch source: **220 tests run remotely, 5 skipped** (Node unavailable), 152.630 s, job `43010307`. The complete 214-test scientific suite had already passed remotely with no skips. All **6 coordinate tests pass locally with Node**, including the five skipped remotely. The final local unified run is 220 tests, 12 skipped because local hifiasm-dependent fixtures are unavailable; those fixtures pass remotely. No single environment is claimed to have run every test without skips.
- The inline JavaScript from the actual final HiFi controlled-error and ONT AutoGapfiller reports passes Node syntax/VM checks. Projected remaining-gap coordinates exactly match N intervals read from their final FASTA files: no interval for the closed HiFi case; ONT intervals **2,647,561–2,667,560** and **3,961,341–4,011,340** (1-based inclusive). These are lightweight report checks, not browser acceptance.
- The final portable archive verifies **99 HTML files and 617 local links**. After download, all **2,615 manifest-listed files** match their SHA-256 values; `README.txt` and the manifest itself are archive-envelope files. Raw reads, BAMs and assembly caches are excluded.

Evidence archive on SSH `180`:

```text
/share/home/xbzhang26/cgat-reads-debug/degap-junction-report-final-20261122/degap-junction-evidence-20261122.tar.gz
SHA-256: 2b519455ce4b77adeeae45c64adebdf7c0aaac9c578d9b6658b36693a246ccd7
Size: 92,709,399 bytes
```

`matrix/` preserves the original 22 workloads, `representatives/` the four reruns, `delivery/` the final read-consensus checks and fresh consumers, and `final-report/` the exact final renderer and report views. Derivative reports remain explicitly marked; they are not relabeled as fresh searches. Source fingerprints and detailed test/provenance records are in the adjacent machine-readable summary.
