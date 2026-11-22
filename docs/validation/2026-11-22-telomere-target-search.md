# TelSeeker shared target search: validation on SSH 180

Date: **2026-11-22**. Branch: `feat/bidirectional-gap-extension`.

## Implemented behavior

Each chromosome end keeps one chromosome-side extension family. Tel reads are static targets, not globally clustered groups or independently extended families. Endpoints share an immutable, content-bound SQLite index; retrieved candidates are refined in batches of at most 256 records or 8 MiB. These are batch bounds, not a top-N candidate limit.

The best supported eligible candidate wins after all retrieved records have been refined. At least two distinct read IDs must corroborate the nonrepeat anchor and up to 200 outward junction bases. Different orientations of one ID are not independent evidence. At the final join, the current chromosome-side sequence is retained and only the selected read's outside tail is added. The existing single-family round/growth budgets remain in force.

See [the implementation contract](../design/telomere-target-search.md) for exact retrieval, alignment, support and cache rules.

## Genuine-read controls: fixed internal coordinates

These are **target-search/integration controls, not natural telomere-completion benchmarks**. Each dataset uses one fixed internal comparator midpoint, checked in L and R directions. Reads are selected from the beginning of the original file, without moving the coordinate after a negative result. Each read is stored forward and reverse-complement under the same ID. Original source and comparator SHA-256 values are recorded in the adjacent JSON.

| Dataset / original IDs | Target records | Direction | Retrieved / eligible | Decision | Distinct-ID support lower bound | Selected tail |
| --- | ---: | --- | ---: | --- | ---: | ---: |
| HiFi / 3,000 | 6,000 | L | 3,718 / 10 | Connected | 2 | 3,429 bp |
| HiFi / 3,000 | 6,000 | R | 4,558 / 11 | Connected | 2 | 4,850 bp |
| ONT / 5,000 | 10,000 | L | 39 / 1 | Insufficient read support | 1 | — |
| ONT / 5,000 | 10,000 | R | 41 / 1 | Insufficient read support | 1 | — |
| ONT / 20,000 | 40,000 | L | 129 / 3 | Connected | 2 | 4,690 bp |
| ONT / 20,000 | 40,000 | R | 171 / 1 | Insufficient read support | 1 | — |

The larger ONT subset contains the smaller one; these are **not independent loci or replicates**. The original negatives are retained, not replaced by the later positive result. The policy does not turn one qualifying read into a supported connection, regardless of total pool size.

First index-plus-two-direction search times were 140.154 s (HiFi), 53.115 s (small ONT), and 161.959 s (larger ONT). Identical-input replay took 0.841 s, 0.295 s, and 0.048 s respectively, with no index or raw-alignment recomputation. These are observed runs, not a standardized cross-dataset speed comparison. Repetitive HiFi still retrieves many coarse candidates; batching does not remove that cost.

Supporting reads corroborate the **local junction**, not the entire tail or telomeric repeat count. Comparator mappings of complete joined sequences also include 30 kb of retained chromosome and are not tail-only accuracy measurements. The larger ONT L join has 133 outer query bases clipped/unmapped by the comparator. No whole-tail polishing or full-tail correctness is claimed.

## Scale and public-CLI controls

- A synthetic pool of **20,002 records** retrieves and refines exactly two relevant records in one batch, selecting the stable best supported candidate. A pure `TTAGGG` query retrieves none and launches no alignment batch.
- Sequential posting writes followed by one covering-index build replace per-seed random B-tree writes. All **601,232 logical postings** exactly match the previous construction, SHA-256 `6c4705295a30180e112b025451b12a856d97ee1ca3f3d85701be4ae762c0ac13`. Selection thresholds were not weakened. This run built in 10.212 s, searched in 1.739 s, and replayed in 0.074 s; Python peak RSS was 106.10 MiB.
- The synthetic public `DEGAP.py --mode telseeker` pipeline runs Part1 through Part3 and reports with `--MaximumExtensionRound 0`: two supported direct joins, **zero extension attempts**, 8 kb added at each end, and a final **120,000 bp sequence exactly matching truth** from a 104,000 bp input. This is an end-to-end controlled fixture, not a natural telomere benchmark.
- The unified suite additionally covers TelSeeker Ctg positive/resume behavior, a single-ID negative, real hifiasm chromosome-only extension, ONT FASTA, both join directions, best-not-first selection, unsupported-leader fallback, concurrent index sharing, bounded batching without candidate loss, ambiguous cuts, cache corruption and operational tool failures.

## Unified checks and report evidence

- Exact delivered source, local: **256 tests, 0 skipped**, 327.210 s.
- Exact delivered source, ARM/SSH 180: **256 tests, 5 skipped**, 152.069 s, job `43012102`. The remote skips are Node-dependent checks; all **six Node fixture tests** pass locally. Do not interpret skipped checks as executed passes.
- The preceding scientific-matrix snapshot also ran 256 tests remotely (five Node skips), job `43012024`. Only two report-directory enumeration functions and their test changed afterwards: hidden `.telomere-target-cache/` is excluded from task rows and per-task reports. All scientific files remain byte-identical; final local and remote suites cover the report-only fix.
- HTML authoring used the requested **AGY `gemini-4-neon-high`**. Reports/logs separate pool, retrieved/refined/eligible records, batches, support lower bounds, witness IDs, rejection reasons, chromosome assembly growth, selected tel-read tail and tested junction span. The global length column correctly identifies full linked-sequence length. Missing legacy evidence stays unavailable.
- The main inspection archive verifies **883 payload files, 16 historical/matrix HTML reports and 61 local links**, with no external assets. Its original CLI report includes a cache-directory pseudo-task; this historical artifact is preserved. The final supplement rerenders the complete original CLI workspace into **three reports for exactly two real tasks**, without changing biological outputs. Manifest hashes, local links and actual inline JavaScript syntax pass lightweight checks. No browser validation is claimed.

Earlier diagnostic failures are preserved: an initial two-CPU memory-limit failure, an old single-read-positive fixture corrected to two IDs, and a deliberately stopped shared-filesystem random-write index build. The incomplete index was never accepted or reused. Completed original ONT decisions match the sequential-build rerun.

## Reproducible evidence

Runtime/test source fingerprints hash the sorted JSON mapping of all 57 file SHA-256 values:

- Scientific matrix: `56d584910f82c16b8e7a56daa8330bbb32973a2d909f4879e6b36d13ffa53b34`.
- Delivered report-final snapshot: `80520ad5279be050a4756c4a8d8f02281f6cde457799456434f20a758b53c6a3`.

On SSH `180`:

```text
/share/home/xbzhang26/cgat-reads-debug/degap-tel-targets-final-20261122/degap-telomere-targets-evidence-20261122.tar.gz
SHA-256: e98a953aea02c8c085c627275f8460544848c0e73a9201ec44aa844a8213046c
Size: 23,793,850 bytes

/share/home/xbzhang26/cgat-reads-debug/degap-tel-targets-report-final-20261122/degap-telomere-targets-report-supplement-20261122.tar.gz
SHA-256: 262dae70fc0233d598102d1ed1ee8efee963f9e0d50fd7b18dcdeda8991c5498
Size: 1,419,717 bytes
```

These archives contain inspection evidence, not relocated resumable workspaces; full raw datasets and shared indexes are not promised. The original inputs and completed workspaces remain unchanged. No GPM bundled runtime refresh, release, merge or tag is included.

Detailed results, source hashes, jobs, dataset provenance and limitations: [machine-readable summary](2026-11-22-telomere-target-search-summary.json).
