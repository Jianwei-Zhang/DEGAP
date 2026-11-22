# TelSeeker: shared target search without global grouping

TelSeeker and TelSeeker Ctg keep **one chromosome-side extension family per endpoint**. Telomeric reads remain static targets: no global subtelomere clustering, one-family-per-read expansion, or reverse assembly is introduced. Round and growth budgets still count chromosome-side extension attempts; the selected tel-read tail is a final connection component, not another assembly round.

## Search and selection

1. Build one immutable, disk-backed SQLite target index per pool content in the endpoint workspace's shared `.telomere-target-cache/`. Construction is locked across workers. Postings are written sequentially and the covering lookup index is built once at the end, avoiding per-seed random B-tree writes on shared filesystems. Store sequences once; identical records with the same read ID are deduplicated. Different orientations/sequence variants of one ID never count as independent support.
2. Sample complex 21-mers every 20 bases in each target. Query **all** complex 21-mers unique within the complete current chromosome sequence. Keeping the complete query avoids losing matches that cross extension-round boundaries. Short tandem repeats and low-complexity keys cannot assign a chromosome. Valid ambiguous IUPAC DNA bases are retained in sequences but cannot contribute exact anchor keys.
3. Retrieve all matching records, without a top-N cap. Process all retrieved records in multi-query MUMmer batches (at most 256 records or 8 MiB per batch), then inspect every alignment. These are batch-size bounds, not a cap on candidates. There is no per-read tool invocation or temporary FASTA file. No alignment tool is run if no index candidate exists.
4. Require forward, end-compatible alignment, at least 10 kb on both sequences, identity at least 80%, and at least three locally unique exact anchor keys spanning 500 bp. These are conservative implementation defaults, not empirically optimized accuracy thresholds. The anchors are unique within the current chromosome query, **not proven unique genome-wide**.
5. Sort eligible records by identity, effective overlap, endpoint overhang, then stable read ID / sequence / coordinate order. Evaluate read corroboration in that order. The first supported candidate wins; a higher-ranked unsupported candidate does not hide a lower-ranked supported one. All retrieved alignments have already been compared before this decision.

The index is an explicit sampled-anchor eligibility filter, not a guarantee of finding every arbitrary 80%-identity alignment. A 40-base exact match includes a stride-20 sampled 21-mer, but that key must also pass the complexity and local-uniqueness rules. Searches without qualifying keys remain unconnected; the program does not conclude that no biological telomere exists. “Best” means best under this retrieval/validation policy, not a mathematical optimum over every possible alignment in the raw pool.

## Read corroboration and composition

A connection needs at least **two distinct read IDs**, sharing a nonrepeat anchor spanning 500 bp and agreeing across up to **200 outward bases immediately beyond the junction**. Shorter additions are tested over their available length. Comparison permits up to 20% edit distance. The selected read is included in the count; evaluation can stop after finding enough evidence, so reported support counts are **lower bounds**, not total depth. The witness IDs and tested spans are recorded.

This supports the local connection; it is **not full-read polishing, consensus of the entire added tail, a telomeric repeat-length estimate, or proof of final telomere completeness**. Divergence beyond the tested junction-proximal span can remain. Do not interpret the number of supporting IDs as proof that every added base is correct.

At the final join, retain the complete current chromosome-side sequence. If MUMmer leaves a small chromosome overhang within `edge`, project it onto the read using bounded edit alignment (at most 1 kb); reject ambiguous or unsupported cuts. Append/prepend only the read sequence outside that projected boundary. The tel read never replaces the current chromosome overlap. This does not claim to polish earlier chromosome-side assemblies.

## Evidence, caching and compatibility

Each direct or extension check publishes:

- `connection-diagnostics.json`: pool, retrieved/refined/eligible counts, rejection reasons, selected read and support scope;
- `candidates.tsv`: all retrieved records and their alignment/selection status;
- `connection-support.json`: ordered support checks, witness IDs and tested spans;
- `telomere-check.json`: content-bound receipt, including negative decisions and copied artifacts.

Raw multi-query alignments are cached under each endpoint's `telomere-alignments/`, shared by round zero and later checks of identical sequence content. Indexes are shared across endpoints; their input/content contracts and files are verified before use. Tool failures, malformed alignments or modified completed artifacts are errors, not biological no-hit results.

Both normal TelSeeker and TelSeeker Ctg use the same matcher. Old first-hit checks are not accepted merely because a log or linker exists: use a new output directory when changing algorithm contracts. Original reads and old result workspaces are not rewritten. No new required CLI parameter is introduced.

Validation evidence and limitations: [2026-11-22 SSH 180 controls](../validation/2026-11-22-telomere-target-search.md).
