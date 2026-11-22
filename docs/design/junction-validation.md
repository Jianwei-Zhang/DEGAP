# Junction validation and terminal flank repair

GapFiller, AutoGapfiller and CtgLinker use the same connection checks. A successful overlap is evidence for a connection, not permission to replace an entire original flank. TelSeeker remains single-ended and does not use this paired-boundary repair policy.

## Connection checks

1. Compare each new cumulative version with every existing opposite-family version, including version 0.
2. Try the existing terminal window first (`20,000 + edge` bp). If it yields no validated candidate, compare complete sequences. A failed cropped comparison is not a final no-overlap decision.
3. Accept a qualifying single alignment or a monotone, nonoverlapping collinear chain. Chains may bridge only small, coordinate-consistent interruptions within `edge`; the difference between the two interruption lengths is at most `max(20 bp, 5% of the larger interruption)`. Unaligned gaps do not count toward the 10 kb minimum effective overlap.
4. Full checks can also recognize containment. Physical sequence coordinates, not cropped-window boundaries, determine terminal geometry. Unique original-flank anchors must still validate the composition.
5. Validate original boundaries and the output composition before ranking candidates. Keep the existing identity / effective-overlap / overhang / stable-version ordering and publish only one winner.

Raw alignment stages, full fallback and pair decisions are cached separately. Their contracts identify content and tools; old v1 no-overlap decisions are not reused under the v2 adapter. Tool failure, malformed alignment and corrupted completed receipts raise errors rather than becoming biological no-overlap stops.

## Original boundaries and gap sequence

`L0` and `R0` are the original flanks, while later versions include cumulative extension. Align both originals to the provisional join using minimap2 `asm20` with a long `cs` path. A unique primary forward alignment must have MAPQ >= 20 and identity >= 95%; it normally covers the gap-facing original boundary. At least 1 kb of alignment is required, or the complete original if shorter.

Project boundaries through substitutions and indels, rather than subtracting original lengths from the provisional output. The sequence between the projected original boundaries is the recovered gap. Unmapped or ambiguous boundaries remain ineligible.

A soft-clipped end can hide the error the inward retry is meant to bypass. With a unique interior anchor and clipping wholly inside the protected-end repair policy's 2 kb window, test a **length-preserving** projection of the clipped tail into the draft. Require ordered, nonoverlapping projected boundaries and independent read support.

For each qualifying spanning read, realign the 32 bp of adjacent context on each side independently, allowing up to four edits per side. Each cut must be uniquely best without scoring the proposed replacement. The sequence between those cuts must then match the **whole** proposed replacement exactly to support it. Reads with different alleles or unresolved cuts remain in the support denominator; unresolved counts are a subset of `other_reads`. Separate edits touching the context are not applied on top of the projection. If support falls short, reject the candidate and record the votes in `terminal-projection-rejected.json` in the candidate evidence directory. This check does not infer a net length change for a clipped segment.

## Bounded end edits

The default repair window is at most **2 kb from each gap-facing original end**, leaving at least 1 kb of original sequence outside the repair window. Differences elsewhere are recorded but the original sequence is retained. These are conservative implementation defaults, not empirically optimized accuracy thresholds or new required CLI arguments.

For proposed changes inside that window:

- Map the available normalized original HiFi/ONT reads to the protected-flank + recovered-gap draft, not just the reads selected for extension assembly.
- Generate candidates from both assembly differences and local read-consensus alleles. This check also runs when the original end error was inherited unchanged by every cumulative version; it does not require an assembly to propose the correction first. Candidate discovery alone never authorizes an edit.
- Read-derived substitutions can be evaluated at a located boundary. Boundary-touching insertions/deletions are not inferred from reads alone, because their assignment to the gap versus the original flank is ambiguous.
- Count each distinct platform/read ID at most once; discard IDs with duplicate qualifying primary records. Ignore secondary, MAPQ < 20 and < 90% identity alignments.
- Require a 500 bp interior original anchor beyond the repair window and at least 200 bp of alignment on each side of the proposed edit.
- Require at least **5 alternate-supporting reads**, alternate fraction >= **80%** among all qualifying spanning reads, and support on both strands. Other alleles count against the fraction. Left-normalize pure indels before comparing support.
- If several changes are accepted on one end, require at least five shared supporting reads across those changes. Do not compose overlapping conflicting edits.

The support-anchor check remains conservative for short inputs: flanks shorter than 2.5 kb cannot supply its current 2 kb + 500 bp interior support interval and are not automatically revised. They can still participate in filling if their connection evidence qualifies.

For an already located boundary, insufficient or conflicting support retains the original allele and records the proposed change as unapplied. A projected clipped boundary instead needs positive evidence before it can be used. A difference from an assembly comparator is not by itself proof that either allele is correct. Reads supporting extension from both directions are not treated as independent merely because they belong to two families.

## One scheduled anchor-backoff retry

If extension from an initial paired-search endpoint stops normally, one subsequent attempt may move the anchor inward by up to 2 kb. It is a **separate scheduled attempt in the shared round budget**, not hidden extra assembly. The original snapshot is never changed. A retry must extend beyond the original cumulative length to create a new version; otherwise the family stops while other families continue. A requested retry may never run if the other side closes the gap or the budget ends first. `anchor_retries` contains only executed retries; `anchor_backoff_requests` also records unexecuted requests. The initial failed attempt already consumes a round, but an unexecuted retry does not consume an extra one. No repeated backoff is applied to later cumulative versions. TelSeeker does not use this retry.

## Composition and consumers

- **Positive/zero gap:** preserve original cores, apply only supported bounded terminal edits, and insert the recovered gap.
- **True original overlap:** remove the duplicate right-original prefix and retain the left-original overlap copy. Record the trim explicitly; this is deduplication, not implicit polishing or preservation of both copies.
- **AutoGapfiller:** integrate only the actual gap and explicitly approved end edits, not the whole flank-containing local sequence. Conflicting bounded replacement intervals are rejected rather than silently clipped.
- **CtgLinker:** original components retain their coordinates/orientation. Gap/repair sequence is represented by separate provenance components. AGP must reconstruct the FASTA exactly.

Gap-only length, terminal-edit length delta and overlap trimming are separate fields. For example, a 5,000 bp recovered gap plus a 5 bp flank deletion is reported as gap length 5,000 bp and flank delta -5 bp, not a 4,995 bp gap.

## Evidence and reports

`family-result.json` adds `junction`, `connection_diagnostics`, `anchor_retries`, `anchor_backoff_requests` and, for validated joins, `integration_sequence`. Existing manifests remain readable; old workspaces are not silently resumed with changed algorithm contracts.

- `connection-diagnostics.json`: fast/full stages and explicit rejection reasons for all checked pairs.
- `junction-evidence.json`: original sources, composed pieces, policy and per-edit support decisions.
- `junction-integration.fa`: the minimal AutoGapfiller replacement payload; it can be empty for a pure overlap removal.
- `process.log`: winner, gap-only length, flank delta, original overlap trim, applied/unapplied edits and budgeted anchor retries.

Family, CtgLinker and AutoGapfiller HTML reports display these fields separately. Detailed alleles/support and rejected comparisons remain available in expandable tables. The new offline HTML helper was implemented by AGY with the user-selected `gemini-4-neon-high` model. New algorithm tests and real-read validation are documented separately; report styling is not biological acceptance.
