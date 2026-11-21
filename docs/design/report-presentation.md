# Report presentation refresh

Date: 2026-11-21. Branch: `feat/bidirectional-gap-extension`.

## Scope

Generated reports share a restrained, offline scientific-report style: neutral
backgrounds, a teal accent, compact headers/metrics, consistent sections, thin
borders, readable tables and clear text-based statuses. Large decorative emoji
headings and gradients are removed. Long tables scroll inside their container on
narrow screens; keyboard focus and print styles are retained.

Applied to family gap reports, CtgLinker family reports, AutoGapfiller overview,
and TelSeeker overview/endpoint reports. `bin/ReportStyle.py` supplies inline CSS;
there are no new remote resources, fonts or JavaScript dependencies.

The design was implemented with the explicitly requested AGY model
`gemini-4-neon-high`. Its run was interrupted after producing the shared style,
family/CtgLinker reports and AutoGapfiller changes. Integration completed the
TelSeeker templates using that shared style and removed redundant summary values.
No alternative design model was invoked.

## Data and compatibility

- Search, ranking, budgets, sequence assembly, result schemas and parsers are not
  changed. An AST comparison checked all algorithm/parser functions in affected
  modules against the prior commit.
- All family versions, eligible comparisons, selected results, strands, counters,
  stop reasons, and artifact/provenance/checkpoint links remain available.
- Existing English/Chinese switches, chromosome filters, table expansion and
  chromosome tracks remain functional.
- Styling applies to newly generated HTML. Existing content-verified result
  directories are not silently rewritten to replace their published reports.

## Validation

- 48 focused report/result/consumer tests pass, including five new presentation
  tests for offline HTML, required evidence, escaping and language hooks.
- Local cached Chromium renders six report variants at 1440 px and 390 px:
  **12 viewport checks, no page-wide horizontal overflow and no JavaScript
  runtime exceptions**. Family tables, AutoGapfiller and TelSeeker screenshots
  were visually reviewed. English/Chinese switching and AutoGapfiller table
  expansion, chromosome filtering and track selection were exercised.
- Browser page requests to HTTP/HTTPS were blocked; no browser/package setup was
  needed. The preview data are controlled fixtures, not biological results.
- Local previews and browser evidence: `/tmp/degap-report-style-preview`.
- No genomic mapping/assembly runs were repeated for this presentation-only change.

## Follow-up alignment (2026-11-22)

- Updated `bin/AutoGapfillerVisualizer.py` task-table headers from `Extended` / `是否延伸` and `Flag` to `Filled` / `是否填补` and `Direction` / `方向/模式` so gap closure is not conflated with extension and deprecated `--flag` wording is removed.
- Updated `bin/TelSeekerVisualizer.py` presentation templates to escape dynamic identifiers/filenames/labels, render unlimited (`None`) round/length limits as `No limit` / `无限制`, and fix the stale `telo_part2_visual.py` footer name to `TelSeekerVisualizer.py`.
- `bin/FamilyResult.py`, `bin/FamilyCtgLinker.py`, and `bin/ReportStyle.py` were inspected and required no changes.
- AGY `gemini-4-neon-high` made these presentation edits. AST comparison confirms only rendering functions changed; no parsers, schemas, algorithms, or JS control flow were modified.
- Final unified regression: **181 tests passed, 0 skipped, 275.821 seconds**, including new coverage for headers, unlimited/zero limits, identifier escaping and the translated HTML layer. Log: `/tmp/degap-docs-final-regression-20261122.log`. No browser was launched for this follow-up.
- Separate real-read validation and report/result checks are recorded in [the new SSH 180 dataset report](../validation/2026-11-22-ssh180-new-dataset.md).

## Junction evidence follow-up (2026-11-22)

AGY `gemini-4-neon-high` implemented `bin/JunctionReport.py` and its focused presentation tests. Integration adds the helper to family, CtgLinker and AutoGapfiller reports. The new sections distinguish fast/full connection diagnostics, rejected pairs, gap-only length, flank delta, original-overlap trimming, supported/unapplied terminal edits, and budgeted anchor retries. They remain offline and expose full evidence through expandable tables. Existing language/chart controls are retained.

Unlike the earlier presentation-only changes above, this follow-up accompanies algorithm and composition changes documented in [junction validation](junction-validation.md). Automated report/schema checks are required; no new browser validation is claimed.


### Final AutoGapfiller coordinate alignment

AGY `gemini-4-neon-high` corrected the inclusive-gap arithmetic, downstream flank-delta shifts, zero/negative gap handling and numeric genomic sorting in `_create_full_html`. Tracks are explicitly labeled as original coordinates and use original gap spans. A source/AST audit confirms that only this rendering method and its tests changed in the final follow-up. Six focused tests pass with Node; actual rendered report projections also match remaining N intervals from the final FASTA. No browser validation is claimed. Results, limitations and the portable evidence archive are recorded in [SSH 180 junction validation](../validation/2026-11-22-junction-repair.md).
