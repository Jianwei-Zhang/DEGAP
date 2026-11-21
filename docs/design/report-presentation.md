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
