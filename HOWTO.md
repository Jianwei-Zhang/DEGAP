![DEGAP](./DEGAP.png)

# DEGAP 2.0 — Users' Guide

Practical reference for installing, running, resuming, and inspecting **DEGAP 2.0** workflows.

## Table of Contents

- [Installation](#installation)
- [Input Reads & ONT Quality Contract](#input-reads--ont-quality-contract)
- [Extension Model, Budgets, Errors & Resume](#extension-model-budgets-errors--resume)
- [Execution Modes & Examples](#execution-modes--examples)
  - [1. GapFiller](#1-gapfiller)
  - [2. CtgLinker](#2-ctglinker)
  - [3. TelSeeker & TelSeekerCheck](#3-telseeker--telseekercheck)
  - [4. TelSeeker Ctg](#4-telseeker-ctg)
  - [5. AutoGapfiller](#5-autogapfiller)
- [Parameters Reference](#parameters-reference)
- [Output Layouts](#output-layouts)
- [Troubleshooting](#troubleshooting)
- [Citation & Contact](#citation--contact)

---

## Installation

DEGAP 2.0 targets Linux. The supplied environment pins Python 3.10; Python packages and external tools are declared in [environment.yml](environment.yml):

```bash
micromamba env create -f environment.yml
micromamba activate degap_env
```

Update after `environment.yml` changes:

```bash
micromamba env update -f environment.yml --prune
```

### Manual Installation (without Micromamba)

Python packages:

```bash
pip install biopython pysam numpy matplotlib scipy
```

External tools required in `PATH` (versions pinned in `environment.yml` where noted):

- **Minimap2** (2.28): <https://github.com/lh3/minimap2>
- **Hifiasm** (0.25.0): <https://github.com/chhylp123/hifiasm>
- **SAMtools** (1.22.1): <http://www.htslib.org/>
- **Seqkit** (2.8.0): <https://github.com/shenwei356/seqkit>
- **MUMmer4**: <https://github.com/mummer4/mummer>
- **GNU Parallel**: <https://www.gnu.org/software/parallel/>
- **Ripgrep** (15.1.0): <https://github.com/BurntSushi/ripgrep>

---

## Input Reads & ONT Quality Contract

Provide at least one of `--hifi` or `--ont`. Both flags accept one or more files (for example `--hifi HIFI1.fq.gz HIFI2.fq.gz`).

| Input Mode | Supported File Formats | Assembly Requirements |
| --- | --- | --- |
| **HiFi-only** (`--hifi`) | `.fa`, `.fasta`, `.fna`, `.fq`, `.fastq` (and `.gz`) | HiFi FASTA and FASTQ are both supported. |
| **ONT-only** (`--ont`) | **Original `.fq` / `.fastq` (and `.gz`) required** | `hifiasm --ont` requires genuine Phred+33 qualities. FASTA-only ONT input is rejected at assembly time. |
| **Mixed** (`--hifi` + `--ont`) | HiFi FASTA/FASTQ; ONT FASTQ recommended | ONT FASTA is supported only as supplementary `--ul` input while HiFi reads remain available in that round. If a round falls back to ONT-only assembly, original ONT FASTQ qualities are required. |

How ONT quality recovery works (`bin/ONTAssemblyReads.py`):

- Although preprocessing converts reads to working FASTA and the CLI parser accepts FASTA paths, ONT-only local assembly never fabricates synthetic quality scores.
- Original ONT FASTQ/FASTQ.gz files (single or multiple) are indexed lazily into a compressed SQLite lookup (`ont-quality-cache/`) shared across families and rounds within a search.
- Selected reads are restored to FASTQ in their selected orientation (reversing qualities for reverse-complement reads) with strict read-ID, sequence, and unique-trim validation. Missing IDs, duplicate IDs, mismatched sequences, ambiguous trims, or malformed FASTQ inputs fail explicitly.
- Original ONT files — including their quality values — are part of the checkpoint contract; quality-only input edits invalidate checkpoint reuse.

---

## Extension Model, Budgets, Errors & Resume

### Bidirectional Historical Families (GapFiller, AutoGapfiller, CtgLinker)

- **GapFiller / AutoGapfiller**: Left (`L0, L1, ...`) and right (`R0, R1, ...`) chains extend in alternating turns. Each version is a cumulative sequence. After round-0 direct closure is checked (`L0` vs `R0`), each new version is aligned against **all** historical versions of the opposite side (including `L0`/`R0`). If one side stops normally (e.g., no extension reads/contigs or no growth), its history stays frozen as an overlap target while the other continues.
- **CtgLinker**: Each outer-end search builds a primary family for the active contig/scaffold plus candidate families for every remaining eligible contig/scaffold in both forward (`+`) and reverse-complement (`-`) orientations. Primary and candidates alternate within that endpoint search, sharing one round and length budget. An original input contig cannot be incorporated twice.
- **TelSeeker / TelSeeker Ctg**: Uses the same persistent single-round engine, but remains **single-ended** with one real endpoint family and no fictitious opposite family. Connection requires matching a telomeric read, not overlapping the initial seed.
- **Deterministic winner selection**: After checking the current comparison frontier, eligible connections are ranked by:
  1. Alignment identity (descending)
  2. Effective overlap length (descending)
  3. Endpoint overhang (ascending)
  4. Fixed family/version/coordinate order

  Only one best eligible pair is selected and emitted. Note that highest alignment identity in a repetitive region does not by itself prove biological correctness; inspect the exported evidence and alignments.

### Shared Budgets

- `--MaximumExtensionRound` counts attempted extension rounds **summed across all families in one search** (both sides in GapFiller/AutoGapfiller; primary + all candidates in one CtgLinker endpoint search). Failed normal attempts still consume rounds.
- `--MaximumExtensionLength` sums newly added bases across all families in that search. Original flanks and historical snapshots are excluded, and final overlap is not subtracted. Winning with an older historical version does not refund already spent rounds or added bases.
- Length is checked **before** launching the next attempt, so the final complete attempt may exceed `--MaximumExtensionLength` and is still kept whole and checked for closure before stopping.
- Setting `--MaximumExtensionRound 0` skips extension attempts while still allowing round-0 direct closure. Negative limits are rejected. Both limits default to `None` (unlimited).

### Outcome States & Tool Errors

Treat **tool ran (exit 0)**, **gap closed / connected**, and **sequence biologically verified** as distinct states:

- Explicit `hifiasm` errors on stdout/stderr fail even if `hifiasm` exits with status `0`.
- An absent primary-contig GFA (`.p_ctg.gfa`) is treated as a tool error, not biological exhaustion.
- A genuinely completed `hifiasm` run that writes a valid empty GFA (0 sequence segments) and yields no usable extension contig is distinct from a tool error.
- Tool errors never write a completed assembler/assembly stage receipt and are never published as normal `exhausted` results. Completed earlier stages (such as `mapping/` and `reads/`) remain reusable after fixing the error.

### Resume, Legacy Directories, `--flag`, and `--remove`

- **Resuming**: Rerun the exact same command (or add `--resume_auto`) to resume from content-verified family checkpoints (`families.json`, `ctg-families.json`, and stage receipts). Uncommitted interrupted stages are archived into `.incomplete.*` directories rather than trusted.
- **Unsupported legacy resume**: `--resume ROUND` is rejected for all four `DEGAP.py` modes (no arbitrary round rewind). Legacy single-sided output directories (from pre-family runs) are not silently resumed or migrated; keep old directories intact and use a **new output directory**.
- **`--flag`**: Deprecated compatibility option in `gapfiller` mode. Passing it prints a notice and ignores it; omit it from new commands.
- **`--remove`**: Legacy cleanup flag (`1`, `2`, `3`; default `2`). In family workflows, prepared reads, split chunks, indices, historical sequences, and verified stages are **always retained** at `--remove 1` and `--remove 2` because recovery depends on them. Remove a finished workspace manually only when histories and resume are no longer needed.

---

## Execution Modes & Examples

### 1. GapFiller

Fill one gap between known left and right flanking sequences.

```bash
python bin/DEGAP.py --mode gapfiller \
    --seqleft ./gap_left.fa \
    --seqright ./gap_right.fa \
    --hifi ./hifi_reads.fq.gz \
    --ont ./ont_reads.fq.gz \
    --out ./output_gapfiller \
    --MaximumExtensionRound 25 \
    -t 20
```

HiFi-only with optional depth filtering and k-mer filtering:

```bash
python bin/DEGAP.py --mode gapfiller \
    --seqleft ./gap_left.fa \
    --seqright ./gap_right.fa \
    --hifi ./hifi_reads.fa \
    -o ./output_gapfiller_hifi \
    --filterDepthHifi 0.3 \
    --kmer_filter \
    -t 20
```

Check `family-result.json`: only `status: "closed"` means the gap was closed into one sequence. Otherwise `<name>.final.fa` contains separate unclosed family records.

### 2. CtgLinker

Link contigs and close gaps from contig edges.

```bash
python bin/DEGAP.py --mode ctglinker \
    --ctgseq ./contigs.fa \
    --hifi ./HIFI1.fq.gz ./HIFI2.fq.gz \
    --ont ./ONT1.fq.gz ./ONT2.fq.gz \
    --out ./output_ctglinker \
    --MaximumExtensionRound 25 \
    -t 20
```

ONT-only example (must use original FASTQ/FASTQ.gz):

```bash
python bin/DEGAP.py --mode ctglinker \
    --ctgseq ./contigs.fa \
    --ont ./ont_reads.fq.gz \
    --out ./output_ctglinker_ont \
    -t 20
```

`DG.scaffold.agp` reconstructs `DG.scaffold.fa` from `DG.scaffold.components.fa` (which includes actual extended snapshots, not just unmodified input contigs).

### 3. TelSeeker & TelSeekerCheck

#### Step A (Optional): Prepare manual chromosome-end review files

```bash
python bin/TelSeekerCheck.py \
    --genome ./genome.fa \
    --motif TTAGGG \
    --out ./telomere_check
```

`TelSeekerCheck.py` writes review files into `<out>/genome.telomere.check/`:

- `genome.telomere.check.left.2kb.fa` and `genome.telomere.check.right.2kb.fa` (default `--extract-length 2000`)
- `all_chromosomes_combined.png` and `<chromosome>_telomere_motif.png` (1 kb plotting window by default, `--plot-window-size 1000`)

It does **not** emit automatic classification lists (`genome.telomere.check.csv`, `need_extension_chr_end.txt`, or `uncertain_chr_end.txt`). Inspect the plots and end sequences, then list the ends you want to extend (one `<chromosome_id>.L` or `<chromosome_id>.R` per line; blank lines and `#` comments are allowed):

```text
# Target ends selected after reviewing TelSeekerCheck plots
Chr01.L
Chr03.R
```

#### Step B: Run TelSeeker

Pass targets directly (`-e Chr01.L Chr03.R`) or via a file (`-e target_ends.txt`):

```bash
python bin/DEGAP.py --mode telseeker \
    --genome ./genome.fa \
    --motif TTAGGG \
    -e ./target_ends.txt \
    --hifi ./hifi_reads.fq.gz \
    --ont ./ont_reads.fq.gz \
    --out ./output_telseeker \
    --work 4 \
    --tel-n 100 \
    --tel-r 0.6 \
    --tel-mm 0 \
    --MaximumExtensionRound 25 \
    -t 20
```

Common uppercase A/C/G/T motifs: `TTAGGG` (vertebrates), `TTTAGGG` (plants), `TTAGG` (insects).

**Terminal telomeric read discovery**: Step 1 scans `--tel-n * len(--motif)` bases at both physical ends of each read using motif rotations and reverse-complement rotations. A read end passes when `hits / --tel-n >= --tel-r`, allowing `--tel-mm` (`0` or `1`) mismatches per motif-length unit. With `TTAGGG`, `--tel-n 100` checks 600 bp and `--tel-r 0.6` requires at least 60 hits. If Step 1 finds too few telomeric reads, relax `--tel-r`, set `--tel-mm 1`, or lower `--tel-n`.

A connected linker in TelSeeker means a linker sequence was generated and integrated into `part3.integration.results/final.genome.fa`; final telomere completeness must still be confirmed by manually reviewing `part3.integration.results/final.genome.telomere.check/`.

### 4. TelSeeker Ctg

Extend one or more standalone single-record contig FASTA endpoints (`L` or `R`) while sharing preprocessed reads and `part1.telo.reads/`.

```bash
python bin/DEGAP.py --mode telseeker_ctg \
    --ctg ./chr1_partial.fa L \
    --ctg ./chr2_partial.fa R \
    --motif TTAGGG \
    --hifi ./hifi_reads.fq.gz \
    --ont ./ont_reads.fq.gz \
    --out ./output_telseeker_ctg \
    --work 4 \
    --tel-n 100 \
    --tel-r 0.6 \
    --tel-mm 0 \
    --MaximumExtensionRound 25 \
    --MaximumExtensionLength 1000000 \
    -t 20
```

Rules for `telseeker_ctg`:

- Each `--ctg <fasta> <L|R>` file must contain **exactly one** FASTA record.
- Reusing `--out` for additional `--ctg` endpoints is permitted only when `--hifi`, `--ont`, `--motif`, `--tel-n`, `--tel-r`, and `--tel-mm` match `telseeker_ctg.input/cache_manifest.json`; otherwise DEGAP stops with a mismatch error.
- `--filterDepthHifi` and `--filterDepthOnt` are not supported in `telseeker_ctg`.
- Aggregate results are written to `result/status.tsv`. Endpoints without matching telomeric reads fail explicitly with `no_tel_reads` or `no_tel_reads_for_target_end`.

### 5. AutoGapfiller

Detect all N gaps in a draft genome (`>= 20 bp`, merging adjacent gaps `< 11,000 bp` apart), create **one bidirectional job per gap**, execute jobs in parallel, integrate closed gaps, and build an HTML report.

```bash
python bin/AutoGapfiller.py \
    --genome ./genome_with_gaps.fa \
    --hifi ./hifi_reads.fq.gz \
    --ont ./ont_reads.fq.gz \
    --output ./output_autogap \
    --work 4 \
    --MaximumExtensionRound 25 \
    --MaximumExtensionLength 500000 \
    --kmer_filter \
    -t 20
```

Use `--setup_only` to run Steps 1–3.1 (preprocessing, gap detection, and job directory / `run.sh` creation) without running gap filling:

```bash
python bin/AutoGapfiller.py \
    --genome ./genome_with_gaps.fa \
    --hifi ./hifi_reads.fq.gz \
    --output ./output_autogap \
    --setup_only
```

---

## Parameters Reference

### `bin/DEGAP.py`

| Parameter | Description | Default | Applies To |
| --- | --- | --- | --- |
| `--mode` | `gapfiller`, `ctglinker`, `telseeker`, or `telseeker_ctg` | Required | All |
| `-o, --out` | Output directory | Required | All |
| `--hifi` | One or more HiFi read files (FASTA/FASTQ, optional `.gz`) | `None` | All (at least one of `--hifi`/`--ont`) |
| `--ont` | One or more ONT read files (original FASTQ/FASTQ.gz required for ONT-only assembly and mixed→ONT-only fallback) | `None` | All (at least one of `--hifi`/`--ont`) |
| `-t, --thread` | Thread count | `20` | All |
| `--edge` | Edge tolerance in bp for boundary/misassembly checks | `500` | All |
| `--MaximumExtensionRound` | Shared attempt cap per search (`0` allows round-0 closure only) | `None` | All |
| `--MaximumExtensionLength` | Shared added-base cap per search in bp (checked before next attempt) | `None` | All |
| `--kmer_filter` | Enable k-mer read pre-filtering before alignment | `False` | All |
| `-ks, --kmer_size` | K-mer size when `--kmer_filter` is enabled | `41` | All |
| `-kn, --kmer_num` | Number of k-mers when `--kmer_filter` is enabled | `20` | All |
| `--filterDepthHifi` | HiFi mapped-depth filter ratio (e.g., `0.3`); implemented in `gapfiller`, skipped in `ctglinker`/`telseeker`, rejected in `telseeker_ctg` | `None` | `gapfiller` |
| `--filterDepthOnt` | ONT mapped-depth filter ratio; same mode support as `--filterDepthHifi` | `None` | `gapfiller` |
| `--remove` | Legacy cleanup level (`1`, `2`, `3`); family histories, prepared reads, indices, and verified stages are retained | `2` | All |
| `--resume_auto` | Resume verified family state (same behavior as rerunning the command) | `False` | All |
| `--resume` | Unsupported legacy single-sided round selection (rejected with an error) | `None` | Unsupported |
| `--seqleft`, `--seqright` | Left and right flanking FASTA files | Required | `gapfiller` |
| `--flag` | Deprecated; prints a notice and is ignored | `None` | `gapfiller` |
| `--ctgseq` | Input contigs FASTA | Required | `ctglinker` |
| `--genome` | Input genome/chromosomes FASTA (`--chr` is deprecated) | Required | `telseeker` |
| `-e, --target_ends` | Target ends (`Chr01.L Chr01.R`) or one text file with one `<chr>.<L\|R>` per line | Required | `telseeker` |
| `--ctg` | Repeatable `<contig.fa> <L\|R>` single-record contig endpoint task | Required | `telseeker_ctg` |
| `-m, --motif` | Uppercase A/C/G/T telomere repeat motif | Required | `telseeker`, `telseeker_ctg` |
| `-tn, --tel-n` | Motif-length units scanned at each read end | `100` | `telseeker`, `telseeker_ctg` |
| `-tr, --tel-r` | Minimum terminal hit ratio (`hits / tel-n`), in `(0, 1]` | `0.6` | `telseeker`, `telseeker_ctg` |
| `-tmm, --tel-mm` | Allowed mismatches per motif-length unit (`0` or `1`) | `0` | `telseeker`, `telseeker_ctg` |
| `-w, --work` | Parallel chromosome/contig endpoint workers | `1` | `telseeker`, `telseeker_ctg` |

### `bin/AutoGapfiller.py`

| Parameter | Description | Default |
| --- | --- | --- |
| `-g, --genome` | Input draft genome FASTA with N gaps | Required |
| `-o, --output` | Output directory (note: `--output`, not `--out`) | Required |
| `--hifi` | One or more HiFi read files (FASTA/FASTQ, optional `.gz`) | `None` |
| `--ont` | One or more ONT read files (original FASTQ/FASTQ.gz required for ONT-only assembly and mixed→ONT-only fallback) | `None` |
| `-w, --work` | Parallel gap-filling workers | `1` |
| `--setup_only` | Prepare preprocessing, gap detection, and job scripts without executing gap filling | `False` |
| `-t, --thread` | Threads passed to each DEGAP job | `20` |
| `--kmer_filter` | Enable k-mer pre-filtering | `False` |
| `-ks, --kmer_size` | K-mer size | `41` |
| `-kn, --kmer_num` | Number of k-mers | `20` |
| `--MaximumExtensionLength` | Shared added-base cap per gap in bp | `None` |
| `--MaximumExtensionRound` | Shared attempt cap per gap | `None` |
| `--filterDepthHifi`, `--filterDepthOnt` | Optional depth filtering thresholds in Step 1 preprocessing | `None` |
| `--edge` | Edge tolerance in bp | `500` |
| `--remove` | Legacy cleanup level (`1`, `2`, `3`); family recovery files are retained | `2` |

### `bin/TelSeekerCheck.py`

| Parameter | Description | Default |
| --- | --- | --- |
| `--genome` | Input genome FASTA | Required |
| `--out` | Output directory (creates `<out>/genome.telomere.check/`) | Required |
| `-m, --motif` | Telomere repeat motif | `TTAGGG` |
| `--extract-length` | Terminal sequence length extracted to left/right FASTAs (bp) | `2000` |
| `--plot-window-size` | Whole-genome motif plot window size (bp) | `1000` |
| `--check-length` | End window length (bp) | `100000` |

> **Depth filtering note**: When `--filterDepthHifi` or `--filterDepthOnt` is active in `DEGAP.py --mode gapfiller`, `selectRawReads.py` retains reads mapping outside `(threshold * mode_depth, (2 - threshold) * mode_depth)`, reads near sequence ends, unmapped reads, and reads with `< 99%` aligned query length, writing FASTA files (`Genome.inputCtg.usefullReads.fasta` and, in mixed mode, `Genome.inputCtg.usefullReads.ont.fasta`).

---

## Output Layouts

### GapFiller

```text
output_directory/
├── <name>.final.fa                 # One closed sequence, or separate unclosed family records
├── family-result.json              # Verified status ("closed" vs unclosed), winner, checksums
├── family-report.html              # Offline bilingual report of both families & eligible pairs
├── family-components.fa            # Cumulative sequences referenced by process.agp
├── families.json                   # Checkpoint, family histories, and shared budget state
├── process.agp                     # AGP reconstructing <name>.final.fa from family-components.fa
├── process.log                     # Human-readable search outcome
├── process.summary                 # Per-family attempts, versions, added bases, stop reason
├── family-overlaps/                # Version-pair MUMmer alignments and receipts
├── ont-quality-cache/              # Shared SQLite quality index (created lazily for ONT-only assembly)
├── processed_reads/                # Normalized working FASTA reads
├── hifi_reads.idx / ont_reads.idx  # BioPython read indices (as applicable)
├── HiFi.reads.stat / ONT.reads.stat
├── hifi_reads_part/ / ont_reads_part/
└── attempts/<family_index>/<attempt>/
    ├── result.json                 # Verified attempt receipt
    ├── mapping/                    # Reusable minimap2/samtools BAM + receipt
    ├── reads/                      # Selected reads + stage receipt
    ├── assembler/                  # Reusable local assembly (hifiasm) + stage receipt
    └── assembly/                   # Contig selection + cumulative sequence + stage receipt
```

### CtgLinker

```text
output_directory/
├── DG.scaffold.fa                  # Final scaffolds (unjoined contigs retained)
├── DG.scaffold.log                 # Scaffold N50 / summary stats
├── DG.scaffold.agp                 # Reconstruction coordinates referencing DG.scaffold.components.fa
├── DG.scaffold.components.fa       # Actual input and extended snapshot sequences used by AGP
├── DG.project.fa                   # Joined assembly snapshots
├── Genome.inputCtg.fa              # N-split / seed-filtered input contigs
├── Genome.inputCtg.log             # Input contig stats
├── used.txt / unplaced.txt         # Original contig membership tracking
├── ctg-families.json               # Project checkpoint and pending endpoint state
├── ctglinker-result.json           # Verified result manifest and artifact checksums
├── ctglinker-family-report.html    # Offline bilingual report linking per-endpoint reports
├── ctg-sources/                    # Immutable input components and .rc.fa views
├── ctg-assemblies/                 # Composed DG*.fa scaffold snapshots
├── ctg-overlaps/                   # MUMmer overlap cache
├── shared-extension-stages/        # Content-addressed reusable single-round stages
└── project/
    └── endpoint-000001/            # One shared-budget endpoint search
        ├── families.json
        ├── attempts/
        ├── family-result.json
        └── family-report.html
```

### TelSeeker

```text
output_directory/
├── genome.telomere.check/
│   └── need_extension_chr_end.txt          # Target ends loaded from -e/--target_ends
├── part1.telo.reads/
│   ├── Global.left.telo.reads.fa
│   ├── Global.right.telo.reads.fa
│   ├── left.telo.reads.fa
│   ├── right.telo.reads.fa
│   └── part1.log
├── part2.chr.end.job/
│   └── <Chr>.<L|R>/                        # e.g., Chr01.L, Chr01.R
│       ├── seed.fa
│       ├── telo.reads.fa
│       ├── direct.check/
│       ├── linker.fa / linker.info         # Present when connected
│       └── extension/
│           ├── families.json               # Single real endpoint family checkpoint
│           ├── attempts/0/<attempt>/       # Attempts with reads/assembly/telomere receipts
│           ├── telomere-result.json        # Verified result manifest used by reports
│           ├── telomere-alignments/        # Cached MUMmer alignments
│           ├── extension.log
│           └── extension.summary
├── part3.integration.results/
│   ├── check_part2_jobs.csv                # Per-end connection summary
│   ├── final.genome.fa                     # Genome with connected linkers integrated
│   └── final.genome.telomere.check/        # Final 2 kb end FASTAs + motif plots for manual review
└── visual.report/
    ├── Global.report.html                  # Global overview report
    └── <Chr>.<L|R>.report.html             # Individual chromosome-end reports
```

### TelSeeker Ctg

```text
output_directory/
├── telseeker_ctg.input/
│   ├── cache_manifest.json                 # Shared read/motif/tel-param compatibility contract
│   ├── ctgs.fa
│   └── ctg_tasks.tsv
├── part1.telo.reads/                       # Shared telomeric read pool
├── telseeker_ctg.jobs/
│   └── <safe_ctg_id>.<L|R>/
│       ├── task-contract.json
│       ├── direct.check/
│       ├── extension/                      # Single-ended family checkpoint & telomere-result.json
│       └── result/
│           ├── status.tsv
│           ├── receipt.json
│           ├── extended.ctg.fa             # On connection success
│           └── failed.reason.txt           # On failure
└── result/
    └── status.tsv                          # Aggregate status across requested endpoints
```

### AutoGapfiller

```text
output_directory/
├── autogapfiller_summary.json              # Overall pipeline summary (or autogapfiller_setup_summary.json with --setup_only)
├── 01.reads_preprocessor/
│   ├── processed_reads/
│   │   ├── hifi_reads.fa                   # Present when --hifi supplied
│   │   └── ont_reads.fa                    # Present when --ont supplied
│   ├── HiFi.reads.stat / ONT.reads.stat
│   ├── hifi_reads.idx / ont_reads.idx
│   ├── hifi_reads_part/ / ont_reads_part/
│   └── preprocessing_config.json
├── 02.gap_detection/
│   ├── gapbase.log                         # All raw N regions
│   ├── gap.log                             # Filtered (>=20 bp) and merged (<11 kb apart) gaps
│   ├── gapDataBase/                        # Raw flanking sequences
│   ├── gapData/                            # Final flanking sequences used by jobs
│   └── gap_detection_config.json
├── 03.gap_filling/
│   ├── job_manager_config.json
│   ├── gap_filling_results.json
│   └── <chrom>.<gap_id>.bidirectional/     # One bidirectional GapFiller workspace per gap
│       ├── run.sh
│       ├── seqleft.fasta / seqright.fasta
│       ├── families.json
│       ├── family-result.json
│       ├── family-report.html
│       ├── <chrom>.<gap_id>.bidirectional.final.fa
│       └── attempts/
├── 04.genome_integration/
│   ├── genome.filled.fasta                 # Integrated genome
│   ├── integration.log
│   └── integration_config.json
└── 05.visualization/
    ├── Global.report.html                  # Interactive overview report
    └── <chrom>.<gap_id>.bidirectional.html # Per-gap family report (with rebased artifact links)
```

---

## Troubleshooting

1. **ONT-only assembly fails with a FASTA / quality error**
   - Supply original ONT `.fq` / `.fastq` (or `.gz`) files with genuine Phred+33 qualities.
   - Verify there are no duplicate read IDs across the original ONT FASTQ files and no truncated/malformed records.
   - If you previously ran an older version before the ONT quality repair (`c0fd8e3`), keep the old directory untouched and rerun in a **new output directory**.

2. **Explicit `hifiasm` error or missing `.p_ctg.gfa` vs. unclosed `exhausted` run**
   - Inspect `attempts/<family>/<attempt>/assembler/hifiasm/hifiasm.stderr.log` and `hifiasm.stdout.log`.
   - Explicit `hifiasm` error lines (even with exit code 0) and missing primary-contig GFA files stop the run as tool errors rather than publishing `exhausted`.
   - If `hifiasm` completed cleanly with a valid empty GFA (0 sequence segments) and no usable extension contigs were produced, the family stops normally (`no_extension_contigs`) and the search may finish `exhausted` without closing the gap.

3. **Legacy output directory or changed-input error on rerun**
   - Legacy single-sided outputs (`checkpoint.info`, `process/`, legacy `round*` directories, or legacy CtgLinker/TelSeeker outputs without family manifests) are rejected rather than overwritten. Use a new output directory.
   - Editing input reads (including ONT quality scores), flanks, parameters, or tool binaries invalidates content-verified checkpoints; use a new output directory for changed runs.

4. **K-mer filtering issues**
   - `--kmer_filter` requires `seqkit`, GNU `parallel`, and `rg` (ripgrep) in `PATH`.
   - If a gap fails to extend with `--kmer_filter` enabled, retry without `--kmer_filter` in a new output directory.

5. **AutoGapfiller detects no gaps**
   - Verify the input genome FASTA contains `N` runs of at least **20 bp** (`AutoGapfillerPart2.py` filters out N regions shorter than 20 bp and merges adjacent gaps closer than 11,000 bp).
   - Check `02.gap_detection/gapbase.log` and `02.gap_detection/gap.log`.

6. **AutoGapfiller job failures**
   - Inspect `03.gap_filling/<chrom>.<gap_id>.bidirectional/run.stdout.log`, `run.stderr.log`, and `process.log`, plus `03.gap_filling/gap_filling_results.json`.
   - Reduce `--work` or `-t` if memory or CPU limits are exceeded.

7. **TelSeeker finds no telomeric reads or fails to connect an end**
   - Verify `-m/--motif` uses uppercase `A/C/G/T` matching the organism (e.g., `TTAGGG`, `TTTAGGG`, `TTAGG`).
   - Verify `-e/--target_ends` chromosome IDs match the headers in `--genome` and use `<chr>.L` or `<chr>.R`.
   - If Step 1 produces empty `part1.telo.reads/Global.left.telo.reads.fa` / `Global.right.telo.reads.fa`, rerun in a new output directory with relaxed discovery parameters such as `--tel-r 0.5 --tel-mm 1` (or smaller `--tel-n`).

8. **TelSeeker Ctg fails with `no_tel_reads`, `no_tel_reads_for_target_end`, or cache mismatch**
   - Check `result/status.tsv` and `telseeker_ctg.jobs/<safe_ctg_id>.<L|R>/result/status.tsv`.
   - `no_tel_reads`: no global telomeric reads passed Part 1.
   - `no_tel_reads_for_target_end`: telomeric reads exist, but none for the requested `L` or `R` end.
   - Cache mismatch when reusing `--out`: `--hifi`, `--ont`, `--motif`, `--tel-n`, `--tel-r`, and `--tel-mm` must match `telseeker_ctg.input/cache_manifest.json`, or use a new `--out`.

9. **Disk usage**
   - Family recovery retains prepared reads, split parts, indices, historical versions, and verified stage receipts regardless of `--remove 1|2`. Plan storage accordingly and archive or delete finished workspaces manually once recovery is no longer needed.

---

## Citation & Contact

Huang, Y. et al. DEGAP: Dynamic Elongation of a Genome Assembly Path. *Briefings in Bioinformatics*, Volume 25, Issue 3, May 2024, bbae194. [DOI: 10.1093/bib/bbae194](https://doi.org/10.1093/bib/bbae194)

For questions, bug reports, or suggestions, contact <jzhang@mail.hzau.edu.cn>.
