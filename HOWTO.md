![DEGAP](./DEGAP.png)

# DEGAP v2 - Users' Guide

## Dynamic Elongation of a Genome Assembly Path

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Modes](#modes)
- [Parameters](#parameters)
- [Usage Examples](#usage-examples)
- [Output Files](#output-files)
- [Performance Optimization](#performance-optimization)
- [Troubleshooting](#troubleshooting)

## Overview

**DEGAP v2** is an enhanced version of DEGAP that supports multiple sequencing technologies (HiFi, ONT, or both) for gap filling and genome assembly improvement. The software uses dynamic elongation strategies to resolve complex gap regions.

### Key Improvements in v2
- Support for ONT reads in addition to HiFi reads
- Mixed mode combining HiFi and ONT data
- Separate depth filtering for different data types
- K-mer based read filtering for improved performance
- Enhanced parallel processing capabilities
- TelSeeker mode for telomere extension
- AutoGapfiller mode for automated whole-genome gap filling

## Installation

### Prerequisites
- Linux operating system
- Python 3.10 
- At least 32GB RAM (64GB+ recommended for large genomes)

### Recommended environment setup (Micromamba)

We recommend using [Micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) to create an isolated environment with all required Python packages and external tools. Micromamba is a lightweight, standalone package manager that is faster than conda.

1. Install Micromamba following the [official installation guide](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html):

   ```bash
   "${SHELL}" <(curl -L micro.mamba.pm/install.sh)
   ```

   **Tip for Micromamba users:** To use `mamba` command instead of `micromamba`, add this line to the end of your `~/.bashrc`:
   ```bash
   alias mamba='micromamba'
   ```
   Then run `source ~/.bashrc` to apply.

2. Clone this repository and change into its directory.
3. Create the environment using Micromamba with the provided `environment.yml` file:

   ```bash
   micromamba env create -f environment.yml
   ```

4. Activate the environment:

   ```bash
   micromamba activate degap_env
   ```

5. If `environment.yml` is updated in the future, you can update the environment with:

   ```bash
   micromamba env update -f environment.yml --prune
   ```

The `environment.yml` file pins the main Python version, core Python packages (Biopython, Pysam, NumPy, Matplotlib, SciPy) and the required external tools (Minimap2, Hifiasm, SAMtools, Seqkit, MUMmer4, GNU Parallel, Ripgrep), ensuring a reproducible runtime environment for **DEGAP v2**.

### Alternative: install Python packages with pip

If you do not use Micromamba, install the required Python packages manually:

```bash
pip install biopython pysam numpy matplotlib scipy
```

### Required External Tools

If you install tools manually (without Micromamba), ensure the following programs are installed and available in your `PATH`:

- **Minimap2** (2.28): https://github.com/lh3/minimap2
- **Hifiasm** (0.25.0): https://github.com/chhylp123/hifiasm
- **SAMtools** (1.22.1): http://www.htslib.org/
- **Seqkit** (2.8.0): https://github.com/shenwei356/seqkit
- **MUMmer4**: https://github.com/mummer4/mummer
- **GNU Parallel**: https://www.gnu.org/software/parallel/
- **Ripgrep** (15.1.0): https://github.com/BurntSushi/ripgrep

## Modes

### 1. GapFiller Mode
**Purpose**: Fill a specific gap when left and right flanking sequences are known.

**Use Cases**:
- Filling known gaps in scaffolds
- Resolving specific problematic regions
- Targeted gap closure

**Key Parameters**:
- `--seqleft`: Left flanking sequence (FASTA)
- `--seqright`: Right flanking sequence (FASTA)
- `--flag`: Direction to elongate (`left` or `right`)

### 2. CtgLinker Mode
**Purpose**: Link contigs and fill gaps in fragmented assemblies.

**Use Cases**:
- Scaffolding unordered contigs
- Filling multiple gaps in draft assemblies
- Improving assembly contiguity

**Key Parameters**:
- `--ctgseq`: Input contigs (FASTA)
- `--edge`: Edge length for error detection (default: 500)

### 3. TelSeeker Mode
**Purpose**: Extend telomeric regions and achieve chromosome-level assembly.

**Use Cases**:
- Telomere extension
- Chromosome-level assembly
- Complete genome assembly

**Key Features**:
- Explicit target chromosome ends via `-e/--target_ends`
- Independent telomere quality check with `TelSeekerCheck.py`
- Parallel processing of chromosome ends
- Four-step automated pipeline:
  1. Target chromosome end loading
  2. Telomeric reads extraction (motif-based filtering)
  3. Chromosome end extension (parallel DEGAP execution)
  4. Result integration and genome assembly
  4. Visualization (HTML reports and statistics)

**Key Parameters**:
- `--genome`: Input genome/chromosomes (FASTA)
- `--motif`: Telomere motif sequence (e.g., TTAGGG for vertebrates)
- `--work`: Number of chromosome ends to process in parallel (default: 1)

### 4. AutoGapfiller Mode
**Purpose**: Automated pipeline for whole-genome gap filling with parallel processing.

**Use Cases**:
- Filling all gaps in a draft genome automatically
- Large-scale gap filling projects
- High-throughput genome finishing

**Key Features**:
- Automatic gap detection from N regions
- Parallel gap filling with multiple workers
- Integrated visualization and statistics
- Five-step automated pipeline:
  1. Reads preprocessing (splitting and indexing)
  2. Gap detection (identifying N regions)
  3. Gap filling (parallel DEGAP execution)
  4. Genome integration (merging filled gaps)
  5. Visualization (HTML reports and statistics)

**Key Parameters**:
- `--genome`: Input genome with gaps (FASTA)
- `--work`: Number of parallel workers (default: 4)
- Supports all DEGAP filtering and optimization parameters

## Input Data Requirements

### Supported Read File Formats

**DEGAP v2** supports the following file formats for sequencing reads:

- `--hifi`: HiFi reads file(s)
  - Supported formats: `.fa`, `.fq`, `.fa.gz`, `.fq.gz`, `.fasta`, `.fastq`
  - Supports one or more files in a single command
  - Can be uncompressed or gzip-compressed
  
- `--ont`: ONT reads file(s)
  - Supported formats: `.fa`, `.fq`, `.fa.gz`, `.fq.gz`, `.fasta`, `.fastq`
  - Supports one or more files in a single command
  - Can be uncompressed or gzip-compressed

### Data Modes

**DEGAP v2** supports three data modes for maximum flexibility:

1. **HiFi-only mode**: Provide only `--hifi` parameter (one or more files)
   - Best for high-accuracy gap filling with HiFi reads
   - Lower computational requirements

2. **ONT-only mode**: Provide only `--ont` parameter (one or more files)
   - Suitable when only ONT reads are available
   - Good for long-range gap filling

3. **Mixed mode** (Recommended): Provide both `--hifi` and `--ont` parameters (each can contain one or more files)
   - Combines strengths of both technologies
   - Best results for complex gap regions
   - Recommended for most applications

**Note**: At least one of `--hifi` or `--ont` must be provided.

## Parameters

### Required Parameters

| Parameter | Description | Applicable Modes |
|-----------|-------------|------------------|
| `--mode` | Execution mode: `gapfiller`, `ctglinker`, or `telseeker` | All |
| `-o, --out` | Output directory | All |

### Input Data Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--hifi` | HiFi reads file(s), supports one or more files (.fa/.fq/.fa.gz/.fq.gz/.fasta/.fastq) | None |
| `--ont` | ONT reads file(s), supports one or more files (.fa/.fq/.fa.gz/.fq.gz/.fasta/.fastq) | None |
| `--seqleft` | Left flanking sequence (GapFiller) | None |
| `--seqright` | Right flanking sequence (GapFiller) | None |
| `--ctgseq` | Input contigs (CtgLinker) | None |
| `--genome` | Input genome/chromosomes (TelSeeker) | None |
| `--motif` | Telomere motif sequence (TelSeeker) | None |
| `-e, --target_ends` | Target chromosome ends for TelSeeker, e.g. `Chr01.L Chr01.R`, or a text file with one target per line | Required |
| `--telo-read-stringency` | Telomeric read discovery strictness: `strict`, `normal`, or `relaxed` | `normal` |

**Note**: At least one of `--hifi` or `--ont` must be provided.

For TelSeeker, each target end must use the format `<chromosome_id>.L` or
`<chromosome_id>.R`. To discover candidate target ends first, run:

```bash
python bin/TelSeekerCheck.py --genome genome.fasta --motif TTAGGG --out check_out
```

Then pass the generated target list:

```bash
-e check_out/genome.telomere.check/need_extension_chr_end.txt
```

Manually review the `TelSeekerCheck.py` result before entering TelSeeker. The
check step generates whole-genome motif distribution plots with the same 1 kb
plotting window used by the final TelSeeker report:
`all_chromosomes_combined.png` and one `<chromosome>_telomere_motif.png` image
per chromosome. Use these plots to decide which ends should be passed to
TelSeeker. The legacy helper files `genome.telomere.check.csv`,
`need_extension_chr_end.txt`, and `uncertain_chr_end.txt` are still written for
compatibility, but should be manually reviewed before choosing targets. If you
already have clear target ends, skip `TelSeekerCheck.py` and pass them directly
with `-e Chr01.L Chr01.R` or a target file.

### Filtering Parameters

| Parameter | Description | Default | Recommended |
|-----------|-------------|---------|-------------|
| `--filterDepthHifi` | HiFi depth filtering threshold | None | 0.3 |
| `--filterDepthOnt` | ONT depth filtering threshold | None | 0.3 |
| `--kmer_filter` | Enable k-mer filtering | False | True for large datasets |
| `--kmer_size` | K-mer size for filtering | 41 | 41 |
| `--kmer_num` | Number of k-mers to extract | 20 | 20 |

**Depth Filtering Explanation**:
- Filters reads based on mapping depth to reference
- Threshold `0.3` means: keep reads with depth outside `[0.3×mode, 1.7×mode]`
- Reduces redundant high-coverage reads and low-quality reads
- Applied once at the beginning; all subsequent steps use filtered data

**K-mer Filtering Explanation**:

The `--kmer_filter` parameter is a switch to enable k-mer based read filtering, which provides significant performance improvements for large datasets:

**How it works**:
- Extracts representative k-mers from target sequence
- Filters reads containing these k-mers
- Significantly reduces data size (typically 50-70% reduction)
- Maintains reads relevant to the gap-filling task

**When to enable** (add `--kmer_filter`):
- Processing large-scale datasets (reads >30GB)
- Significantly reduces memory usage and processing time
- Recommended for most production workflows

**When to disable** (omit `--kmer_filter`):
- For highest accuracy requirements
- Small to medium datasets (reads <30GB)
- When computational resources are sufficient
- When every read might be critical

**Performance Impact**:
- Can reduce processing time by 50-80% on large datasets
- Memory usage typically reduced by 50-70%
- Minimal impact on gap-filling accuracy (usually <1% difference)

### Performance Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-t, --thread` | Number of threads | 20 |
| `--remove` | Cleanup level (1/2/3) | 2 |

### Advanced Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--flag` | Elongation direction (GapFiller) | `left` |
| `--edge` | Edge length for error detection | 500 |
| `--MaximumExtensionLength` | Maximum extension length | None |
| `--seedLen` | Seed length for alignment | Auto-calculated |
| `--work` | Parallel chromosome ends (TelSeeker) | 1 |

### AutoGapfiller Specific Parameters

**Required Parameters:**

| Parameter | Description |
|-----------|-------------|
| `--genome` | Input genome FASTA file (with N gaps) |
| `--out` | Output directory |

**Reads Files (at least one required):**

| Parameter | Description |
|-----------|-------------|
| `--hifi` | HiFi reads file(s) (FASTA/FASTQ, .gz supported; multiple files supported) |
| `--ont` | ONT reads file(s) (FASTA/FASTQ, .gz supported; multiple files supported) |

**Common Parameters:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--work` | Number of parallel workers | 4 |
| `--kmer_filter` | Enable k-mer filtering | False |

**Advanced Parameters (inherited from DEGAP):**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-t, --thread` | Threads for parallel processing (controls both threading and parallel jobs) | 20 |
| `-ks, --kmer_size` | K-mer size | 41 |
| `-kn, --kmer_num` | Number of k-mers | 20 |
| `--MaximumExtensionLength` | Maximum extension length | 1000000 |
| `--MaximumExtensionRound` | Maximum extension rounds | None |
| `--filterDepthHifi` | HiFi depth filtering threshold | None |
| `--filterDepthOnt` | ONT depth filtering threshold | None |
| `--edge` | Edge length for misassembly detection | 500 |
| `--remove` | File cleanup level (1/2/3) | 2 |

## Usage Examples

### Quick Start (Recommended Configurations)

These examples demonstrate the recommended configurations for each mode using mixed data (HiFi + ONT) with k-mer filtering enabled:

**GapFiller Mode - Mixed Mode with K-mer Filtering:**
```bash
python bin/DEGAP.py --mode gapfiller \
    --seqleft ./path/gapLeftSequence.fasta \
    --seqright ./path/gapRightSequence.fasta \
    --hifi ./path/HiFiReads.fq.gz \
    --ont ./path/ONTReads.fq.gz \
    --out ./path/Output \
    --flag left \
    --kmer_filter \
    -t 20
```

**Multi-file input (Issue #8 style)**
```bash
python bin/DEGAP.py --mode ctglinker \
    --ctgseq ./path/contigs.fasta \
    --hifi ./path/HIFI1.fq.gz ./path/HIFI2.fq.gz ./path/HIFI3.fq.gz \
    --ont ./path/ONT1.fq.gz ./path/ONT2.fq.gz \
    --out ./path/Output \
    --kmer_filter \
    -t 20
```

**CtgLinker Mode - Mixed Mode with K-mer Filtering:**
```bash
python bin/DEGAP.py --mode ctglinker \
    --ctgseq ./path/contigs.fasta \
    --hifi ./path/HiFiReads.fq.gz \
    --ont ./path/ONTReads.fq.gz \
    --out ./path/Output \
    --kmer_filter \
    -t 20
```

**TelSeeker Mode - Mixed Mode with K-mer Filtering:**
```bash
python bin/DEGAP.py --mode telseeker \
    --genome ./path/genome.fasta \
    --motif TTAGGG \
    -e ./path/target_ends.txt \
    --hifi ./path/HiFiReads.fq.gz \
    --ont ./path/ONTReads.fq.gz \
    --out ./path/Output \
    --work 4 \
    --telo-read-stringency normal \
    --kmer_filter \
    -t 20
```

**AutoGapfiller Mode - Mixed Mode with K-mer Filtering:**
```bash
python bin/AutoGapfiller.py \
    --genome ./path/genome.fasta \
    --hifi ./path/HiFiReads.fq.gz \
    --ont ./path/ONTReads.fq.gz \
    --out ./path/Output \
    --work 4 \
    --kmer_filter \
    -t 20
```

### Additional Examples

Below are additional examples demonstrating various configurations and use cases:

### Example 1: Basic GapFiller (HiFi only)
```bash
python bin/DEGAP.py --mode gapfiller \
    --seqleft gap_left.fasta \
    --seqright gap_right.fasta \
    --hifi hifi_reads.fq \
    -o ./output_gapfiller \
    --flag left \
    -t 20
```

### Example 2: CtgLinker with Depth Filtering
```bash
python bin/DEGAP.py --mode ctglinker \
    --ctgseq contigs.fasta \
    --hifi hifi_reads.fq \
    -o ./output_ctglinker \
    --filterDepthHifi 0.3 \
    -t 20
```

### Example 3: CtgLinker with K-mer Filtering
```bash
python bin/DEGAP.py --mode ctglinker \
    --ctgseq contigs.fasta \
    --hifi hifi_reads.fq \
    -o ./output_ctglinker \
    --filterDepthHifi 0.3 \
    --kmer_filter \
    --kmer_size 41 \
    --kmer_num 20 \
    -t 20 \
```

### Example 4: Mixed Mode (HiFi + ONT)
```bash
python bin/DEGAP.py --mode ctglinker \
    --ctgseq contigs.fasta \
    --hifi hifi_reads.fq \
    --ont ont_reads.fq \
    -o ./output_mixed \
    --filterDepthHifi 0.3 \
    --filterDepthOnt 0.3 \
    --kmer_filter \
    -t 20
```

### Example 5: ONT-only Mode
```bash
python bin/DEGAP.py --mode ctglinker \
    --ctgseq contigs.fasta \
    --ont ont_reads.fq \
    -o ./output_ont \
    --filterDepthOnt 0.3 \
    -t 20
```

### Example 6: TelSeeker Mode
```bash
python bin/DEGAP.py --mode telseeker \
    --genome genome.fasta \
    --motif TTAGGG \
    -e Chr01.L Chr01.R \
    --hifi hifi_reads.fq \
    -o ./output_telseeker \
    -t 20
```

### Example 7: TelSeeker Mode (Recommended with Filtering)
```bash
python bin/DEGAP.py --mode telseeker \
    --genome genome.fasta \
    --motif TTAGGG \
    -e target_ends.txt \
    --hifi hifi_reads.fq \
    --ont ont_reads.fq \
    -o ./output_telseeker \
    --work 2 \
    --filterDepthHifi 0.3 \
    --filterDepthOnt 0.3 \
    --kmer_filter \
    -t 40
```

### Example 8: Depth Filtering Only (No Gap Filling)
```bash
python bin/DEGAP.py --mode ctglinker \
    --ctgseq contigs.fasta \
    --hifi hifi_reads.fq \
    -o ./output_filter \
    --filterDepthOnly 0.3 \
    -t 30
```

### Example 9: AutoGapfiller - Basic Usage (HiFi only)
```bash
python bin/AutoGapfiller.py \
    --genome genome_with_gaps.fasta \
    --hifi hifi_reads.fq \
    --out ./output_autogap \
    --work 4
```

### Example 10: AutoGapfiller - Recommended (HiFi + ONT with k-mer filter)
```bash
python bin/AutoGapfiller.py \
    --genome genome_with_gaps.fasta \
    --hifi hifi_reads.fq \
    --ont ont_reads.fq \
    --out ./output_autogap \
    --work 4 \
    --kmer_filter
```

### Example 11: AutoGapfiller - Advanced (Custom Parameters)
```bash
python bin/AutoGapfiller.py \
    --genome genome_with_gaps.fasta \
    --hifi hifi_reads.fq \
    --ont ont_reads.fq \
    --out ./output_autogap \
    --work 10 \
    -t 30 \
    --kmer_filter \
    -ks 51 \
    -kn 30 \
    --MaximumExtensionLength 500000 \
    --filterDepthHifi 0.3 \
    --filterDepthOnt 0.3 \
    -j 100
```

## Output Files

### Main Output Files

| File | Description | Mode |
|------|-------------|------|
| `DG.scaffold.fa` | Final scaffold sequences | CtgLinker, TelSeeker |
| `DG.scaffold.log` | Assembly statistics (N50, etc.) | CtgLinker, TelSeeker |
| `DG.project.fa` | Project sequences | CtgLinker, TelSeeker |
| `Genome.inputCtg.fa` | Processed input contigs | All |
| `linkedSequence.*.fa` | Linked sequences per round | GapFiller |
| `genome.filled.fasta` | Final filled genome | AutoGapfiller |
| `autogapfiller_summary.json` | Overall pipeline summary | AutoGapfiller |
| `gap_filling_report.html` | Interactive visualization report | AutoGapfiller |

### Intermediate Files

| File/Directory | Description |
|----------------|-------------|
| `Genome.inputCtg.usefullReads.fastq` | Filtered HiFi reads (if depth filtering enabled) |
| `Genome.inputCtg.usefullReads.ont.fastq` | Filtered ONT reads (if depth filtering enabled) |
| `hifi_reads_part/` | Split HiFi reads for parallel processing |
| `ont_reads_part/` | Split ONT reads for parallel processing |
| `project/` | Per-round intermediate results |
| `selectedReads.idx` | Index of filtered reads |

### K-mer Filtering Output (if enabled)

| File/Directory | Description |
|----------------|-------------|
| `temp_kmer_filter_1/` | First-level k-mer filtering results |
| `temp_kmer_filter_2/` | Second-level k-mer filtering results |
| `hifi.kmer.Reads.*.fa` | K-mer filtered HiFi reads |
| `ont.kmer.Reads.*.fa` | K-mer filtered ONT reads |
| `kmer_filtering.log.json` | K-mer filtering statistics |

### GapFiller Output Structure

```
output_directory/
├── <name>.final.fa                 # Final filled sequence
├── process.log                     # Detailed processing log
├── process.summary                 # Summary of each round
├── process.agp                     # AGP format assembly info
├── checkpoint.info                 # Checkpoint for resume
├── process/                        # Round-by-round results
│   ├── round1/
│   │   ├── linkedSequence.fasta   # Linked sequence for this round
│   │   ├── outputExtensionSequence.fasta
│   │   ├── potentialExtensionReads.*.sam
│   │   └── ...
│   ├── round2/
│   └── ...
└── USED_DIRECT_CONNECTION          # Marker if direct connection was used
```

**Key Output Files:**
- `<name>.final.fa`: Final gap-filled sequence (main result)
- `process.log`: Detailed log including gap length, alignment info, extension details
- `process.summary`: Tab-separated summary of each extension round
- `linkedSequence.fasta`: Successfully linked sequence (if gap closed)

### CtgLinker Output Structure

```
output_directory/
├── DG.scaffold.fa                  # Final scaffold sequences (main result)
├── DG.scaffold.log                 # Assembly statistics (N50, etc.)
├── DG.scaffold.agp                 # AGP format scaffold assembly
├── DG.project.fa                   # All project sequences
├── DG.agp.path.txt                 # AGP path information
├── Genome.inputCtg.fa              # Processed input contigs
├── Genome.inputCtg.log             # Input contigs statistics
├── used.txt                        # Contig usage tracking
├── unplaced.txt                    # Unplaced contigs list
├── project/                        # Project-by-project results
│   ├── DG1/                       # Project 1
│   │   ├── DG1.fa                 # Project sequence
│   │   ├── DG1.log                # Project log
│   │   ├── DG1.agp                # Project AGP
│   │   ├── DG1.agp.path.txt       # AGP path
│   │   ├── DG1.used.txt           # Used contigs
│   │   ├── DG1.input.log          # Input log
│   │   └── data/                  # Project data
│   │       ├── terminalCtg.fa     # Terminal contigs
│   │       ├── left/              # Left extension
│   │       │   └── process/       # Extension rounds
│   │       └── right/             # Right extension
│   │           └── process/       # Extension rounds
│   ├── DG2/                       # Project 2
│   └── ...
├── hifi_reads_part/                # Split HiFi reads (if preprocessing done)
├── ont_reads_part/                 # Split ONT reads (if preprocessing done)
└── Genome.inputCtg.usefullReads.*  # Filtered reads (if depth filtering enabled)
```

**Key Output Files:**
- `DG.scaffold.fa`: Final scaffold sequences (main result)
- `DG.scaffold.log`: Assembly statistics including N50, total length, number of scaffolds
- `DG.scaffold.agp`: AGP format file describing scaffold assembly
- `DG.project.fa`: All intermediate project sequences
- `project/DG*/DG*.fa`: Individual project sequences for each linking iteration

### TelSeeker Output Structure

```
output_directory/
├── genome.telomere.check/          # Target end list for TelSeekerPart2
│   └── need_extension_chr_end.txt # Target ends provided by -e/--target_ends
├── part1.telo.reads/               # Telomeric reads extraction
│   ├── Global.left.telo.reads.fa  # All left telomeric reads
│   ├── Global.right.telo.reads.fa # All right telomeric reads
│   ├── left.telo.reads.fa         # Filtered left reads
│   ├── right.telo.reads.fa        # Filtered right reads
│   └── part1.log                  # Extraction log
├── part2.chr.end.job/              # Chromosome end extension jobs
│   ├── Chr01.Left/                # Job for Chr01 left end
│   │   ├── linker.fa              # Final linker sequence
│   │   ├── direct.check/          # Direct connection check
│   │   ├── extension/             # Extension results
│   │   │   └── round*/            # Round-by-round extension
│   │   └── workspace/             # Working files
│   ├── Chr01.Right/               # Job for Chr01 right end
│   └── ...
├── part3.integration.results/      # Integration and summary
│   ├── check_part2_jobs.csv       # Summary of all chromosome ends
│   ├── final.genome.fa            # Final integrated genome
│   └── final.genome.telomere.check/  # Final genome telomere check
│       ├── genome.telomere.check.csv
│       └── ...
├── visual.report/                  # Visualization reports
│   ├── Global.report.html         # Global overview report
│   ├── Chr01.Left.report.html     # Individual chromosome end reports
│   └── ...
├── hifi_reads_part/                # Split HiFi reads (if preprocessing done)
├── ont_reads_part/                 # Split ONT reads (if preprocessing done)
└── Genome.inputCtg.usefullReads.*  # Filtered reads (if depth filtering enabled)
```

**Key Output Files:**
- `part3.integration.results/final.genome.fa`: Final genome with extended telomeres
- `genome.telomere.check/need_extension_chr_end.txt`: Target chromosome ends used by TelSeeker
- `part3.integration.results/final.genome.telomere.check/uncertain_chr_end.txt`: Final-genome ends needing manual telomere review, if any
- `part2.chr.end.job/<Chr>.<End>/linker.fa`: Extended sequence for each chromosome end
- `visual.report/Global.report.html`: Interactive visualization of all results

### AutoGapfiller Output Structure

```
output_directory/
├── autogapfiller_summary.json      # Overall pipeline summary
├── 01.reads_preprocessor/          # Reads preprocessing
│   ├── processed_reads/           # Processed and converted reads
│   │   ├── HiFi.reads.fa         # HiFi reads (FASTA format)
│   │   ├── ONT.reads.fa          # ONT reads (FASTA format)
│   │   ├── HiFi.reads.stat       # HiFi reads statistics
│   │   └── ONT.reads.stat        # ONT reads statistics
│   ├── hifi_reads_part/           # Split HiFi reads for parallel processing
│   ├── ont_reads_part/            # Split ONT reads for parallel processing
│   └── preprocessing_config.json  # Preprocessing configuration
├── 02.gap_detection/               # Gap detection
│   ├── gaps.bed                   # Gap coordinates (BED format)
│   ├── gap_sequences/             # Extracted gap sequences
│   │   ├── left/                 # Left flanking sequences
│   │   ├── right/                # Right flanking sequences
│   │   └── gap_info.txt          # Gap information
│   └── gap_detection_config.json  # Gap detection configuration
├── 03.gap_filling/                 # Gap filling results
│   ├── gap_filling_results.json   # Gap filling statistics
│   └── jobs/                       # Individual gap filling jobs
│       ├── gap_001/               # Job for gap 1
│       │   ├── <gap_name>.final.fa          # Final filled sequence
│       │   ├── process.log                   # Processing log
│       │   ├── process.summary               # Summary
│       │   └── process/                      # Round-by-round results
│       ├── gap_002/               # Job for gap 2
│       └── ...
├── 04.genome_integration/          # Genome integration
│   ├── genome.filled.fasta        # Final filled genome (main result)
│   └── integration_report.json    # Integration statistics
└── 05.visualization/               # Visualization
    └── gap_filling_report.html    # Interactive HTML report
```

**Key Output Files:**
- `04.genome_integration/genome.filled.fasta`: Final genome with filled gaps (main result)
- `05.visualization/gap_filling_report.html`: Interactive visualization report
- `03.gap_filling/gap_filling_results.json`: Detailed gap filling statistics
- `autogapfiller_summary.json`: Overall pipeline summary

## Performance Optimization

### Memory Optimization
- Use `--kmer_filter` for large datasets (>30GB) to reduce memory usage by 50-70%
- Adjust `-t` (threads) parameter based on available CPU cores
- Use depth filtering (`--filterDepthHifi` or `--filterDepthOnt`) to reduce data size before processing
- Split large read files if memory constraints are tight

### Speed Optimization
- **Enable k-mer filtering** (`--kmer_filter`) for datasets >30GB - provides 50-80% speed improvement
- Increase `-t` (threads) based on CPU cores (recommended: 20-40 threads)
- Use SSD storage for intermediate files to reduce I/O bottleneck
- Pre-filter low-quality reads before running DEGAP
- Use mixed mode (HiFi + ONT) for best results with reasonable resource usage

### Recommended Settings by Dataset Size

| Genome Size | Reads Size | Recommended Configuration | Example Command Additions |
|-------------|------------|---------------------------|---------------------------|
| Small (<500MB) | <30GB | Basic mode without k-mer filtering | `-t 20` |
| Medium (500MB-2GB) | 30-100GB | Enable k-mer filtering | `--kmer_filter -t 20` |
| Large (>2GB) | >100GB | K-mer filtering + parallel workers | `--kmer_filter --work 10 -t 40` |

**General Recommendations**:
- **Small datasets (<30GB)**: K-mer filtering optional, focus on accuracy
- **Medium datasets (30-100GB)**: Enable `--kmer_filter` for better performance
- **Large datasets (>100GB)**: Always use `--kmer_filter` + parallel processing

### AutoGapfiller Performance Tips

**Worker Count (`--work`)**:
- Small genomes (<1GB, <100 gaps): `--work 4`
- Medium genomes (1-3GB, 100-500 gaps): `--work 10`
- Large genomes (>3GB, >500 gaps): `--work 20`

**Memory Considerations**:

- Each worker runs a DEGAP instance
- Estimate: `Memory = workers × (2-4GB per worker)`
- Example: 10 workers ≈ 20-40GB RAM

**Disk Space**:
- Temporary files can be large
- Estimate: `Disk = genome_size × 10-20`
- Use `--remove 3` to minimize disk usage

### TelSeeker Performance Tips

**Worker Count (`--work`)**:
- Small genomes (<10 chromosomes): `--work 1-2`
- Medium genomes (10-30 chromosomes): `--work 4-8`
- Large genomes (>30 chromosomes): `--work 10-20`

**Memory Considerations**:
- Each worker processes one chromosome end
- Estimate: `Memory = workers × (2-4GB per worker)`
- Example: 4 workers ≈ 8-16GB RAM

**Processing Strategy**:
- Run `TelSeekerCheck.py` separately if you need target discovery from telomere motif windows
- Manually confirm the `TelSeekerCheck.py` output before using it as TelSeeker input
- Review `uncertain_chr_end.txt` as ambiguous evidence, not as an automatic TelSeeker target list
- Skip `TelSeekerCheck.py` when target ends are already known
- Pass target ends explicitly with `-e Chr01.L Chr01.R` or `-e target_ends.txt`
- Use `--telo-read-stringency strict` for legacy marker+TRC read discovery
- Use `--telo-read-stringency relaxed` when Step 1 finds too few or no telo reads
- Use `--kmer_filter` for large genomes to reduce processing time

## Troubleshooting

### Common Issues

**1. Out of Memory Error**
- Solution: Enable `--kmer_filter` and depth filtering
- Reduce `-j` parameter
- Process on a machine with more RAM

**2. No Extension Reads Found**
- Check if input sequences are correct
- Try adjusting `--seedLen` parameter
- Verify read quality and coverage

**3. K-mer Filtering Fails**
- Ensure GNU Parallel is installed
- Check if seqkit is in PATH
- Verify sufficient disk space for temporary files

**4. Slow Performance**
- Enable k-mer filtering
- Increase thread count (`-t`)
- Use faster storage (SSD)

**5. Depth Filtering Removes Too Many Reads**
- Adjust threshold (try 0.2 or 0.4 instead of 0.3)
- Check read coverage distribution
- Verify reference sequence quality

**6. AutoGapfiller Worker Failures**
- Check individual job logs in `03.gap_filling/jobs/gap_XXX/`
- Reduce worker count (`--work`) if memory issues
- Verify gap sequences are valid
- Check gap filling results JSON for failed jobs

**7. AutoGapfiller No Gaps Detected**
- Verify input genome contains N regions
- Check gap detection log in `02.gap_detection/`
- Ensure gaps are at least 1bp long
- Verify FASTA format is correct

**8. TelSeeker No Chromosome Ends Extended**
- Verify `-e/--target_ends` contains the intended chromosome ends
- If targets came from `TelSeekerCheck.py`, manually confirm the check results before rerunning
- If using a target file, verify it has one non-empty `Chr.L` or `Chr.R` target per line
- If the log reports `No telomeric reads found after Step 1`, or `part1.telo.reads/Global.left.telo.reads.fa` and `Global.right.telo.reads.fa` are empty, rerun with `--telo-read-stringency relaxed`
- Check `--motif` parameter matches your species (e.g., TTAGGG for vertebrates)

Example rerun with relaxed telo-read discovery:

```bash
python bin/DEGAP.py --mode telseeker \
    --genome genome.fasta \
    --motif TTAGGG \
    -e target_ends.txt \
    --hifi reads.fa \
    --out output_telseeker_relaxed \
    --telo-read-stringency relaxed
```

**9. TelSeeker Motif Not Found**
- Verify motif sequence is correct for your species
- Use uppercase letters only (A/T/C/G)
- Common motifs: TTAGGG (vertebrates), TTTAGGG (plants), TTAGG (insects)
- Check if genome already has complete telomeres

### Log Files
Check the following log files for debugging:

**GapFiller mode:**
- `process.log`: Detailed processing log with gap length and alignment info
- `process.summary`: Summary of each extension round
- `checkpoint.info`: Checkpoint information for resume

**CtgLinker mode:**
- `project/round*/process.log`: Per-round processing logs
- `DG.scaffold.log`: Final scaffold statistics

**TelSeeker mode:**
- `genome.telomere.check/need_extension_chr_end.txt`: Target ends loaded from `-e/--target_ends`
- `part3.integration.results/final.genome.telomere.check/uncertain_chr_end.txt`: Ambiguous final telomere-check calls
- `part1.telo.reads/part1.log`: Telomeric reads extraction log
- `part2.chr.end.job/<Chr>.<End>/extension/process.log`: Extension log for each chromosome end
- `part3.integration.results/check_part2_jobs.csv`: Summary of all chromosome ends
- `visual.report/Global.report.html`: Interactive visualization report

**AutoGapfiller mode:**
- `autogapfiller_summary.json`: Overall pipeline statistics
- `03.gap_filling/gap_filling_results.json`: Gap filling statistics
- `03.gap_filling/jobs/gap_XXX/process.log`: Individual gap filling logs
- `05.visualization/gap_filling_report.html`: Interactive visualization report

## Citation

Huang, Y. et al. DEGAP: Dynamic Elongation of a Genome Assembly Path. Briefings in Bioinformatics, Volume 25, Issue 3, May 2024, bbae194. [DOI: 10.1093/bib/bbae194](https://doi.org/10.1093/bib/bbae194)

## Contact

For questions, bug reports, or feature requests:
- Email: <jzhang@mail.hzau.edu.cn>

## License

Please refer to the LICENSE file for licensing information.
