![DEGAP](./DEGAP.png)

# DEGAP 2.0 — Dynamic Elongation of a Genome Assembly Path

**DEGAP 2.0** resolves genome assembly gaps, links contigs, and extends chromosome ends using HiFi and/or Oxford Nanopore (ONT) long reads.

## Installation

DEGAP 2.0 runs on Linux. The supplied environment pins Python 3.10; packages and external tools are declared in [environment.yml](environment.yml):

- **Python packages**: Python 3.10, Biopython, Pysam, NumPy, Matplotlib, SciPy
- **External tools**: Minimap2 (2.28), Hifiasm (0.25.0), SAMtools (1.22.1), Seqkit (2.8.0), MUMmer4, GNU Parallel, Ripgrep (15.1.0)

Recommended setup with [Micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html):

```bash
micromamba env create -f environment.yml
micromamba activate degap_env
```

To update an existing environment after `environment.yml` changes:

```bash
micromamba env update -f environment.yml --prune
```

## Choosing a Mode

| Mode | Command Entry | Purpose | Key Inputs |
| --- | --- | --- | --- |
| **GapFiller** | `python bin/DEGAP.py --mode gapfiller` | Fill one gap between known left and right flanks using bidirectional extension families | `--seqleft`, `--seqright`, `--hifi`/`--ont`, `-o/--out` |
| **CtgLinker** | `python bin/DEGAP.py --mode ctglinker` | Link contigs and close gaps from contig edges using primary and forward/reverse candidate families | `--ctgseq`, `--hifi`/`--ont`, `-o/--out` |
| **TelSeeker** | `python bin/DEGAP.py --mode telseeker` | Single-ended extension of specified chromosome ends (`<chr>.L` / `<chr>.R`) toward telomeric reads | `--genome`, `-m/--motif`, `-e/--target_ends`, `--hifi`/`--ont`, `-o/--out` |
| **TelSeeker Ctg** | `python bin/DEGAP.py --mode telseeker_ctg` | Single-ended extension of one or more single-record contig endpoints (`L` or `R`) with shared read/telomere caching | `--ctg <fasta> <L\|R>`, `-m/--motif`, `--hifi`/`--ont`, `-o/--out` |
| **AutoGapfiller** | `python bin/AutoGapfiller.py` | Detect genome N gaps and run one bidirectional gap-filling job per gap | `--genome`, `--hifi`/`--ont`, `-o/--output` (optional `--setup_only`) |

Standalone chromosome-end review plots and 2 kb end FASTAs can be generated beforehand with `python bin/TelSeekerCheck.py --genome <genome.fa> --motif <MOTIF> --out <dir>`.

## Quick Start

### 1. Fill a specific gap (GapFiller)

```bash
python bin/DEGAP.py --mode gapfiller \
    --seqleft ./left_flank.fa \
    --seqright ./right_flank.fa \
    --hifi ./hifi_reads.fq.gz \
    --ont ./ont_reads.fq.gz \
    --out ./gap_out \
    --MaximumExtensionRound 25 \
    -t 20
```

### 2. Fill all detected N gaps in a genome (AutoGapfiller)

Note that `AutoGapfiller.py` uses `--output` / `-o` (not `--out`):

```bash
python bin/AutoGapfiller.py \
    --genome ./genome_with_gaps.fa \
    --hifi ./hifi_reads.fq.gz \
    --ont ./ont_reads.fq.gz \
    --output ./autogap_out \
    --work 4 \
    --MaximumExtensionRound 25 \
    -t 20
```

Add `--setup_only` to prepare job directories and `run.sh` scripts without executing them.

### 3. Extend selected chromosome ends (TelSeeker)

Optional review step to inspect chromosome ends before choosing targets:

```bash
python bin/TelSeekerCheck.py \
    --genome ./genome.fa \
    --motif TTAGGG \
    --out ./telomere_check
```

Then run TelSeeker with explicit target ends (`-e Chr01.L Chr01.R` or `-e target_ends.txt`):

```bash
python bin/DEGAP.py --mode telseeker \
    --genome ./genome.fa \
    --motif TTAGGG \
    -e Chr01.L Chr01.R \
    --hifi ./hifi_reads.fq.gz \
    --ont ./ont_reads.fq.gz \
    --out ./telseeker_out \
    --work 2 \
    --MaximumExtensionRound 25 \
    -t 20
```

## Before You Run

- **Use original FASTQ for ONT.** ONT-only assembly needs genuine quality scores, including rounds where mixed input falls back to ONT-only. HiFi accepts FASTA or FASTQ. ONT FASTA works only as supplementary `--ul` input while HiFi remains available; it is not a replacement for ONT FASTQ.
- **Both sides share one budget.** GapFiller alternates left and right extension, checking each new sequence against the opposite side's history. If one side stops normally, the other can continue. Round and length limits count the total across both sides, not a separate allowance for each. Limits are unlimited by default; set them for bounded runs.
- **Completion is not closure.** Check `status: "closed"` in `family-result.json` or the HTML report. An unclosed run exports separate family sequences. A closed result still needs sequence validation; a zero process exit alone proves neither closure nor biological accuracy.
- **Resume without deleting evidence.** Rerun the same command, optionally with `--resume_auto`. Keep inputs unchanged, including ONT quality scores. Old single-sided directories and `--resume ROUND` are unsupported; use a new output directory. Family histories and recovery files are retained even with `--remove 1|2`.
- **Omit `--flag`.** It is accepted only for compatibility, prints a notice, and has no effect.

See [HOWTO.md](HOWTO.md) for budget boundaries, historical-pair selection, ONT quality recovery, and tool-error handling.

## Documentation

- **Practical Users' Guide**: [HOWTO.md](HOWTO.md) — full parameter tables, mode workflows, output layouts, and troubleshooting.
- **Design Notes**: [docs/design/bidirectional-extension.md](docs/design/bidirectional-extension.md), [docs/design/report-presentation.md](docs/design/report-presentation.md)
- **New-data validation**: [SSH 180: two masked HiFi gaps closed](docs/validation/2026-11-22-ssh180-new-dataset.md) — measured gap-only comparisons and ONT suitability limits.
- **Earlier Validation Records**: [docs/validation/2026-11-21-ssh180-real-ont.md](docs/validation/2026-11-21-ssh180-real-ont.md), [docs/validation/2026-11-21-ont-quality-repair.md](docs/validation/2026-11-21-ont-quality-repair.md)
- **Historical v1 Guide**: [v1/README.md](v1/README.md)

## Citation

Huang, Y. et al. DEGAP: Dynamic Elongation of a Genome Assembly Path. *Briefings in Bioinformatics*, Volume 25, Issue 3, May 2024, bbae194. [DOI: 10.1093/bib/bbae194](https://doi.org/10.1093/bib/bbae194)

## Contact

For comments, bug reports, or questions, please contact <jzhang@mail.hzau.edu.cn>.
