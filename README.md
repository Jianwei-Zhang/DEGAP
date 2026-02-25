![DEGAP](./DEGAP.png)

## Dynamic Elongation of a Genome Assembly Path

**DEGAP degaps gaps!**

**DEGAP v2** is an enhanced gap-filling software that resolves gap regions by utilizing high-fidelity (HiFi) and/or Oxford Nanopore (ONT) long reads. It supports multiple modes for different assembly scenarios.

## System Requirements

**DEGAP v2** has been developed with Python under a Linux environment. The following packages/tools are required (versions aligned with `environment.yml`):

**Python Packages:**
- Python (3.10)
- Biopython
- Pysam
- NumPy
- Matplotlib
- SciPy

**External Tools:**
- Minimap2 (2.28)
- Hifiasm (0.25.0)
- SAMtools (1.22.1)
- Seqkit (2.8.0)
- MUMmer4
- GNU Parallel
- Ripgrep (15.1.0)

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

## Tutorial

**DEGAP v2** provides 4 modes (GapFiller, CtgLinker, TelSeeker, AutoGapfiller) to fill gaps and extend sequences.

### Input Data Requirements

**Supported Read Files:**
- `--hifi`: HiFi reads file(s), supports one or more files (formats: `.fa`, `.fq`, `.fa.gz`, `.fq.gz`, `.fasta`, `.fastq`)
- `--ont`: ONT reads file(s), supports one or more files (formats: `.fa`, `.fq`, `.fa.gz`, `.fq.gz`, `.fasta`, `.fastq`)

**Note:** At least one of `--hifi` or `--ont` must be provided.
You can provide multiple files directly, for example:
`--hifi HIFI1.fq.gz HIFI2.fq.gz HIFI3.fq.gz`

### K-mer Filtering Option

The `--kmer_filter` parameter is a switch to enable k-mer based read filtering:

- **When to enable** (add `--kmer_filter`):
  - Processing large-scale datasets (reads.fa >20GB)
  - Significantly reduces memory usage and processing time
- **When to disable** (omit `--kmer_filter`):
  - For highest accuracy requirements
  - Small to medium datasets (reads.fa <20GB)
  - When computational resources are sufficient

**Note:** K-mer filtering provides a good balance between speed and accuracy for most applications.

### GapFiller Mode
Designed for filling a specific gap when the exact left and right sequences are known.

```bash
python bin/DEGAP.py --mode gapfiller \
    --seqleft ./path/gapLeftSequence.fasta \
    --seqright ./path/gapRightSequence.fasta \
    --hifi ./path/HiFiReads.fa \
    --ont ./path/ONTReads.fa \
    --out ./path/Output \
    --flag left \
    --kmer_filter \
    -t 20
```

### CtgLinker Mode
Designed for assemblies with gaps. Links contigs and fills gaps by elongating edges.

```bash
python bin/DEGAP.py --mode ctglinker \
    --ctgseq ./path/contigs.fasta \
    --hifi ./path/HiFiReads.fa \
    --ont ./path/ONTReads.fa \
    --out ./path/Output \
    --kmer_filter \
    -t 20
```

### TelSeeker Mode
Designed for telomere extension and chromosome-level assembly.

```bash
python bin/DEGAP.py --mode telseeker \
    --genome ./path/genome.fasta \
    --motif TTAGGG \
    --hifi ./path/HiFiReads.fa \
    --ont ./path/ONTReads.fa \
    --out ./path/Output \
    --work 4 \
    --kmer_filter \
    -t 20
```

### AutoGapfiller Mode
Automated pipeline for whole-genome gap filling. Automatically detects gaps and fills them in parallel.

```bash
python bin/AutoGapfiller.py \
    --genome ./path/genome.fasta \
    --hifi ./path/HiFiReads.fa \
    --ont ./path/ONTReads.fa \
    --out ./path/Output \
    --work 4 \
    --kmer_filter \
    -t 20
```

**Note:** You can add the `--setup_only` parameter to only create the gap-filling job directories without executing them. This is useful for inspecting the tasks before running or for manual execution control.

## Key Features in v2

- **Multi-data support**: HiFi-only, ONT-only, or Mixed (HiFi+ONT) modes
- **K-mer filtering**: `--kmer_filter` for faster processing with large datasets
- **Enhanced performance**: Optimized parallel processing and memory management

## Quick Start

For detailed usage, parameter descriptions, and advanced features, please refer to [Users' Guide](HOWTO.md).

## How to Cite?

Huang, Y. et al. DEGAP: Dynamic Elongation of a Genome Assembly Path. Briefings in Bioinformatics, Volume 25, Issue 3, May 2024, bbae194. [DOI: 10.1093/bib/bbae194](https://doi.org/10.1093/bib/bbae194)

## Contact

For comments and suggestions, please contact <jzhang@mail.hzau.edu.cn>.
