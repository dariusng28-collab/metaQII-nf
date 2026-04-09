# metaQII-nf

> A reproducible Nextflow pipeline for QIIME2-based amplicon sequencing analysis (16S & ITS).

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04-brightgreen)](https://www.nextflow.io/)
[![QIIME2](https://img.shields.io/badge/QIIME2-2026.1-blue)](https://qiime2.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20WSL2-lightgrey)]()

---

## Overview

`metaQII-nf` is an end-to-end amplicon metagenomics pipeline built on QIIME2. It processes raw amplicon sequencing reads — from quality control and denoising through to taxonomy, diversity analysis, and differential abundance testing — making it suitable for microbial community profiling studies using either **16S rRNA** or **ITS** markers.

The pipeline handles the full analytical stack for amplicon-based metagenomics:

| Stage | What it does |
|---|---|
| QC | FastQC reports on raw reads |
| Import & demultiplex | Ingest reads via manifest or Casava convention |
| Primer trimming | Dedicated `q2-itsxpress` path for ITS; standard DADA2 trimming for 16S |
| Denoising | DADA2 for ASV inference, feature table construction, and representative sequence recovery |
| Phylogeny | MAFFT alignment + FastTree rooted tree (16S default; optional for ITS) |
| Diversity | Alpha and beta diversity metrics, rarefaction curves, and group significance tests |
| Taxonomy | Feature classification and interactive taxa barplots |
| Differential abundance | ANCOM across multiple taxonomic levels |
| Visualisation | Reproducible PDF plots for alpha diversity, PCoA, and taxonomic composition |

`metaQII-nf` converts the original monolithic `metaQII.py` script into a staged, reproducible Nextflow pipeline, eliminating hard-coded paths, interactive prompts, and desktop-only plotting logic.

---

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Inputs](#inputs)
- [Parameters](#parameters)
- [Outputs](#outputs)
- [Project Structure](#project-structure)
- [Notes](#notes)
- [References](#references)

---

## Requirements

| Dependency | Version |
|---|---|
| Nextflow | ≥ 23.04 |
| Conda / Mamba | any recent |
| QIIME2 amplicon | 2026.1 |
| OS | Linux or WSL2 |

> **Windows PowerShell is not supported.** QIIME2 and the bash-based Nextflow processes in this pipeline require a Unix-like environment.

---

## Installation

### 1. Install QIIME2 (recommended: use the official environment)

```bash
conda env create \
  --name qiime2-amplicon-2026.1 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml
```

Activate it and run the pipeline with `-profile standard` for exact parity with the upstream distribution.

### 2. Clone this repository

```bash
git clone https://github.com/<your-org>/metaqii-nf.git
cd metaqii-nf
```

### 3. (Alternative) Use the bundled environment

The pipeline ships `envs/qiime2-amplicon-2026.1-metaqii.yml`, which extends the official QIIME2 amplicon environment with `fastqc`, `itsxpress`, and `q2-itsxpress`. Use `-profile conda` to let Nextflow resolve this automatically.

---

## Quick Start

### Paired-end ITS

```bash
nextflow run main.nf -profile conda \
  --input_mode manifest \
  --manifest assets/example_manifest.tsv \
  --metadata assets/example_metadata.tsv \
  --classifier /path/to/classifier.qza \
  --paired_end true \
  --sequence_type ITS \
  --its_region ITS1 \
  --sampling_depth 10000 \
  --metadata_column state_of_disease \
  --trim_left_f 20 \
  --trim_left_r 20 \
  --trunc_len_f 250 \
  --trunc_len_r 220 \
  --dada2_n_reads_learn 1000000 \
  --outdir results_its
```

### Paired-end 16S

```bash
nextflow run main.nf -profile conda \
  --input_mode manifest \
  --manifest samples.tsv \
  --metadata metadata.tsv \
  --classifier /path/to/classifier.qza \
  --paired_end true \
  --sequence_type 16S \
  --sampling_depth 12000 \
  --metadata_column state_of_disease \
  --trim_left_f 17 \
  --trim_left_r 21 \
  --trunc_len_f 250 \
  --trunc_len_r 220 \
  --dada2_n_reads_learn 1000000 \
  --build_phylogeny true \
  --outdir results_16s
```

---

## Inputs

### Manifest mode (`--input_mode manifest`)

Provide a QIIME2-format manifest TSV. For paired-end data:

```tsv
sample-id	absolute-filepath	direction
SampleA	/data/SampleA_R1.fastq.gz	forward
SampleA	/data/SampleA_R2.fastq.gz	reverse
```

### Casava mode (`--input_mode casava`)

If reads already follow QIIME2 Casava naming conventions, point `--input_dir` at the folder directly. No manifest is needed.

Example files are provided in `assets/`:

- `assets/example_manifest.tsv`
- `assets/example_metadata.tsv`

---

## Parameters

Full parameter documentation is in [`docs/parameters.md`](docs/parameters.md).

### Input / output

| Parameter | Description | Default |
|---|---|---|
| `--input_mode` | `manifest` or `casava` | `manifest` |
| `--manifest` | Path to manifest TSV (manifest mode) | `null` |
| `--input_dir` | Path to Casava read directory (casava mode) | `null` |
| `--metadata` | Path to sample metadata TSV **[required]** | `null` |
| `--classifier` | Path to trained QIIME2 classifier `.qza` **[required]** | `null` |
| `--outdir` | Output directory | `results` |

### Sequencing layout

| Parameter | Description | Default |
|---|---|---|
| `--paired_end` | `true` for paired-end, `false` for single-end | `true` |
| `--sequence_type` | `ITS` or `16S` | `ITS` |
| `--its_region` | Target ITS region: `ITS1` or `ITS2` | `ITS1` |
| `--its_taxa` | ITSxpress taxa group (e.g. `F` for Fungi) | `F` |

### DADA2 denoising

| Parameter | Description | Default |
|---|---|---|
| `--trim_left_f` | Bases to trim from 5′ of forward reads | `0` |
| `--trim_left_r` | Bases to trim from 5′ of reverse reads | `0` |
| `--trunc_len_f` | Truncate forward reads at this position | `0` |
| `--trunc_len_r` | Truncate reverse reads at this position | `0` |
| `--dada2_pooling` | Pooling strategy: `independent`, `pseudo`, or `pooled` | `independent` |
| `--dada2_n_reads_learn` | Reads used for DADA2 error model learning | `1000000` |

### Diversity & taxonomy

| Parameter | Description | Default |
|---|---|---|
| `--sampling_depth` | Rarefaction depth for core diversity **[required]** | `null` |
| `--metadata_column` | Metadata column used for group comparisons **[required]** | `null` |
| `--build_phylogeny` | Build a rooted phylogenetic tree | `false` for ITS, `true` for 16S |
| `--taxonomy_start_level` | Taxonomic level to begin ANCOM collapse (1–7) | `2` |
| `--alpha_rarefaction_min_depth` | Minimum depth for alpha rarefaction curves | `1` |
| `--alpha_rarefaction_steps` | Number of steps in rarefaction curves | `10` |
| `--plot_top_n_taxa` | Top N taxa shown in composition plots | `15` |

### Execution

| Parameter | Description | Default |
|---|---|---|
| `--run_fastqc` | Run FastQC on raw reads before import | `true` |
| `--threads` | CPU threads per QIIME2 process | `4` |

---

## Outputs

All outputs are written under `--outdir` in the following structure:

```
results/
├── 00_qc/fastqc/            # Raw-read FastQC reports
├── 01_import/               # Imported demultiplexed QIIME2 artifact
├── 02_preprocessing/        # ITS-trimmed data, DADA2 feature table,
│                            # representative sequences, denoising stats
├── 03_summaries/            # Feature-table and denoising visualizations
├── 04_phylogeny/            # Alignment and rooted tree (when enabled)
├── 05_diversity/            # Core metrics, alpha rarefaction, significance tests
├── 06_taxonomy/             # Taxonomy assignments and taxa barplots
├── 07_ancom/                # Collapsed tables and ANCOM results per taxonomic level
├── 08_plots/                # PDF plots: alpha diversity, PCoA, taxonomic composition
└── pipeline_info/           # Nextflow trace, timeline, report, and DAG
```

---

## Project Structure

```
.
├── main.nf
├── nextflow.config
├── modules/local/metaqii.nf
├── envs/
│   ├── qiime2.yml                              # Legacy
│   ├── qiime2-amplicon-2026.1-upstream.yml
│   └── qiime2-amplicon-2026.1-metaqii.yml      # Active environment
├── bin/
│   ├── plot_alpha.py
│   ├── plot_pcoa.py
│   └── plot_taxa.py
├── assets/
│   ├── example_manifest.tsv
│   ├── example_metadata.tsv
│   └── NO_PHYLOGENY
└── docs/
    └── parameters.md
```

---

## Notes

**ITS vs 16S phylogeny behaviour** — `--build_phylogeny` defaults to `false` for ITS and `true` for 16S. ITS analyses are commonly run without phylogenetic diversity metrics. Pass `--build_phylogeny true` explicitly if you need them for ITS.

**`q2-itsxpress` availability** — the pipeline assumes this plugin is installed and importable when `--sequence_type ITS` is used. If you use the bundled environment (`envs/qiime2-amplicon-2026.1-metaqii.yml`), it is included. Installation details can vary by QIIME2 release; see the [q2-itsxpress repository](https://github.com/USDA-ARS-GBRU/q2_itsxpress) if you encounter import errors.

**Official QIIME2 environment vs bundled environment** — using the official install command produces exact parity with the upstream distribution. The bundled `metaqii.yml` adds only `fastqc`, `itsxpress`, and `q2-itsxpress` on top of that base.

---

## References

- [QIIME2 amplicon quickstart](https://library.qiime2.org/quickstart/amplicon)
- [QIIME2 diversity plugin reference](https://amplicon-docs.qiime2.org/en/latest/references/plugins/diversity.html)
- [QIIME2 news feed](https://qiime2.org/news)
- [q2-itsxpress usage](https://github.com/USDA-ARS-GBRU/q2_itsxpress)
- [scikit-bio ordination format](https://scikit.bio/docs/latest/generated/skbio.io.format.ordination.html)
