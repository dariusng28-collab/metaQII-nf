# metaQII-nf

> A reproducible Nextflow pipeline for QIIME2-based amplicon sequencing analysis (16S & ITS).

[![CI](https://github.com/dariusng28-collab/metaQII-nf/actions/workflows/ci.yml/badge.svg)](https://github.com/dariusng28-collab/metaQII-nf/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04-brightgreen)](https://www.nextflow.io/)
[![QIIME2](https://img.shields.io/badge/QIIME2-2026.1%20%7C%202026.4-blue)](https://qiime2.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20WSL2-lightgrey)]()

---

## Overview

`metaQII-nf` is an end-to-end amplicon metagenomics pipeline built on QIIME2. It processes raw amplicon sequencing reads — from quality control and denoising through to taxonomy, diversity analysis, and differential abundance testing — making it suitable for microbial community profiling studies using either **16S rRNA** or **ITS** markers.

The pipeline handles the full analytical stack for amplicon-based metagenomics:

| Stage | What it does |
|---|---|
| QC | FastQC reports on raw reads, aggregated with MultiQC |
| Import & demultiplex | Ingest reads via manifest or Casava convention |
| Read-quality summaries | `qiime demux summarize` interactive quality plots (pre- and post-trimming) to guide truncation |
| Primer trimming | Optional `q2-cutadapt` primer removal (16S); dedicated `q2-itsxpress` path for ITS |
| Denoising | DADA2 for ASV inference, feature table construction, and representative sequence recovery |
| Taxonomy | `classify-sklearn` **or** reference-based `classify-consensus-vsearch`, plus interactive taxa barplots |
| Feature filtering | Optional control-based `decontam`, host contaminant removal (mitochondria/chloroplast), and low-prevalence filtering |
| Phylogeny | MAFFT + FastTree de-novo tree, or SEPP fragment insertion against a reference (16S default; optional for ITS) |
| Diversity | Alpha and beta diversity metrics, rarefaction curves, group significance, and optional longitudinal volatility |
| Differential abundance | ANCOM-BC2 across multiple taxonomic levels, with bar-plot visualisations |
| Visualisation | Reproducible PDF plots (alpha diversity, PCoA, taxonomic composition) plus a consolidated HTML report |

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
| QIIME2 | 2026.1 (bundled/tested) or 2026.4+ |
| Container runtime *(optional, recommended)* | Docker (local) and/or Apptainer/Singularity (HPC) |
| OS | Linux or WSL2 |

> **QIIME 2 versions** — the bundled environment pins the **2026.1** amplicon distribution for reproducibility, and the pipeline is tested against it. It is also compatible with **2026.4+**, in which the `amplicon` distribution was renamed to `qiime2`. All commands used here (`demux summarize`, `taxa filter-table`, `feature-table filter-features`, `composition ancombc2` / `ancombc2-visualizer`) are available in both.

> **Windows PowerShell is not supported.** QIIME2 and the bash-based Nextflow processes in this pipeline require a Unix-like environment.

---

## Installation

### 1. Install QIIME2 (recommended: use the official environment)

Current release (2026.4+, distribution renamed to `qiime2`):

```bash
conda env create \
  --name qiime2-2026.4 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.4/qiime2/released/rachis-qiime2-linux-64-conda.yml
```

Or the pinned/tested release (2026.1):

```bash
conda env create \
  --name qiime2-amplicon-2026.1 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml
```

The ITS path additionally needs `fastqc`, `itsxpress`, and `q2-itsxpress` (add them to the environment, or use the bundled env below). Activate the environment and run with `-profile standard`.

### 2. Clone this repository

```bash
git clone https://github.com/dariusng28-collab/metaQII-nf.git
cd metaQII-nf
```

### 3. (Alternative) Use the bundled environment

The pipeline ships `envs/qiime2-amplicon-2026.1-metaqii.yml`, which extends the official QIIME2 amplicon environment with `fastqc`, `itsxpress`, and `q2-itsxpress`. Use `-profile conda` to let Nextflow resolve this automatically.

### 4. (Recommended) Containers

For the most reproducible runs, bake the pinned environment into an image once and reuse it everywhere. Build locally with Docker, then pull it with Apptainer/Singularity on your HPC.

```bash
# Build locally (uses the Dockerfile in the repo root)
docker build -t ghcr.io/dariusng28-collab/metaqii-nf:2026.1.0 .

# Option A — push to a registry so compute nodes can pull it
docker push ghcr.io/dariusng28-collab/metaqii-nf:2026.1.0
apptainer pull metaqii-nf_2026.1.0.sif docker://ghcr.io/dariusng28-collab/metaqii-nf:2026.1.0

# Option B — air-gapped: ship the image as a tarball
docker save ghcr.io/dariusng28-collab/metaqii-nf:2026.1.0 -o metaqii.tar
apptainer build metaqii-nf_2026.1.0.sif docker-archive://metaqii.tar
```

Point the pipeline at the image with `--container` (a registry ref or a local `.sif` path) and pick the runtime profile. Tagged releases (`v*`) also build and push the image automatically via GitHub Actions.

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

## Profiles

Profiles fall into two independent groups that you **compose** with a comma. Pick one from each:

| Group | Profiles | Purpose |
|---|---|---|
| Software | `docker`, `singularity`, `conda`, `prebuilt` | How QIIME 2 is provided |
| Executor | `standard` (local), `slurm`, `sge`, `lsf`, `pbs` | Where jobs run |

```bash
# Local development with Docker
nextflow run main.nf -profile docker ...

# Reproducible HPC run (Apptainer image + SLURM scheduler)
nextflow run main.nf -profile singularity,slurm --container /shared/images/metaqii-nf_2026.1.0.sif ...

# Cluster with a pre-activated conda environment
nextflow run main.nf -profile prebuilt,sge ...
```

Queue names can be overridden with `--slurm_queue` / `--sge_queue` / `--lsf_queue` / `--pbs_queue`.

There are also two check profiles: `-profile test` (ITS example, compile/DAG check via `-preview`) and `-profile smoke` (runnable 16S toy dataset — run `python bin/make_toy_16s_dataset.py` first to regenerate the manifest for your machine).

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

### Feature filtering

| Parameter | Description | Default |
|---|---|---|
| `--filter_contaminants` | Remove contaminant lineages (`taxa filter-table`) | `true` for 16S, `false` for ITS |
| `--contaminant_taxa` | Comma-separated substrings to exclude when filtering | `mitochondria,chloroplast` |
| `--filter_min_samples` | Drop features seen in fewer than N samples (`1` = off) | `1` |
| `--filter_min_frequency` | Drop features with total frequency below N (`0` = off) | `0` |

### Primer removal (cutadapt, optional)

| Parameter | Description | Default |
|---|---|---|
| `--fwd_primer` | Forward primer; enables cutadapt when set | `null` |
| `--rev_primer` | Reverse primer (required for paired-end) | `null` |
| `--cutadapt_error_rate` | Allowed primer mismatch rate | `0.1` |
| `--cutadapt_discard_untrimmed` | Drop reads without the primer | `true` |

### Control-based decontamination (optional)

| Parameter | Description | Default |
|---|---|---|
| `--decontam_control_column` | Metadata column flagging controls; enables decontam | `null` |
| `--decontam_control_indicator` | Value marking a control sample | `control` |
| `--decontam_threshold` | Keep features with decontam score `p >` this | `0.1` |

### Taxonomy classification

| Parameter | Description | Default |
|---|---|---|
| `--classification_method` | `sklearn` (trained classifier) or `vsearch` (reference-based) | `sklearn` |
| `--classifier` | Trained classifier `.qza` (required for sklearn) | `null` |
| `--reference_reads` / `--reference_taxonomy` | Reference `.qza`s (required for vsearch) | `null` |
| `--classify_confidence` | `classify-sklearn` confidence threshold | `0.7` |
| `--classify_reads_per_batch` | Reads per classification batch | `auto` |
| `--vsearch_perc_identity` / `--vsearch_maxaccepts` | vsearch identity / max hits | `0.8` / `10` |

### Phylogeny & longitudinal (optional)

| Parameter | Description | Default |
|---|---|---|
| `--phylogeny_method` | `denovo` (MAFFT+FastTree) or `sepp` (fragment insertion) | `denovo` |
| `--sepp_reference` | `SeppReferenceDatabase` `.qza` (required for sepp) | `null` |
| `--longitudinal_state_column` / `--longitudinal_individual_column` | Enable `longitudinal volatility` when both set | `null` |
| `--run_multiqc` | Aggregate FastQC with MultiQC | `true` |

### Diversity & differential abundance

| Parameter | Description | Default |
|---|---|---|
| `--sampling_depth` | Rarefaction depth for core diversity **[required]** | `null` |
| `--metadata_column` | Metadata column used for group comparisons **[required]** | `null` |
| `--build_phylogeny` | Build a rooted phylogenetic tree | `false` for ITS, `true` for 16S |
| `--taxonomy_start_level` | First taxonomic level for the ANCOM-BC2 sweep (1–7) | `2` |
| `--da_significance_threshold` | q-value cutoff shown in DA bar plots | `0.05` |
| `--alpha_rarefaction_min_depth` | Minimum depth for alpha rarefaction curves | `1` |
| `--alpha_rarefaction_max_depth` | Maximum depth for rarefaction curves | `null` → `--sampling_depth` (set higher; see notes) |
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
├── 00_qc/                     # Raw-read FastQC reports + aggregated MultiQC report
├── 01_import/                 # Imported demultiplexed QIIME2 artifact
├── 02_preprocessing/          # ITS-trimmed data, DADA2 feature table,
│   └── filtered/              #   representative sequences, and filtered table/seqs
├── 03_summaries/              # demux quality (demux.qzv), feature-table & denoising visualizations
├── 04_phylogeny/              # Alignment and rooted tree (when enabled)
├── 05_diversity/              # Core metrics, alpha rarefaction, significance tests
├── 06_taxonomy/               # Taxonomy assignments and taxa barplots
├── 07_differential_abundance/ # ANCOM-BC2 results & bar plots per taxonomic level
├── 08_plots/                  # PDF plots: alpha diversity, PCoA, taxonomic composition
├── 09_report/                 # Consolidated index.html linking every .qzv and figure
└── pipeline_info/             # trace, timeline, report, DAG + software_versions.yml,
                               # params.json and run_info.txt (provenance)
```

Open `09_report/index.html` in a browser for a single entry point; interactive `.qzv` files render at [view.qiime2.org](https://view.qiime2.org) (drag-and-drop, nothing is uploaded).

---

## Project Structure

```
.
├── main.nf
├── nextflow.config
├── modules/local/metaqii.nf
├── Dockerfile
├── .github/workflows/          # CI (compile/lint) + image build
├── LICENSE
├── envs/
│   ├── qiime2-amplicon-2026.1-upstream.yml
│   └── qiime2-amplicon-2026.1-metaqii.yml      # Active environment
├── bin/
│   ├── plot_alpha.py
│   ├── plot_pcoa.py
│   ├── plot_taxa.py
│   ├── build_report.py
│   └── make_toy_16s_dataset.py
├── assets/
│   ├── example_manifest.tsv
│   ├── example_metadata.tsv
│   └── NO_PHYLOGENY
└── docs/
    └── parameters.md
```

---

## Notes

**Choosing truncation lengths** — the pipeline now emits `03_summaries/demux.qzv` (and `trimmed-demux.qzv` for ITS). Open it at [view.qiime2.org](https://view.qiime2.org) to read per-base quality before setting `--trunc_len_f/r`.

**Alpha-rarefaction max depth** — `--alpha_rarefaction_max_depth` controls how far the rarefaction curve extends. When unset it falls back to `--sampling_depth` (and warns), which hides whether richness has plateaued. For a useful curve set it near the maximum per-sample frequency shown in `03_summaries/table.qzv`.

**Contaminant filtering** — for 16S, host mitochondria and chloroplast reads are removed before diversity and differential abundance (`--filter_contaminants`, on by default for 16S, off for ITS). Adjust the lineages with `--contaminant_taxa`, and optionally drop rare features with `--filter_min_samples` / `--filter_min_frequency`.

**Differential abundance** — differential abundance uses **ANCOM-BC2** (`qiime composition ancombc2` + `ancombc2-visualizer`), the method QIIME 2 now recommends over the legacy ANCOM. Results and bar plots are written per taxonomic level to `07_differential_abundance/`; a level is skipped (with a warning) rather than aborting the run if its collapsed table is too sparse to model.

**Primer removal** — 16S reads usually still carry primers; set `--fwd_primer`/`--rev_primer` to run `q2-cutadapt` before denoising. The ITS path relies on `q2-itsxpress` for primer/region trimming, so cutadapt is normally only needed for 16S.

**Control-based decontamination** — if your run includes negative controls, set `--decontam_control_column` (and `--decontam_control_indicator`) to remove reagent contaminants with `decontam` (prevalence method); the control samples are then dropped before diversity. This is more principled than the name-based mitochondria/chloroplast filter and can be combined with it.

**Classifier compatibility** — a trained `classify-sklearn` classifier must match the QIIME 2 / scikit-learn version it was built with. If you hit version-mismatch errors, switch to `--classification_method vsearch` with `--reference_reads`/`--reference_taxonomy`, which is not tied to the sklearn version.

**SEPP phylogeny** — `--phylogeny_method sepp` places ASVs into a reference tree (`--sepp_reference`, e.g. `sepp-refs-gg-13-8.qza` from the QIIME 2 data resources) instead of building a de-novo tree, which is more defensible for short reads.

**PICRUSt2 functional prediction is not wired in.** It requires the `q2-picrust2` plugin and a large reference that are not part of the bundled distribution/image. Running it would need its own environment/image, so it is intentionally left out rather than half-implemented; add it as a separate module if your study needs functional inference.

**ITS vs 16S phylogeny behaviour** — `--build_phylogeny` defaults to `false` for ITS and `true` for 16S. ITS analyses are commonly run without phylogenetic diversity metrics. Pass `--build_phylogeny true` explicitly if you need them for ITS.

**`q2-itsxpress` availability** — the pipeline assumes this plugin is installed and importable when `--sequence_type ITS` is used. If you use the bundled environment (`envs/qiime2-amplicon-2026.1-metaqii.yml`), it is included. Installation details can vary by QIIME2 release; see the [q2-itsxpress repository](https://github.com/USDA-ARS-GBRU/q2_itsxpress) if you encounter import errors.

**Official QIIME2 environment vs bundled environment** — using the official install command produces exact parity with the upstream distribution. The bundled `metaqii.yml` adds only `fastqc`, `itsxpress`, and `q2-itsxpress` on top of that base.

---

## References

- [QIIME2 amplicon documentation (stable)](https://amplicon-docs.qiime2.org/en/stable/)
- [QIIME2 amplicon quickstart](https://library.qiime2.org/quickstart/amplicon)
- [QIIME2 diversity plugin reference](https://amplicon-docs.qiime2.org/en/latest/references/plugins/diversity.html)
- [q2-composition (ANCOM-BC2) reference](https://library.qiime2.org/plugins/qiime2/q2-composition/overview)
- [QIIME2 news feed](https://qiime2.org/news)
- [q2-itsxpress usage](https://github.com/USDA-ARS-GBRU/q2_itsxpress)
- [scikit-bio ordination format](https://scikit.bio/docs/latest/generated/skbio.io.format.ordination.html)
