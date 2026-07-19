# Parameters

## Required

| Parameter | Description |
| --- | --- |
| `--metadata` | QIIME2-compatible sample metadata TSV. |
| `--classifier` | Pre-trained QIIME2 taxonomic classifier (`.qza`). |
| `--sampling_depth` | Rarefaction depth used for diversity workflows. |
| `--metadata_column` | Categorical metadata column used for beta significance and ANCOM. |
| `--manifest` or `--input_dir` | Input reads as a QIIME2 manifest or a Casava directory, depending on `--input_mode`. |

## Core workflow switches

| Parameter | Default | Notes |
| --- | --- | --- |
| `--input_mode` | `manifest` | Use `manifest` or `casava`. |
| `--paired_end` | `true` | Set `false` for single-end reads. |
| `--sequence_type` | `ITS` | Use `ITS` or `16S`. |
| `--run_fastqc` | `true` | Optional raw-read quality control. |
| `--build_phylogeny` | `null` | Auto-disables for ITS and auto-enables otherwise unless explicitly set. |
| `--its_region` | `ITS1` | `ITS1`, `ITS2`, or `ALL`. |
| `--its_taxa` | `F` | q2-itsxpress taxa code, `F` for fungi by default. |

## Primer removal (cutadapt)

Runs only when `--fwd_primer` is set (before ITS trimming / denoising). Typically used for 16S where primers remain on the reads.

| Parameter | Default | Notes |
| --- | --- | --- |
| `--fwd_primer` | `null` | Forward primer sequence; enables cutadapt when set. |
| `--rev_primer` | `null` | Reverse primer; **required** for paired-end. |
| `--cutadapt_error_rate` | `0.1` | Allowed primer mismatch rate. |
| `--cutadapt_discard_untrimmed` | `true` | Drop reads where the primer was not found. |

## DADA2 trimming

| Parameter | Default |
| --- | --- |
| `--trim_left_f` | `0` |
| `--trim_left_r` | `0` |
| `--trunc_len_f` | `0` |
| `--trunc_len_r` | `0` |
| `--dada2_pooling` | `independent` |
| `--dada2_n_reads_learn` | `1000000` |

Inspect `03_summaries/demux.qzv` (and `trimmed-demux.qzv` for ITS) to choose truncation positions.

## Feature filtering

| Parameter | Default | Notes |
| --- | --- | --- |
| `--filter_contaminants` | `null` | Auto: `true` for 16S, `false` for ITS. Removes contaminant lineages via `qiime taxa filter-table`. |
| `--contaminant_taxa` | `mitochondria,chloroplast` | Comma-separated substrings excluded when contaminant filtering is on. |
| `--filter_min_samples` | `1` | Drop features present in fewer than N samples (`1` = no filtering). |
| `--filter_min_frequency` | `0` | Drop features whose total frequency is below N (`0` = no filtering). |

## Control-based decontamination (decontam)

Runs only when `--decontam_control_column` is set. Uses the prevalence method against negative controls, then drops both contaminant features and the control samples themselves.

| Parameter | Default | Notes |
| --- | --- | --- |
| `--decontam_control_column` | `null` | Metadata column flagging control vs. true samples. |
| `--decontam_control_indicator` | `control` | Value in that column marking a control. |
| `--decontam_threshold` | `0.1` | decontam score cutoff; features with `p > threshold` are kept. |

## Taxonomy classification

| Parameter | Default | Notes |
| --- | --- | --- |
| `--classification_method` | `sklearn` | `sklearn` (trained classifier) or `vsearch` (reference-based). |
| `--classifier` | `null` | Trained classifier `.qza` (**required** for sklearn). |
| `--reference_reads` | `null` | `FeatureData[Sequence]` (**required** for vsearch). |
| `--reference_taxonomy` | `null` | `FeatureData[Taxonomy]` (**required** for vsearch). |
| `--classify_confidence` | `0.7` | sklearn confidence threshold. |
| `--classify_reads_per_batch` | `auto` | sklearn reads per batch. |
| `--vsearch_perc_identity` | `0.8` | vsearch minimum identity. |
| `--vsearch_maxaccepts` | `10` | vsearch max hits per query. |

sklearn classification runs with `--p-n-jobs`/`--p-threads` set to `--threads`. A trained sklearn classifier must match the QIIME 2 (scikit-learn) version — use `vsearch` if you hit version-mismatch errors.

## Phylogeny

| Parameter | Default | Notes |
| --- | --- | --- |
| `--phylogeny_method` | `denovo` | `denovo` (MAFFT + FastTree) or `sepp` (fragment insertion). |
| `--sepp_reference` | `null` | `SeppReferenceDatabase` `.qza` (**required** for sepp; e.g. `sepp-refs-gg-13-8.qza`). |

## Longitudinal (optional)

Runs `qiime longitudinal volatility` only when both columns are set.

| Parameter | Default | Notes |
| --- | --- | --- |
| `--longitudinal_state_column` | `null` | Time/state metadata column. |
| `--longitudinal_individual_column` | `null` | Subject/individual id column. |
| `--longitudinal_metric` | `shannon_entropy` | Default metric plotted. |

## Read QC aggregation

| Parameter | Default | Notes |
| --- | --- | --- |
| `--run_multiqc` | `true` | Aggregate FastQC with MultiQC (requires `--run_fastqc`; skipped if MultiQC is absent). |

## Differential abundance (ANCOM-BC2)

| Parameter | Default | Notes |
| --- | --- | --- |
| `--taxonomy_start_level` | `2` | First taxonomic level collapsed for the ANCOM-BC2 sweep (1–7). |
| `--da_significance_threshold` | `0.05` | q-value cutoff shown in the DA bar plots. |

## Reporting

| Parameter | Default |
| --- | --- |
| `--alpha_rarefaction_min_depth` | `1` |
| `--alpha_rarefaction_max_depth` | `null` (falls back to `--sampling_depth`; a warning is emitted) |
| `--alpha_rarefaction_steps` | `10` |
| `--plot_top_n_taxa` | `15` |
| `--threads` | `4` |
| `--outdir` | `results` |
