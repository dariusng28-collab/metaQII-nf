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

## DADA2 trimming

| Parameter | Default |
| --- | --- |
| `--trim_left_f` | `0` |
| `--trim_left_r` | `0` |
| `--trunc_len_f` | `0` |
| `--trunc_len_r` | `0` |
| `--dada2_pooling` | `independent` |
| `--dada2_n_reads_learn` | `1000000` |

## Reporting

| Parameter | Default |
| --- | --- |
| `--taxonomy_start_level` | `2` |
| `--alpha_rarefaction_min_depth` | `1` |
| `--alpha_rarefaction_steps` | `10` |
| `--plot_top_n_taxa` | `15` |
| `--threads` | `4` |
| `--outdir` | `results` |
