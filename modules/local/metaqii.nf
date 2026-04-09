process RUN_FASTQC {
    tag { input_mode }
    label 'qiime'
    publishDir "${params.outdir}/00_qc/fastqc", mode: 'copy'

    input:
    tuple val(input_mode), path(input_source)

    output:
    path 'fastqc'

    script:
    """
    mkdir -p fastqc

    if [[ "${input_mode}" == "manifest" ]]; then
        python - <<'PY'
import csv
import subprocess

files = []
with open("${input_source}", newline="") as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        fastq = row.get("absolute-filepath")
        if fastq:
            files.append(fastq)

files = sorted(set(files))
if not files:
    raise SystemExit("No FASTQ files found in manifest ${input_source}")

subprocess.run(
    ["fastqc", "--threads", "${task.cpus}", "--outdir", "fastqc", *files],
    check=True,
)
PY
    else
        shopt -s nullglob
        fastqs=( "${input_source}"/*.fastq.gz "${input_source}"/*.fq.gz )
        if [[ \${#fastqs[@]} -eq 0 ]]; then
            echo "No FASTQ files found in ${input_source}" >&2
            exit 1
        fi
        fastqc --threads ${task.cpus} --outdir fastqc "\${fastqs[@]}"
    fi
    """
}

process IMPORT_QIIME {
    tag { input_mode }
    label 'qiime'
    publishDir "${params.outdir}/01_import", mode: 'copy', pattern: '*.qza'

    input:
    tuple val(input_mode), path(input_source)

    output:
    path 'demux.qza', emit: demux

    script:
    def qiimeType = params.paired_end ? 'SampleData[PairedEndSequencesWithQuality]' : 'SampleData[SequencesWithQuality]'
    def qiimeFormat = input_mode == 'manifest'
        ? (params.paired_end ? 'PairedEndFastqManifestPhred33V2' : 'SingleEndFastqManifestPhred33V2')
        : 'CasavaOneEightSingleLanePerSampleDirFmt'

    """
    qiime tools import \\
      --type '${qiimeType}' \\
      --input-path "${input_source}" \\
      --input-format '${qiimeFormat}' \\
      --output-path demux.qza
    """
}

process ITSXPRESS_TRIM {
    label 'qiime'
    publishDir "${params.outdir}/02_preprocessing", mode: 'copy', pattern: '*.qza'

    input:
    path demux_qza

    output:
    path 'trimmed_demux.qza', emit: trimmed_demux

    script:
    def commandName = params.paired_end ? 'trim-pair-output-unmerged' : 'trim-single'

    """
    qiime itsxpress ${commandName} \\
      --i-per-sample-sequences "${demux_qza}" \\
      --p-region ${params.its_region} \\
      --p-taxa ${params.its_taxa} \\
      --p-threads ${task.cpus} \\
      --o-trimmed trimmed_demux.qza
    """
}

process DADA2_DENOISE {
    label 'qiime'
    publishDir "${params.outdir}/02_preprocessing", mode: 'copy', pattern: '*.qza'

    input:
    path demux_qza

    output:
    path 'table.qza', emit: table
    path 'rep-seqs.qza', emit: rep_seqs
    path 'denoising-stats.qza', emit: denoising_stats
    path 'base-transition-stats.qza', emit: base_transition_stats

    script:
    def pooling = (params.dada2_pooling ?: 'independent').toString()
    def nReadsLearn = (params.dada2_n_reads_learn ?: 1000000) as int
    def trimLeftF = (params.trim_left_f ?: 0) as int
    def trimLeftR = (params.trim_left_r ?: 0) as int
    def truncLenF = (params.trunc_len_f ?: 0) as int
    def truncLenR = (params.trunc_len_r ?: 0) as int

    if (params.paired_end) {
        """
        qiime dada2 denoise-paired \\
          --i-demultiplexed-seqs "${demux_qza}" \\
          --p-trim-left-f ${trimLeftF} \\
          --p-trim-left-r ${trimLeftR} \\
          --p-trunc-len-f ${truncLenF} \\
          --p-trunc-len-r ${truncLenR} \\
          --p-pooling-method ${pooling} \\
          --p-n-threads ${task.cpus} \\
          --p-n-reads-learn ${nReadsLearn} \\
          --o-table table.qza \\
          --o-representative-sequences rep-seqs.qza \\
          --o-denoising-stats denoising-stats.qza \\
          --o-base-transition-stats base-transition-stats.qza
        """
    } else {
        """
        qiime dada2 denoise-single \\
          --i-demultiplexed-seqs "${demux_qza}" \\
          --p-trim-left ${trimLeftF} \\
          --p-trunc-len ${truncLenF} \\
          --p-pooling-method ${pooling} \\
          --p-n-threads ${task.cpus} \\
          --p-n-reads-learn ${nReadsLearn} \\
          --o-table table.qza \\
          --o-representative-sequences rep-seqs.qza \\
          --o-denoising-stats denoising-stats.qza \\
          --o-base-transition-stats base-transition-stats.qza
        """
    }
}

process SUMMARIZE_FEATURES {
    label 'qiime'
    publishDir "${params.outdir}/03_summaries", mode: 'copy', pattern: '*.qzv'

    input:
    path table_qza
    path rep_seqs_qza
    path denoising_stats_qza
    path metadata_tsv

    output:
    path 'denoising-stats.qzv'
    path 'table.qzv'
    path 'rep-seqs.qzv'

    script:
    """
    qiime metadata tabulate \\
      --m-input-file "${denoising_stats_qza}" \\
      --o-visualization denoising-stats.qzv

    qiime feature-table summarize \\
      --i-table "${table_qza}" \\
      --m-sample-metadata-file "${metadata_tsv}" \\
      --o-visualization table.qzv

    qiime feature-table tabulate-seqs \\
      --i-data "${rep_seqs_qza}" \\
      --o-visualization rep-seqs.qzv
    """
}

process BUILD_PHYLOGENY {
    label 'qiime'
    publishDir "${params.outdir}/04_phylogeny", mode: 'copy'

    input:
    path rep_seqs_qza

    output:
    path 'phylogeny/rooted_tree.qza', emit: rooted_tree
    path 'phylogeny/unrooted_tree.qza'
    path 'phylogeny/aligned_rep_seqs.qza'
    path 'phylogeny/masked_aligned_rep_seqs.qza'

    script:
    """
    qiime phylogeny align-to-tree-mafft-fasttree \\
      --i-sequences "${rep_seqs_qza}" \\
      --p-n-threads ${task.cpus} \\
      --output-dir phylogeny
    """
}

process DIVERSITY_ANALYSIS {
    label 'qiime'
    publishDir "${params.outdir}/05_diversity", mode: 'copy'

    input:
    path table_qza
    path metadata_tsv
    path rooted_tree
    val build_phylogeny

    output:
    path 'core_metrics', emit: core_metrics_dir

    script:
    def samplingDepth = params.sampling_depth as int

    """
    mkdir -p core_metrics

    if [[ "${build_phylogeny}" == "true" ]]; then
        qiime diversity core-metrics-phylogenetic \\
          --i-phylogeny "${rooted_tree}" \\
          --i-table "${table_qza}" \\
          --p-sampling-depth ${samplingDepth} \\
          --m-metadata-file "${metadata_tsv}" \\
          --output-dir core_metrics
    else
        qiime diversity core-metrics \\
          --i-table "${table_qza}" \\
          --p-sampling-depth ${samplingDepth} \\
          --m-metadata-file "${metadata_tsv}" \\
          --output-dir core_metrics
    fi
    """
}

process ALPHA_RAREFACTION {
    label 'qiime'
    publishDir "${params.outdir}/05_diversity", mode: 'copy', pattern: 'alpha-rarefaction.qzv'

    input:
    path table_qza
    path metadata_tsv
    path rooted_tree
    val build_phylogeny

    output:
    path 'alpha-rarefaction.qzv'

    script:
    def minDepth = (params.alpha_rarefaction_min_depth ?: 1) as int
    def maxDepth = params.sampling_depth as int
    def steps = (params.alpha_rarefaction_steps ?: 10) as int

    """
    if [[ "${build_phylogeny}" == "true" ]]; then
        qiime diversity alpha-rarefaction \\
          --i-table "${table_qza}" \\
          --i-phylogeny "${rooted_tree}" \\
          --p-min-depth ${minDepth} \\
          --p-max-depth ${maxDepth} \\
          --p-steps ${steps} \\
          --m-metadata-file "${metadata_tsv}" \\
          --o-visualization alpha-rarefaction.qzv
    else
        qiime diversity alpha-rarefaction \\
          --i-table "${table_qza}" \\
          --p-min-depth ${minDepth} \\
          --p-max-depth ${maxDepth} \\
          --p-steps ${steps} \\
          --m-metadata-file "${metadata_tsv}" \\
          --o-visualization alpha-rarefaction.qzv
    fi
    """
}

process ALPHA_GROUP_SIGNIFICANCE {
    label 'qiime'
    publishDir "${params.outdir}/05_diversity/alpha_significance", mode: 'copy'

    input:
    path core_metrics_dir
    path metadata_tsv

    output:
    path 'alpha_significance'

    script:
    """
    mkdir -p alpha_significance

    for metric in shannon_vector observed_features_vector evenness_vector faith_pd_vector; do
        if [[ -f "${core_metrics_dir}/\${metric}.qza" ]]; then
            base="\${metric%_vector}"
            qiime diversity alpha-group-significance \\
              --i-alpha-diversity "${core_metrics_dir}/\${metric}.qza" \\
              --m-metadata-file "${metadata_tsv}" \\
              --o-visualization "alpha_significance/\${base}-group-significance.qzv"
        fi
    done
    """
}

process BETA_GROUP_SIGNIFICANCE {
    label 'qiime'
    publishDir "${params.outdir}/05_diversity/beta_significance", mode: 'copy'

    input:
    path core_metrics_dir
    path metadata_tsv

    output:
    path 'beta_significance'

    script:
    """
    mkdir -p beta_significance

    for metric in bray_curtis jaccard weighted_unifrac unweighted_unifrac; do
        if [[ -f "${core_metrics_dir}/\${metric}_distance_matrix.qza" ]]; then
            qiime diversity beta-group-significance \\
              --i-distance-matrix "${core_metrics_dir}/\${metric}_distance_matrix.qza" \\
              --m-metadata-file "${metadata_tsv}" \\
              --m-metadata-column "${params.metadata_column}" \\
              --p-pairwise \\
              --o-visualization "beta_significance/\${metric}-group-significance.qzv"
        fi
    done
    """
}

process CLASSIFY_TAXONOMY {
    label 'qiime'
    publishDir "${params.outdir}/06_taxonomy", mode: 'copy'

    input:
    path rep_seqs_qza
    path table_qza
    path metadata_tsv
    path classifier_qza

    output:
    path 'taxonomy.qza', emit: taxonomy_qza
    path 'taxonomy.qzv'
    path 'taxa-bar-plots.qzv'

    script:
    """
    qiime feature-classifier classify-sklearn \\
      --i-classifier "${classifier_qza}" \\
      --i-reads "${rep_seqs_qza}" \\
      --o-classification taxonomy.qza

    qiime metadata tabulate \\
      --m-input-file taxonomy.qza \\
      --o-visualization taxonomy.qzv

    qiime taxa barplot \\
      --i-table "${table_qza}" \\
      --i-taxonomy taxonomy.qza \\
      --m-metadata-file "${metadata_tsv}" \\
      --o-visualization taxa-bar-plots.qzv
    """
}

process ANCOM_SWEEP {
    label 'qiime'
    publishDir "${params.outdir}/07_ancom", mode: 'copy'

    input:
    path table_qza
    path taxonomy_qza
    path metadata_tsv

    output:
    path 'ancom'

    script:
    def startLevel = (params.taxonomy_start_level ?: 2) as int

    """
    mkdir -p ancom

    for level in \$(seq ${startLevel} 7); do
        qiime taxa collapse \\
          --i-table "${table_qza}" \\
          --i-taxonomy "${taxonomy_qza}" \\
          --p-level "\${level}" \\
          --o-collapsed-table "ancom/collapsed-taxonomy-table-level\${level}.qza"

        qiime composition add-pseudocount \\
          --i-table "ancom/collapsed-taxonomy-table-level\${level}.qza" \\
          --o-composition-table "ancom/comp-ancom-table-level\${level}.qza"

        qiime composition ancom \\
          --i-table "ancom/comp-ancom-table-level\${level}.qza" \\
          --m-metadata-file "${metadata_tsv}" \\
          --m-metadata-column "${params.metadata_column}" \\
          --o-visualization "ancom/ancom-results-level\${level}.qzv"
    done
    """
}

process PLOT_REPORTS {
    label 'qiime'
    publishDir "${params.outdir}/08_plots", mode: 'copy'

    input:
    path core_metrics_dir
    path table_qza
    path taxonomy_qza
    path metadata_tsv

    output:
    path 'plots'

    script:
    """
    mkdir -p plots/exports/alpha plots/exports/pcoa plots/exports/table plots/exports/taxonomy

    for metric in shannon_vector observed_features_vector evenness_vector faith_pd_vector; do
        if [[ -f "${core_metrics_dir}/\${metric}.qza" ]]; then
            qiime tools export \\
              --input-path "${core_metrics_dir}/\${metric}.qza" \\
              --output-path "plots/exports/alpha/\${metric}"
        fi
    done

    for metric in bray_curtis_pcoa_results jaccard_pcoa_results weighted_unifrac_pcoa_results unweighted_unifrac_pcoa_results; do
        if [[ -f "${core_metrics_dir}/\${metric}.qza" ]]; then
            qiime tools export \\
              --input-path "${core_metrics_dir}/\${metric}.qza" \\
              --output-path "plots/exports/pcoa/\${metric}"
        fi
    done

    qiime tools export \\
      --input-path "${table_qza}" \\
      --output-path plots/exports/table

    qiime tools export \\
      --input-path "${taxonomy_qza}" \\
      --output-path plots/exports/taxonomy

    python "${projectDir}/bin/plot_alpha.py" \\
      --metadata "${metadata_tsv}" \\
      --metadata-column "${params.metadata_column}" \\
      --alpha-dir plots/exports/alpha \\
      --outdir plots

    python "${projectDir}/bin/plot_pcoa.py" \\
      --metadata "${metadata_tsv}" \\
      --metadata-column "${params.metadata_column}" \\
      --pcoa-dir plots/exports/pcoa \\
      --outdir plots

    python "${projectDir}/bin/plot_taxa.py" \\
      --metadata "${metadata_tsv}" \\
      --metadata-column "${params.metadata_column}" \\
      --biom plots/exports/table/feature-table.biom \\
      --taxonomy plots/exports/taxonomy/taxonomy.tsv \\
      --outdir plots \\
      --top-n ${params.plot_top_n_taxa}
    """
}
