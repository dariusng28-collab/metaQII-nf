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

# QIIME 2 manifests may be comma- or tab-separated. Sniff the delimiter
# instead of assuming CSV so tab-separated manifests (the documented
# default) are parsed correctly.
path = "${input_source}"
with open(path, newline="") as handle:
    sample = handle.read(4096)
    handle.seek(0)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\\t")
    except csv.Error:
        dialect = csv.excel_tab if "\\t" in sample.splitlines()[0] else csv.excel
    reader = csv.DictReader(handle, dialect=dialect)
    files = []
    for row in reader:
        fastq = row.get("absolute-filepath")
        if fastq:
            files.append(fastq.strip())

files = sorted(set(files))
if not files:
    raise SystemExit(f"No FASTQ files found in manifest {path}")

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

process MULTIQC {
    label 'qiime'
    publishDir "${params.outdir}/00_qc/multiqc", mode: 'copy'

    input:
    path fastqc_dir

    output:
    path 'multiqc_report.html', emit: report, optional: true
    path 'multiqc_data', optional: true

    script:
    // multiqc is an optional extra; guard so the run does not fail when it is
    // absent (e.g. -profile standard against an environment without it).
    """
    if command -v multiqc >/dev/null 2>&1; then
        multiqc "${fastqc_dir}" --filename multiqc_report.html
    else
        echo "[metaQII] multiqc not found in environment; skipping QC aggregation." >&2
    fi
    """
}

process IMPORT_QIIME {
    tag { input_mode }
    label 'qiime'
    publishDir "${params.outdir}/01_import", mode: 'copy', pattern: '*.qza'

    input:
    tuple val(input_mode), path(input_source)
    val paired_end

    output:
    path 'demux.qza', emit: demux

    script:
    def qiimeType = paired_end ? 'SampleData[PairedEndSequencesWithQuality]' : 'SampleData[SequencesWithQuality]'
    def qiimeFormat = input_mode == 'manifest'
        ? (paired_end ? 'PairedEndFastqManifestPhred33V2' : 'SingleEndFastqManifestPhred33V2')
        : 'CasavaOneEightSingleLanePerSampleDirFmt'

    """
    qiime tools import \\
      --type '${qiimeType}' \\
      --input-path "${input_source}" \\
      --input-format '${qiimeFormat}' \\
      --output-path demux.qza
    """
}

process DEMUX_SUMMARIZE {
    tag { label_name }
    label 'qiime'
    publishDir "${params.outdir}/03_summaries", mode: 'copy', pattern: '*.qzv'

    input:
    tuple val(label_name), path(demux_qza)

    output:
    path "${label_name}.qzv", emit: qzv

    script:
    """
    qiime demux summarize \\
      --i-data "${demux_qza}" \\
      --o-visualization "${label_name}.qzv"
    """
}

process CUTADAPT_TRIM {
    label 'qiime'
    publishDir "${params.outdir}/02_preprocessing", mode: 'copy', pattern: '*.qza'

    input:
    path demux_qza
    val paired_end

    output:
    path 'primer-trimmed.qza', emit: trimmed

    script:
    def errRate = params.cutadapt_error_rate ?: 0.1
    def discard = ((params.cutadapt_discard_untrimmed == null ? true : params.cutadapt_discard_untrimmed)
                    ? '--p-discard-untrimmed' : '--p-no-discard-untrimmed')
    def fwd = params.fwd_primer
    def rev = params.rev_primer

    if (paired_end) {
        """
        qiime cutadapt trim-paired \\
          --i-demultiplexed-sequences "${demux_qza}" \\
          --p-front-f ${fwd} \\
          --p-front-r ${rev} \\
          --p-error-rate ${errRate} \\
          ${discard} \\
          --p-cores ${task.cpus} \\
          --o-trimmed-sequences primer-trimmed.qza
        """
    } else {
        """
        qiime cutadapt trim-single \\
          --i-demultiplexed-sequences "${demux_qza}" \\
          --p-front ${fwd} \\
          --p-error-rate ${errRate} \\
          ${discard} \\
          --p-cores ${task.cpus} \\
          --o-trimmed-sequences primer-trimmed.qza
        """
    }
}

process ITSXPRESS_TRIM {
    label 'qiime'
    publishDir "${params.outdir}/02_preprocessing", mode: 'copy', pattern: '*.qza'

    input:
    path demux_qza
    val paired_end

    output:
    path 'trimmed_demux.qza', emit: trimmed_demux

    script:
    def commandName = paired_end ? 'trim-pair-output-unmerged' : 'trim-single'

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
    label 'qiime_heavy'
    publishDir "${params.outdir}/02_preprocessing", mode: 'copy', pattern: '*.qza'

    input:
    path demux_qza
    val paired_end

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

    if (paired_end) {
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

process CLASSIFY_TAXONOMY {
    label 'qiime_heavy'
    publishDir "${params.outdir}/06_taxonomy", mode: 'copy', pattern: '*.qz*'

    input:
    path rep_seqs_qza
    path classifier_qza
    path ref_reads_qza
    path ref_taxonomy_qza

    output:
    path 'taxonomy.qza', emit: taxonomy_qza
    path 'taxonomy.qzv', emit: taxonomy_qzv
    path 'vsearch-hits.qza', optional: true

    script:
    def method = (params.classification_method ?: 'sklearn').toString().toLowerCase()
    def confidence = (params.classify_confidence ?: 0.7)
    def readsPerBatch = (params.classify_reads_per_batch ?: 'auto').toString()
    def percId = (params.vsearch_perc_identity ?: 0.8)
    def maxAccepts = (params.vsearch_maxaccepts ?: 10)

    if (method == 'vsearch') {
        """
        qiime feature-classifier classify-consensus-vsearch \\
          --i-query "${rep_seqs_qza}" \\
          --i-reference-reads "${ref_reads_qza}" \\
          --i-reference-taxonomy "${ref_taxonomy_qza}" \\
          --p-perc-identity ${percId} \\
          --p-maxaccepts ${maxAccepts} \\
          --p-threads ${task.cpus} \\
          --o-classification taxonomy.qza \\
          --o-search-results vsearch-hits.qza

        qiime metadata tabulate \\
          --m-input-file taxonomy.qza \\
          --o-visualization taxonomy.qzv
        """
    } else {
        """
        qiime feature-classifier classify-sklearn \\
          --i-classifier "${classifier_qza}" \\
          --i-reads "${rep_seqs_qza}" \\
          --p-confidence ${confidence} \\
          --p-reads-per-batch ${readsPerBatch} \\
          --p-n-jobs ${task.cpus} \\
          --o-classification taxonomy.qza

        qiime metadata tabulate \\
          --m-input-file taxonomy.qza \\
          --o-visualization taxonomy.qzv
        """
    }
}

process FILTER_TABLE {
    label 'qiime'
    publishDir "${params.outdir}/02_preprocessing/filtered", mode: 'copy', pattern: '*.qza'

    input:
    path table_qza
    path rep_seqs_qza
    path taxonomy_qza
    path metadata_tsv
    val filter_contaminants

    output:
    path 'filtered-table.qza', emit: table
    path 'filtered-rep-seqs.qza', emit: rep_seqs
    path 'decontam-scores.qza', optional: true

    script:
    def contaminantTaxa = (params.contaminant_taxa ?: 'mitochondria,chloroplast').toString()
    def minSamples = (params.filter_min_samples ?: 1) as int
    def minFrequency = (params.filter_min_frequency ?: 0) as int
    def controlCol = params.decontam_control_column
    def controlInd = (params.decontam_control_indicator ?: 'control').toString()
    def decontamThreshold = (params.decontam_threshold ?: 0.1)

    """
    TABLE="${table_qza}"

    # 0. Control-based decontamination (prevalence method) — optional.
    if [[ -n "${controlCol ?: ''}" ]]; then
        qiime quality-control decontam-identify \\
          --i-table "\$TABLE" \\
          --m-metadata-file "${metadata_tsv}" \\
          --p-method prevalence \\
          --p-prev-control-column "${controlCol}" \\
          --p-prev-control-indicator "${controlInd}" \\
          --o-decontam-scores decontam-scores.qza

        # Keep non-contaminant features (decontam score p > threshold, or unscored).
        qiime feature-table filter-features \\
          --i-table "\$TABLE" \\
          --m-metadata-file decontam-scores.qza \\
          --p-where '[p]>${decontamThreshold} OR [p] IS NULL' \\
          --o-filtered-table decontam-table.qza

        # Drop the control samples themselves before downstream analysis.
        qiime feature-table filter-samples \\
          --i-table decontam-table.qza \\
          --m-metadata-file "${metadata_tsv}" \\
          --p-where "[${controlCol}]!='${controlInd}'" \\
          --o-filtered-table decontam-samples-table.qza
        TABLE=decontam-samples-table.qza
    fi

    # 1. Remove contaminant lineages (e.g. host mitochondria / chloroplast).
    if [[ "${filter_contaminants}" == "true" ]]; then
        qiime taxa filter-table \\
          --i-table "\$TABLE" \\
          --i-taxonomy "${taxonomy_qza}" \\
          --p-exclude "${contaminantTaxa}" \\
          --p-mode contains \\
          --o-filtered-table taxa-filtered-table.qza
    else
        cp "\$TABLE" taxa-filtered-table.qza
    fi

    # 2. Drop low-prevalence / low-abundance features.
    qiime feature-table filter-features \\
      --i-table taxa-filtered-table.qza \\
      --p-min-samples ${minSamples} \\
      --p-min-frequency ${minFrequency} \\
      --o-filtered-table filtered-table.qza

    # 3. Keep representative sequences in sync with the filtered table.
    qiime feature-table filter-seqs \\
      --i-data "${rep_seqs_qza}" \\
      --i-table filtered-table.qza \\
      --o-filtered-data filtered-rep-seqs.qza
    """
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
    path 'denoising-stats.qzv', emit: denoising_stats_qzv
    path 'table.qzv', emit: table_qzv
    path 'rep-seqs.qzv', emit: rep_seqs_qzv

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

process TAXA_BARPLOT {
    label 'qiime'
    publishDir "${params.outdir}/06_taxonomy", mode: 'copy', pattern: '*.qzv'

    input:
    path table_qza
    path taxonomy_qza
    path metadata_tsv

    output:
    path 'taxa-bar-plots.qzv', emit: qzv

    script:
    """
    qiime taxa barplot \\
      --i-table "${table_qza}" \\
      --i-taxonomy "${taxonomy_qza}" \\
      --m-metadata-file "${metadata_tsv}" \\
      --o-visualization taxa-bar-plots.qzv
    """
}

process BUILD_PHYLOGENY {
    label 'qiime_heavy'
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

process BUILD_PHYLOGENY_SEPP {
    label 'qiime_heavy'
    publishDir "${params.outdir}/04_phylogeny", mode: 'copy'

    input:
    path rep_seqs_qza
    path sepp_reference

    output:
    path 'sepp_tree.qza', emit: rooted_tree
    path 'sepp_placements.qza'

    script:
    """
    qiime fragment-insertion sepp \\
      --i-representative-sequences "${rep_seqs_qza}" \\
      --i-reference-database "${sepp_reference}" \\
      --p-threads ${task.cpus} \\
      --o-tree sepp_tree.qza \\
      --o-placements sepp_placements.qza
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
    val max_depth

    output:
    path 'alpha-rarefaction.qzv', emit: qzv

    script:
    def minDepth = (params.alpha_rarefaction_min_depth ?: 1) as int
    def maxDepth = max_depth as int
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

process LONGITUDINAL {
    label 'qiime'
    publishDir "${params.outdir}/05_diversity/longitudinal", mode: 'copy'

    input:
    path core_metrics_dir
    path metadata_tsv

    output:
    path 'volatility.qzv', optional: true

    script:
    def stateCol = params.longitudinal_state_column
    def indivCol = params.longitudinal_individual_column
    def metric = (params.longitudinal_metric ?: 'shannon_entropy').toString()

    """
    # Feed the Shannon vector as extra metadata so it is plottable, then guard
    # the call so an unsuitable design warns instead of aborting the run.
    EXTRA=""
    if [[ -f "${core_metrics_dir}/shannon_vector.qza" ]]; then
        EXTRA="--m-metadata-file ${core_metrics_dir}/shannon_vector.qza"
    fi

    set +e
    qiime longitudinal volatility \\
      --m-metadata-file "${metadata_tsv}" \\
      \$EXTRA \\
      --p-state-column "${stateCol}" \\
      --p-individual-id-column "${indivCol}" \\
      --p-default-metric "${metric}" \\
      --p-default-group-column "${params.metadata_column}" \\
      --o-visualization volatility.qzv
    status=\$?
    set -e
    if [[ \$status -ne 0 ]]; then
        echo "[metaQII] longitudinal volatility skipped (exit \${status}); check state/individual columns." >&2
    fi
    """
}

process ANCOMBC2_SWEEP {
    label 'qiime'
    publishDir "${params.outdir}/07_differential_abundance", mode: 'copy'

    input:
    path table_qza
    path taxonomy_qza
    path metadata_tsv

    output:
    path 'ancombc2'

    script:
    def startLevel = (params.taxonomy_start_level ?: 2) as int
    def column = params.metadata_column
    def sigThreshold = (params.da_significance_threshold ?: 0.05)

    """
    mkdir -p ancombc2

    for level in \$(seq ${startLevel} 7); do
        # ANCOM-BC2 needs at least two taxa to model; guard each level so a
        # sparse collapse at one rank does not abort the whole sweep.
        set +e
        qiime taxa collapse \\
          --i-table "${table_qza}" \\
          --i-taxonomy "${taxonomy_qza}" \\
          --p-level "\${level}" \\
          --o-collapsed-table "ancombc2/collapsed-table-level\${level}.qza" \\
        && qiime composition ancombc2 \\
          --i-table "ancombc2/collapsed-table-level\${level}.qza" \\
          --m-metadata-file "${metadata_tsv}" \\
          --p-fixed-effects-formula "${column}" \\
          --o-ancombc2-output "ancombc2/ancombc2-level\${level}.qza" \\
        && qiime composition ancombc2-visualizer \\
          --i-data "ancombc2/ancombc2-level\${level}.qza" \\
          --p-significance-threshold ${sigThreshold} \\
          --o-visualization "ancombc2/ancombc2-level\${level}.qzv"
        status=\$?
        set -e
        if [[ \$status -ne 0 ]]; then
            echo "[metaQII] ANCOM-BC2 skipped at level \${level} (exit \${status}); continuing." >&2
            rm -f "ancombc2/ancombc2-level\${level}.qza"
        fi
    done

    if ! ls ancombc2/*.qzv >/dev/null 2>&1; then
        echo "[metaQII] WARNING: ANCOM-BC2 produced no results for any taxonomic level." >&2
    fi
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
    path 'plots', emit: plots

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

process BUILD_REPORT {
    label 'qiime'
    publishDir "${params.outdir}/09_report", mode: 'copy'

    input:
    path qzv_files
    path plots_dir

    output:
    path 'index.html'
    path '*.qzv'
    path "${plots_dir}"

    script:
    // The collected .qzv files and the plots directory are staged into the
    // task directory under their own names, so build_report.py scans them in
    // place; publishing re-exports them into a self-contained report folder.
    """
    python "${projectDir}/bin/build_report.py" \\
      --outdir . \\
      --title "metaQII-nf report"
    """
}

process DUMP_VERSIONS {
    label 'qiime'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    output:
    path 'software_versions.yml'

    script:
    // Each tool is guarded so a missing binary records "n/a" rather than
    // failing the run; this file makes results traceable to exact versions.
    """
    {
      echo "metaQII-nf: '${workflow.manifest.version ?: 'unknown'}'"
      echo "nextflow: '${workflow.nextflow.version}'"
      echo "qiime2: '\$(qiime --version 2>/dev/null | head -n1 | sed 's/^q2cli version //;s/\\r//' || echo n/a)'"
      echo "fastqc: '\$(fastqc --version 2>/dev/null | sed 's/^FastQC //' || echo n/a)'"
      echo "itsxpress: '\$(itsxpress --version 2>/dev/null || echo n/a)'"
      echo "python: '\$(python --version 2>&1 | sed 's/^Python //' || echo n/a)'"
    } > software_versions.yml
    """
}
