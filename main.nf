nextflow.enable.dsl=2

include {
    RUN_FASTQC
    MULTIQC
    IMPORT_QIIME
    DEMUX_SUMMARIZE
    CUTADAPT_TRIM
    ITSXPRESS_TRIM
    DADA2_DENOISE
    CLASSIFY_TAXONOMY
    FILTER_TABLE
    SUMMARIZE_FEATURES
    TAXA_BARPLOT
    BUILD_PHYLOGENY
    BUILD_PHYLOGENY_SEPP
    DIVERSITY_ANALYSIS
    ALPHA_RAREFACTION
    ALPHA_GROUP_SIGNIFICANCE
    BETA_GROUP_SIGNIFICANCE
    LONGITUDINAL
    ANCOMBC2_SWEEP
    PLOT_REPORTS
    BUILD_REPORT
    DUMP_VERSIONS
} from './modules/local/metaqii'

import groovy.json.JsonOutput

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

def parseBool(value) {
    if (value == null)               return null
    if (value instanceof Boolean)    return value
    def s = value.toString().trim().toLowerCase()
    if (s in ['1', 'true',  't', 'yes', 'y']) return true
    if (s in ['0', 'false', 'f', 'no',  'n']) return false
    throw new IllegalArgumentException(
        "Cannot parse '${value}' as a boolean. Use true/false, yes/no, 1/0."
    )
}

def requireExistingPath(String name, def value) {
    if (!value) {
        throw new IllegalArgumentException("Missing required parameter: --${name}")
    }
    def target = file(value)
    if (!target.exists()) {
        throw new IllegalArgumentException(
            "Parameter --${name} points to a missing path: ${value}"
        )
    }
    return target
}

def requirePositiveInt(String name, def value) {
    if (value == null) {
        throw new IllegalArgumentException("Missing required parameter: --${name}")
    }
    def i = value as int
    if (i <= 0) {
        throw new IllegalArgumentException(
            "--${name} must be a positive integer. Received: ${value}"
        )
    }
    return i
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    // --- Input mode --------------------------------------------------------
    def inputMode = (params.input_mode ?: 'manifest').toString().trim().toLowerCase()
    if (!(inputMode in ['manifest', 'casava'])) {
        throw new IllegalArgumentException(
            "Unsupported --input_mode '${params.input_mode}'. Use 'manifest' or 'casava'."
        )
    }

    // --- Required paths ----------------------------------------------------
    def metadata   = requireExistingPath('metadata',   params.metadata)

    def inputSource = inputMode == 'manifest'
        ? requireExistingPath('manifest',  params.manifest)
        : requireExistingPath('input_dir', params.input_dir)

    // --- Classification method + its inputs --------------------------------
    def placeholder = file("${projectDir}/assets/NO_PHYLOGENY")
    def classifyMethod = (params.classification_method ?: 'sklearn').toString().toLowerCase()
    def classifierFile, refReadsFile, refTaxonomyFile
    if (classifyMethod == 'sklearn') {
        classifierFile  = requireExistingPath('classifier', params.classifier)
        refReadsFile    = placeholder
        refTaxonomyFile = placeholder
    }
    else if (classifyMethod == 'vsearch') {
        refReadsFile    = requireExistingPath('reference_reads',    params.reference_reads)
        refTaxonomyFile = requireExistingPath('reference_taxonomy', params.reference_taxonomy)
        classifierFile  = placeholder
    }
    else {
        throw new IllegalArgumentException(
            "Unsupported --classification_method '${params.classification_method}'. Use 'sklearn' or 'vsearch'."
        )
    }

    // --- Required scalars --------------------------------------------------
    def samplingDepth = requirePositiveInt('sampling_depth', params.sampling_depth)

    if (!params.metadata_column) {
        throw new IllegalArgumentException('Missing required parameter: --metadata_column')
    }

    // --- Taxonomy level ----------------------------------------------------
    def startLevel = (params.taxonomy_start_level ?: 2) as int
    if (startLevel < 1 || startLevel > 7) {
        throw new IllegalArgumentException(
            "--taxonomy_start_level must be between 1 and 7. Received: ${startLevel}"
        )
    }

    // --- Boolean flags -----------------------------------------------------
    def pairedEnd      = parseBool(params.paired_end)      ?: false
    def runFastqc      = parseBool(params.run_fastqc)      ?: false
    def isIts          = (params.sequence_type ?: 'ITS').toString().equalsIgnoreCase('ITS')
    def buildPhylogeny = params.build_phylogeny == null
        ? !isIts
        : (parseBool(params.build_phylogeny) ?: false)

    // Contaminant filtering (mitochondria/chloroplast) is meaningful for 16S
    // but not for ITS, so default it on only for non-ITS runs.
    def filterContaminants = params.filter_contaminants == null
        ? !isIts
        : (parseBool(params.filter_contaminants) ?: false)

    // --- Phylogeny method + SEPP reference ---------------------------------
    def phylogenyMethod = (params.phylogeny_method ?: 'denovo').toString().toLowerCase()
    if (!(phylogenyMethod in ['denovo', 'sepp'])) {
        throw new IllegalArgumentException(
            "Unsupported --phylogeny_method '${params.phylogeny_method}'. Use 'denovo' or 'sepp'."
        )
    }
    def seppReferenceFile = placeholder
    if (phylogenyMethod == 'sepp') {
        if (!buildPhylogeny) {
            log.warn '[metaQII] --phylogeny_method sepp requires build_phylogeny; no tree will be built while it is disabled.'
        } else {
            seppReferenceFile = requireExistingPath('sepp_reference', params.sepp_reference)
        }
    }

    // --- Primer removal (cutadapt) -----------------------------------------
    def runCutadapt = (params.fwd_primer != null && params.fwd_primer.toString().trim())
    if (runCutadapt && pairedEnd && !(params.rev_primer && params.rev_primer.toString().trim())) {
        throw new IllegalArgumentException(
            'Paired-end primer removal requires both --fwd_primer and --rev_primer.'
        )
    }

    // --- Longitudinal (optional) -------------------------------------------
    def runLongitudinal = (params.longitudinal_state_column && params.longitudinal_individual_column)

    // --- Alpha-rarefaction max depth ---------------------------------------
    // The rarefaction curve should extend well past the diversity sampling
    // depth so users can see whether richness saturates. Fall back to the
    // sampling depth when unset, but warn because that hides the plateau.
    def rarefactionMaxDepth = (params.alpha_rarefaction_max_depth ?: samplingDepth) as int
    if (params.alpha_rarefaction_max_depth == null) {
        log.warn """\
            [metaQII] --alpha_rarefaction_max_depth is unset; defaulting to the
            sampling depth (${samplingDepth}). For a useful rarefaction curve set
            it near the maximum per-sample frequency shown in 03_summaries/table.qzv.
        """.stripIndent()
    }

    // --- DADA2 truncation guard (paired-end) --------------------------------
    if (pairedEnd) {
        def tf = (params.trunc_len_f ?: 0) as int
        def tr = (params.trunc_len_r ?: 0) as int
        if (tf == 0 || tr == 0) {
            log.warn """\
                [metaQII] WARNING: --trunc_len_f (${tf}) or --trunc_len_r (${tr}) is 0
                for a paired-end run. DADA2 will not truncate reads, which can reduce
                denoising quality. Set both values explicitly unless your reads are
                already trimmed to a fixed length. Inspect 03_summaries/demux.qzv to
                choose truncation positions.
            """.stripIndent()
        }
    }

    // --- Runtime log -------------------------------------------------------
    log.info """
    =====================================================
     metaQII-nf
    =====================================================
     input mode          : ${inputMode}
     input source        : ${inputSource}
     paired-end          : ${pairedEnd}
     sequence type       : ${params.sequence_type}
     ITS trimming        : ${isIts}
     ITS region          : ${isIts ? params.its_region : 'N/A'}
     ITS taxa            : ${isIts ? params.its_taxa   : 'N/A'}
     build phylogeny     : ${buildPhylogeny}
     phylogeny method    : ${buildPhylogeny ? phylogenyMethod : 'N/A'}
     primer removal      : ${runCutadapt ? "fwd=${params.fwd_primer}${pairedEnd ? " rev=${params.rev_primer}" : ''}" : 'off'}
     classification      : ${classifyMethod}
     decontam controls   : ${params.decontam_control_column ?: 'off'}
     filter contaminants : ${filterContaminants}
     longitudinal        : ${runLongitudinal ? "${params.longitudinal_state_column} / ${params.longitudinal_individual_column}" : 'off'}
     run FastQC          : ${runFastqc}
     sampling depth      : ${samplingDepth}
     rarefaction max     : ${rarefactionMaxDepth}
     metadata column     : ${params.metadata_column}
     taxonomy start      : ${startLevel}
     DADA2 pooling       : ${params.dada2_pooling}
     trunc_len f/r       : ${params.trunc_len_f} / ${params.trunc_len_r}
     trim_left f/r       : ${params.trim_left_f} / ${params.trim_left_r}
     threads             : ${params.threads}
     output directory    : ${params.outdir}
    =====================================================
    """.stripIndent()

    // --- Channels ----------------------------------------------------------
    reads_input_ch     = Channel.of([inputMode, inputSource])
    metadata_ch        = Channel.value(metadata)
    classifier_ch      = Channel.value(classifierFile)
    refReads_ch        = Channel.value(refReadsFile)
    refTaxonomy_ch     = Channel.value(refTaxonomyFile)
    seppReference_ch   = Channel.value(seppReferenceFile)
    pairedEndCh        = Channel.value(pairedEnd)
    buildPhylogenyCh   = Channel.value(buildPhylogeny)
    filterContamCh     = Channel.value(filterContaminants)
    rarefactionMaxCh   = Channel.value(rarefactionMaxDepth)

    // --- FastQC + MultiQC (optional) ---------------------------------------
    if (runFastqc) {
        fastqc = RUN_FASTQC(reads_input_ch)
        if (parseBool(params.run_multiqc) ?: false) {
            MULTIQC(fastqc)
        }
    }

    // --- Import ------------------------------------------------------------
    imported = IMPORT_QIIME(reads_input_ch, pairedEndCh)

    // --- Primer removal + ITS trimming (both optional) ---------------------
    demuxForDenoising = imported.demux
    demuxSummaryCh    = imported.demux.map { qza -> tuple('demux', qza) }
    if (runCutadapt) {
        primerTrimmed     = CUTADAPT_TRIM(demuxForDenoising, pairedEndCh)
        demuxForDenoising = primerTrimmed.trimmed
        demuxSummaryCh    = demuxSummaryCh.mix(
            primerTrimmed.trimmed.map { qza -> tuple('primer-trimmed', qza) }
        )
    }
    if (isIts) {
        trimmed           = ITSXPRESS_TRIM(demuxForDenoising, pairedEndCh)
        demuxForDenoising = trimmed.trimmed_demux
        demuxSummaryCh    = demuxSummaryCh.mix(
            trimmed.trimmed_demux.map { qza -> tuple('trimmed-demux', qza) }
        )
    }

    // --- Quality summary of (trimmed) demultiplexed reads ------------------
    demuxSummaries = DEMUX_SUMMARIZE(demuxSummaryCh)

    // --- Denoising ---------------------------------------------------------
    denoised = DADA2_DENOISE(demuxForDenoising, pairedEndCh)

    // --- Taxonomy (needed before contaminant filtering) --------------------
    taxonomy = CLASSIFY_TAXONOMY(
        denoised.rep_seqs,
        classifier_ch,
        refReads_ch,
        refTaxonomy_ch
    )

    // --- Feature filtering (decontam + contaminants + prevalence) ----------
    filtered = FILTER_TABLE(
        denoised.table,
        denoised.rep_seqs,
        taxonomy.taxonomy_qza,
        metadata_ch,
        filterContamCh
    )

    // --- Summaries on the filtered feature table ---------------------------
    summaries = SUMMARIZE_FEATURES(
        filtered.table,
        filtered.rep_seqs,
        denoised.denoising_stats,
        metadata_ch
    )
    barplot = TAXA_BARPLOT(filtered.table, taxonomy.taxonomy_qza, metadata_ch)

    // --- Phylogeny (optional; de-novo or SEPP) -----------------------------
    rootedTreeCh = Channel.value(file("${projectDir}/assets/NO_PHYLOGENY"))
    if (buildPhylogeny) {
        if (phylogenyMethod == 'sepp') {
            phylogeny    = BUILD_PHYLOGENY_SEPP(filtered.rep_seqs, seppReference_ch)
            rootedTreeCh = phylogeny.rooted_tree
        } else {
            phylogeny    = BUILD_PHYLOGENY(filtered.rep_seqs)
            rootedTreeCh = phylogeny.rooted_tree
        }
    }

    // --- Diversity ---------------------------------------------------------
    diversity = DIVERSITY_ANALYSIS(
        filtered.table,
        metadata_ch,
        rootedTreeCh,
        buildPhylogenyCh
    )
    rarefaction = ALPHA_RAREFACTION(
        filtered.table,
        metadata_ch,
        rootedTreeCh,
        buildPhylogenyCh,
        rarefactionMaxCh
    )
    ALPHA_GROUP_SIGNIFICANCE(diversity.core_metrics_dir, metadata_ch)
    BETA_GROUP_SIGNIFICANCE(diversity.core_metrics_dir,  metadata_ch)

    // --- Longitudinal (optional) -------------------------------------------
    if (runLongitudinal) {
        LONGITUDINAL(diversity.core_metrics_dir, metadata_ch)
    }

    // --- Differential abundance (ANCOM-BC2) --------------------------------
    ANCOMBC2_SWEEP(filtered.table, taxonomy.taxonomy_qza, metadata_ch)

    // --- Plots -------------------------------------------------------------
    plots = PLOT_REPORTS(
        diversity.core_metrics_dir,
        filtered.table,
        taxonomy.taxonomy_qza,
        metadata_ch
    )

    // --- Consolidated HTML report ------------------------------------------
    report_qzvs = demuxSummaries.qzv
        .mix(summaries.denoising_stats_qzv)
        .mix(summaries.table_qzv)
        .mix(summaries.rep_seqs_qzv)
        .mix(taxonomy.taxonomy_qzv)
        .mix(barplot.qzv)
        .mix(rarefaction.qzv)
        .collect()

    BUILD_REPORT(report_qzvs, plots.plots)

    // --- Provenance: record tool versions ----------------------------------
    DUMP_VERSIONS()
}

// ---------------------------------------------------------------------------
// Run summary + resolved parameters (for reproducibility)
// ---------------------------------------------------------------------------

workflow.onComplete {
    def infoDir = file("${params.outdir}/pipeline_info")
    infoDir.mkdirs()

    file("${infoDir}/params.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(params))

    file("${infoDir}/run_summary.txt").text = """\
        metaQII-nf run summary
        ======================
        pipeline version : ${workflow.manifest.version}
        run name         : ${workflow.runName}
        session id       : ${workflow.sessionId}
        started          : ${workflow.start}
        completed        : ${workflow.complete}
        duration         : ${workflow.duration}
        success          : ${workflow.success}
        exit status      : ${workflow.exitStatus}
        command line     : ${workflow.commandLine}
        nextflow version : ${workflow.nextflow.version}
        container        : ${params.container}
        output directory : ${params.outdir}
    """.stripIndent()

    log.info "[metaQII] ${workflow.success ? 'Completed successfully' : 'Finished with errors'} — outputs in ${params.outdir}"
}
