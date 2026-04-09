nextflow.enable.dsl=2

include {
    RUN_FASTQC
    IMPORT_QIIME
    ITSXPRESS_TRIM
    DADA2_DENOISE
    SUMMARIZE_FEATURES
    BUILD_PHYLOGENY
    DIVERSITY_ANALYSIS
    ALPHA_RAREFACTION
    ALPHA_GROUP_SIGNIFICANCE
    BETA_GROUP_SIGNIFICANCE
    CLASSIFY_TAXONOMY
    ANCOM_SWEEP
    PLOT_REPORTS
} from './modules/local/metaqii'

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
    def classifier = requireExistingPath('classifier', params.classifier)

    def inputSource = inputMode == 'manifest'
        ? requireExistingPath('manifest',  params.manifest)
        : requireExistingPath('input_dir', params.input_dir)

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

    // --- DADA2 truncation guard (paired-end) --------------------------------
    if (pairedEnd) {
        def tf = (params.trunc_len_f ?: 0) as int
        def tr = (params.trunc_len_r ?: 0) as int
        if (tf == 0 || tr == 0) {
            log.warn """\
                [metaQII] WARNING: --trunc_len_f (${tf}) or --trunc_len_r (${tr}) is 0
                for a paired-end run. DADA2 will not truncate reads, which can reduce
                denoising quality. Set both values explicitly unless your reads are
                already trimmed to a fixed length.
            """.stripIndent()
        }
    }

    // Expose resolved paired_end back to params so modules can read it
    params.paired_end = pairedEnd

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
     run FastQC          : ${runFastqc}
     sampling depth      : ${samplingDepth}
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
    reads_input_ch = Channel.of([inputMode, inputSource])
    metadata_ch    = Channel.value(metadata)
    classifier_ch  = Channel.value(classifier)

    // --- FastQC (optional) -------------------------------------------------
    if (runFastqc) {
        RUN_FASTQC(reads_input_ch)
    }

    // --- Import ------------------------------------------------------------
    imported = IMPORT_QIIME(reads_input_ch)

    // --- ITS trimming or passthrough ---------------------------------------
    demuxForDenoising = imported.demux
    if (isIts) {
        trimmed           = ITSXPRESS_TRIM(imported.demux)
        demuxForDenoising = trimmed.trimmed_demux
    }

    // --- Denoising & summaries --------------------------------------------
    denoised = DADA2_DENOISE(demuxForDenoising)
    SUMMARIZE_FEATURES(
        denoised.table,
        denoised.rep_seqs,
        denoised.denoising_stats,
        metadata_ch
    )

    // --- Phylogeny (optional) ---------------------------------------------
    rootedTreeCh = Channel.value(file("${projectDir}/assets/NO_PHYLOGENY"))
    if (buildPhylogeny) {
        phylogeny    = BUILD_PHYLOGENY(denoised.rep_seqs)
        rootedTreeCh = phylogeny.rooted_tree
    }
    buildPhylogenyCh = Channel.value(buildPhylogeny)

    // --- Diversity ---------------------------------------------------------
    diversity = DIVERSITY_ANALYSIS(
        denoised.table,
        metadata_ch,
        rootedTreeCh,
        buildPhylogenyCh
    )
    ALPHA_RAREFACTION(
        denoised.table,
        metadata_ch,
        rootedTreeCh,
        buildPhylogenyCh
    )
    ALPHA_GROUP_SIGNIFICANCE(diversity.core_metrics_dir, metadata_ch)
    BETA_GROUP_SIGNIFICANCE(diversity.core_metrics_dir,  metadata_ch)

    // --- Taxonomy & ANCOM -------------------------------------------------
    taxonomy = CLASSIFY_TAXONOMY(
        denoised.rep_seqs,
        denoised.table,
        metadata_ch,
        classifier_ch
    )
    ANCOM_SWEEP(denoised.table, taxonomy.taxonomy_qza, metadata_ch)

    // --- Plots ------------------------------------------------------------
    PLOT_REPORTS(
        diversity.core_metrics_dir,
        denoised.table,
        taxonomy.taxonomy_qza,
        metadata_ch
    )
}
