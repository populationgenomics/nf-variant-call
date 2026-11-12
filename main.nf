/*
 * Genotyping of one sequencing type: Metamist cohorts or SGs, from any number of datasets -> sample sheet ->
 * GATK HaplotypeCaller over N genotyping intervals -> one reblocked gVCF per SG, with gVCF QC, published to and
 * registered in each SG's own dataset. Run records go to the analysis dataset (params.dataset). With
 * --qc_genotyped, already-genotyped SGs without gVCF QC get gVCF QC alone. See docs/DESIGN.md.
 */

include { BUILD_SAMPLESHEET } from './modules/local/build_samplesheet'
include { MAKE_INTERVALS } from './modules/local/make_intervals'
include { PREFLIGHT } from './modules/local/preflight'
include { GENOTYPE } from './subworkflows/local/genotype'
include { GVCF_QC } from './subworkflows/local/gvcf_qc'

/*
 * Types only: the defaults are in nextflow.config, which conf/base.config reads before this script is parsed.
 * Nextflow converts command-line values to these types (`--register false` is a Boolean), and fails the run if a
 * param without `?` has no value.
 */
params {
    dataset: String
    sequencing_type: String
    accessLevel: String
    arGuid: String
    cohorts: String? // or a list, from a config file
    sgids: String?   // or a list, from a config file
    output_root: String?

    register: Boolean
    qc_genotyped: Boolean
    register_wait_minutes: Integer
    genotyping_intervals: Integer
    downsample_reads: Integer
    pairhmm_threads: Integer
    reblock_gq_bands: List<Integer>

    aligned_ref_fasta: String
    genotyped_ref_fasta: String
    calling_regions_genome: String
    calling_regions_exome: String
    calling_exclusions: String
    noalt_bed: String
    dbsnp_vcf: String
    genome_evaluation_regions: String
    exome_evaluation_regions: String

    gatk_container: String
    picard_container: String
    cpg_container: String
}

/*
 * Output destinations. Each SG's are worked out here and carried in its meta: publishDir (conf/base.config),
 * the pre-flight check and registration all read them from there, so they can't drift apart.
 * Exomes live under exome/, as in cpg-flow.
 */
def typePrefix() {
    params.sequencing_type == 'exome' ? 'exome/' : ''
}

def metamistProject(ds) {
    params.accessLevel == 'test' ? "${ds}-test" : ds
}

// a dataset's bucket, or <output_root>/<ds> for local runs (absolute, whatever the launch dir)
def datasetRoot(ds) {
    params.output_root
        ? uriOf(file(params.output_root).resolve(ds))
        : "gs://cpg-${ds}-${params.accessLevel == 'test' ? 'test' : 'main'}"
}

// per-run records, in the analysis dataset: the sample sheet
def runDir() {
    "${datasetRoot(params.dataset)}/genotyping_samplesheets/${params.arGuid}"
}

// the files of an SG's gvcf analysis, by registration key: the gVCF, then its index
def gvcfFiles(meta) {
    [gvcf: meta.gvcf, index: "${meta.gvcf}.tbi".toString()]
}

// the gVCF QC files this run publishes for an SG, the secondary files of its gvcf_qc analysis
def gvcfQcFiles(meta) {
    def names = [
        variant_calling_summary_metrics: 'variant_calling_summary_metrics',
        variant_calling_detail_metrics: 'variant_calling_detail_metrics',
        bcftools_stats: 'bcftools-stats',
    ]
    names.collectEntries { key, suffix -> [(key): "${meta.qcDir}/${meta.id}.${suffix}".toString()] }
}

// every file this run publishes for an SG, which the pre-flight check makes sure don't exist yet. A genotyped
// SG's gVCF is already registered, so at most its gVCF QC is published
def outputs(meta) {
    if (meta.kind == 'to_genotype') {
        return gvcfFiles(meta) + gvcfQcFiles(meta)
    }
    meta.kind == 'qc_only' ? gvcfQcFiles(meta) : [:]
}

// the .dict of a reference: the <name>.dict beside <name>.fasta, where Picard and GATK look for it. null for a
// FASTA without a FASTA extension
def refDict(fasta) {
    def dict = fasta.replaceFirst(/\.(fasta|fa|fna)(\.gz)?$/, '.dict')
    dict == fasta ? null : dict
}

// gs://bucket/path for cloud paths (Path.toString() drops the scheme), a plain path locally
def uriOf(path) {
    path.scheme == 'file' ? path.toString() : path.toUriString()
}

// a list from config, or a comma/space separated string from the command line
def idList(param) {
    (param instanceof List) ? param : param.toString().tokenize(', ')
}

// what the params block can't check: allowed values, and params that are needed together or not at all
def validateParams() {
    if (!(params.sequencing_type in ['genome', 'exome'])) {
        error("--sequencing_type must be genome or exome, got '${params.sequencing_type}'")
    }
    if (!(params.accessLevel in ['full', 'standard', 'test'])) {
        error("--accessLevel must be full, standard or test, got '${params.accessLevel}'")
    }
    if ([params.cohorts, params.sgids].count { input -> input } != 1) {
        error('Provide exactly one of --cohorts or --sgids')
    }
    if (params.genotyping_intervals < 0) {
        error("--genotyping_intervals must be 0 (the default for the sequencing type) or more, got ${params.genotyping_intervals}")
    }
    ['downsample_reads', 'pairhmm_threads', 'register_wait_minutes'].each { name ->
        if (params[name] < 1) {
            error("--${name} must be at least 1, got ${params[name]}")
        }
    }
    if (!params.reblock_gq_bands || params.reblock_gq_bands.any { band -> band < 1 }) {
        error("--reblock_gq_bands must be a list of positive GQs, got ${params.reblock_gq_bands}")
    }
}

/*
 * Fail before the sample sheet is built if this arGuid's run records already exist: publishDir would keep the
 * old ones silently. arGuid is unique in production, but set by hand while testing. A resumed run finds its
 * own records there, so -resume skips the check.
 */
def checkNoRunRecords() {
    if (workflow.resume) {
        return
    }
    def sheet = "${runDir()}/samplesheet.tsv"
    if (file(sheet).exists()) {
        error("Run records for arGuid ${params.arGuid} already exist. Launch with a new --arGuid:\n  ${sheet}")
    }
}

/*
 * A sample sheet row -> meta: the SG's destinations, and what its tasks read. kind is to_genotype, genotyped, or
 * qc_only (a genotyped SG whose gVCF has no gvcf_qc analysis, with --qc_genotyped). ref is the reference the SG's
 * tasks read: the CRAM's for an SG to genotype, which its new gVCF is then written against, or the registered
 * gVCF's otherwise, and refDict its .dict. gvcf is the registered gVCF, or where this run publishes the new one.
 */
def parseRow(r) {
    // decides where the SG is published and registered, and which dataset its compute is billed to
    if (!r.dataset) {
        fail("${r.sg_id}: a row needs dataset, the dataset the SG belongs to")
    }
    def genotyped = r.gvcf as boolean
    if (!genotyped) {
        def missing = ['cram', 'cram_ref', 'selfsm'].findAll { column -> !r[column] }
        if (missing) {
            fail("${r.sg_id}: a row to genotype needs ${missing.join(', ')}")
        }
    }
    def kind = !genotyped ? 'to_genotype' : (params.qc_genotyped && !r.gvcf_qc ? 'qc_only' : 'genotyped')
    def ref = genotyped ? r.gvcf_ref : r.cram_ref
    // GATK and Picard need the <name>.dict beside the reference. Only checked where a task reads it
    if (kind != 'genotyped' && !refDict(ref)) {
        fail("${r.sg_id}: reference '${ref}' must end in .fasta, .fa or .fna, with its .dict beside it")
    }
    def base = "${datasetRoot(r.dataset)}/${typePrefix()}"
    [
        id: r.sg_id,
        kind: kind,
        dataset: r.dataset,
        project: metamistProject(r.dataset),
        cram: r.cram ?: null,
        ref: ref,
        refDict: refDict(ref),
        selfsm: r.selfsm ?: null,
        gvcf: genotyped ? r.gvcf : "${base}gvcf/${r.sg_id}.g.vcf.gz".toString(),
        gvcfDir: "${base}gvcf".toString(),
        qcDir: "${base}gvcf_qc/${r.sg_id}".toString(),
    ]
}

// error() inside a channel operator reaches the console only as "Unexpected error", so log the reason first
def fail(message) {
    log.error(message)
    error(message)
}

workflow {
    validateParams()
    checkNoRunRecords()
    def genome = params.sequencing_type == 'genome'
    // matches {genome,exome}_scatter_count in cpg-flow's config_template.toml
    def interval_count = params.genotyping_intervals ?: (genome ? 30 : 1)

    // every run starts from a sheet built from Metamist: SGs are discovered through cohorts, or queried directly,
    // and one of the two lists is empty. A resumed run takes the sheet from the cache, not from Metamist, which
    // would by now see the SGs its first attempt registered as genotyped. So it carries on with the same SGs, and
    // the cached pre-flight check lets its own published files through
    def cohorts = params.cohorts ? idList(params.cohorts) : []
    def sgids = params.sgids ? idList(params.sgids) : []
    def sheet = BUILD_SAMPLESHEET(
        params.accessLevel,
        params.sequencing_type,
        cohorts,
        sgids,
        params.aligned_ref_fasta,
        params.genotyped_ref_fasta,
        runDir(),
    )

    // hold every SG that publishes anything until the whole sheet has passed the pre-flight check, so a clash fails
    // the run before any work starts. Genotyped SGs that aren't QC-only publish nothing, and go no further
    def rows = PREFLIGHT(
        sheet
            .splitCsv(header: true, sep: '\t')
            .map { r -> parseRow(r) }
            .filter { meta -> meta.kind != 'genotyped' }
            .toList()
            .map { metas -> [metas, metas.collectMany { meta -> outputs(meta).values() }.unique()] }
    )
        .flatMap { metas -> metas }
        .branch { meta ->
            to_genotype: meta.kind == 'to_genotype'
            qc_only: meta.kind == 'qc_only'
        }

    // once per run, and only if there's anything to genotype: every SG gets the same lists
    def intervals = MAKE_INTERVALS(
        file(genome ? params.calling_regions_genome : params.calling_regions_exome),
        file(params.calling_exclusions),
        rows.to_genotype.first().map { _meta -> interval_count },
        genome ? 100000 : 0,
    )

    // each CRAM is read from GCS through its index, as nf-align publishes it, not staged
    def genotyped = GENOTYPE(
        rows.to_genotype.map { meta -> [meta, meta.cram, "${meta.cram}.crai".toString()] },
        intervals,
        file(params.noalt_bed),
        params.sequencing_type,
        params.downsample_reads,
        params.pairhmm_threads,
        params.reblock_gq_bands,
    )

    // every new gVCF, against the reference its CRAM was written with, and each QC-only SG's registered gVCF,
    // against the one its gvcf analysis records (meta.ref, from the row's gvcf_ref)
    GVCF_QC(
        genotyped.gvcf.mix(rows.qc_only.map { meta -> [meta, file(meta.gvcf), file("${meta.gvcf}.tbi")] }),
        file(genome ? params.genome_evaluation_regions : params.exome_evaluation_regions),
        [file(params.dbsnp_vcf), file("${params.dbsnp_vcf}.idx")],
    )
}
