/*
 * Genotyping of one sequencing type: Metamist cohorts or SGs, from any number of datasets -> sample sheet ->
 * GATK HaplotypeCaller over N genotyping intervals -> one reblocked gVCF per SG, with gVCF QC, published to and
 * registered in each SG's own dataset. Run records go to the analysis dataset (params.dataset). With
 * --qc_genotyped, already-genotyped SGs without gVCF QC get gVCF QC alone. See docs/DESIGN.md.
 */

include { BUILD_SAMPLESHEET } from './modules/local/build_samplesheet'

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

// error() inside a channel operator reaches the console only as "Unexpected error", so log the reason first
def fail(message) {
    log.error(message)
    error(message)
}

workflow {
    validateParams()
    checkNoRunRecords()

    // every run starts from a sheet built from Metamist: SGs are discovered through cohorts, or queried directly,
    // and one of the two lists is empty. A resumed run takes the sheet from the cache, not from Metamist, which
    // would by now see the SGs its first attempt registered as genotyped. So it carries on with the same SGs, and
    // the cached pre-flight check lets its own published files through
    def cohorts = params.cohorts ? idList(params.cohorts) : []
    def sgids = params.sgids ? idList(params.sgids) : []
    BUILD_SAMPLESHEET(
        params.accessLevel,
        params.sequencing_type,
        cohorts,
        sgids,
        params.aligned_ref_fasta,
        params.genotyped_ref_fasta,
        runDir(),
    )
}
