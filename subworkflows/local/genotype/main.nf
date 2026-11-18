/*
 * CRAM -> single-sample, reblocked gVCF, for each SG to genotype: HaplotypeCaller over each genotyping interval,
 * then the SG's interval gVCFs merged and postprocessed. Ported from cpg-flow's genotype.py, without DRAGEN mode
 * and with contamination from CRAM QC (ADR 0001).
 *
 * Self-contained, so another pipeline (nf-align) can include it: everything comes in through meta and the inputs.
 */

include { MAKE_INTERVALS } from '../../../modules/local/make_intervals'
include { GENOTYPE_INTERVAL } from '../../../modules/local/genotype_interval'
include { POSTPROC_GVCF } from '../../../modules/local/postproc_gvcf'

workflow GENOTYPE {
    take:
    // [meta, cram, crai]: the CRAM and its index, as URI strings or Paths, read by GATK, not staged. meta needs id,
    // ref (the reference the CRAM was written against, with its .fai beside it), refDict (its .dict) and selfsm
    // (VerifyBamID's output, for the contamination)
    crams
    calling_regions  // the sequencing type's calling regions, as an interval list
    exclusions       // regions never genotyped, as an interval list, subtracted whatever N is
    interval_count   // N, the number of genotyping intervals, already resolved (never 0)
    noalt_bed        // the primary contigs plus chrM
    sequencing_type  // genome | exome
    downsample_reads // HaplotypeCaller's --max-reads-per-alignment-start, for genomes only
    pairhmm_threads
    gq_bands         // ReblockGVCF's GQ bands

    main:
    // GATK reads the CRAM by URI. A cloud Path's toString() drops its gs:// scheme, so a Path is turned into its URI
    def by_uri = crams.map { meta, cram, crai -> [meta, uriOf(cram), uriOf(crai)] }

    // once per call, and only if there's anything to genotype: every SG gets the same lists, named 0001, 0002...
    // so they sort into genomic order
    def intervals = MAKE_INTERVALS(
        calling_regions,
        exclusions,
        by_uri.first().map { _row -> interval_count },
        sequencing_type == 'genome' ? 100000 : 0,
    )

    def per_interval = GENOTYPE_INTERVAL(
        by_uri
            .combine(intervals.flatten())
            .map { meta, cram, crai, interval -> [meta, cram, crai, interval] + staged(meta) },
        sequencing_type == 'genome' ? downsample_reads : 0,
        pairhmm_threads,
    )

    // groupKey releases each SG as soon as its last interval finishes. An SG with a failed interval never
    // completes its group, and remainder: false drops it, so it fails rather than giving a partial gVCF. MAKE_INTERVALS
    // fails unless it made exactly N lists. The interval gVCFs are sorted by name, which is genomic order
    def grouped = per_interval
        .map { meta, gvcf, tbi -> [groupKey(meta, interval_count), gvcf, tbi] }
        .groupTuple(sort: { a, b -> a.name <=> b.name }, remainder: false)
        .map { key, gvcfs, tbis ->
            def meta = key.getGroupTarget()
            [meta, gvcfs, tbis] + staged(meta)
        }

    // the GATK version and contamination, recorded by the task, go into meta for the gvcf analysis
    def postprocessed = POSTPROC_GVCF(grouped, noalt_bed, gq_bands)
        .map { meta, gvcf, tbi, version, contamination ->
            [meta + [gatk_version: version.text.trim(), contamination: contamination.text.trim()], gvcf, tbi]
        }

    emit:
    gvcf = postprocessed // [meta, gvcf, tbi]
}

// what each task stages beside its gVCFs: the CRAM's own reference (FASTA, .fai, .dict) and the .selfSM. The .fai
// is resolved from the FASTA's Path, not by appending to a string, so a cloud Path keeps its scheme
def staged(meta) {
    def fasta = file(meta.ref)
    [fasta, fasta.resolveSibling("${fasta.name}.fai"), file(meta.refDict), file(meta.selfsm)]
}

// gs://bucket/path for a cloud Path (whose toString() drops the scheme), a plain path for a local one, and a
// string as it is
def uriOf(location) {
    if (!(location instanceof Path)) {
        return location.toString()
    }
    location.scheme == 'file' ? location.toString() : location.toUriString()
}
