/*
 * CRAM -> single-sample, reblocked gVCF, for each SG to genotype: HaplotypeCaller over each genotyping interval,
 * then the SG's interval gVCFs merged and postprocessed. Ported from cpg-flow's genotype.py, without DRAGEN mode
 * and with contamination from CRAM QC (ADR 0001).
 *
 * Self-contained, so another pipeline (nf-align) can include it: everything comes in through meta and the inputs.
 */

include { GENOTYPE_INTERVAL } from '../../../modules/local/genotype_interval'
include { POSTPROC_GVCF } from '../../../modules/local/postproc_gvcf'

workflow GENOTYPE {
    take:
    // [meta, cram, crai]: the CRAM and its index as URIs, read by GATK, not staged. meta needs id, ref (the
    // reference the CRAM was written against, with its .fai beside it), refDict (its .dict) and selfsm
    // (VerifyBamID's output, for the contamination)
    crams
    intervals        // value: the genotyping interval lists, named so they sort into genomic order
    noalt_bed        // the primary contigs plus chrM
    sequencing_type  // genome | exome
    downsample_reads // HaplotypeCaller's --max-reads-per-alignment-start, for genomes only
    pairhmm_threads
    gq_bands         // ReblockGVCF's GQ bands

    main:
    def per_interval = GENOTYPE_INTERVAL(
        crams
            .combine(intervals.flatten())
            .map { meta, cram, crai, interval -> [meta, cram, crai, interval] + staged(meta) },
        sequencing_type == 'genome' ? downsample_reads : 0,
        pairhmm_threads,
    )

    // groupKey releases each SG as soon as its last interval finishes. An SG with a failed interval never
    // completes its group, and remainder: false drops it, so it fails rather than giving a partial gVCF. The
    // interval gVCFs are sorted by name, which is genomic order
    def interval_count = intervals.flatten().count()
    def grouped = per_interval
        .combine(interval_count)
        .map { meta, gvcf, tbi, count -> [groupKey(meta, count), gvcf, tbi] }
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

// what each task stages beside its gVCFs: the CRAM's own reference (FASTA, .fai, .dict) and the .selfSM
def staged(meta) {
    [file(meta.ref), file("${meta.ref}.fai"), file(meta.refDict), file(meta.selfsm)]
}
