/*
 * CRAM -> single-sample, reblocked gVCF, for each SG to genotype: HaplotypeCaller over each genotyping interval,
 * then the SG's interval gVCFs merged and postprocessed. Ported from cpg-flow's genotype.py, without DRAGEN mode
 * and with contamination from CRAM QC (ADR 0001).
 *
 * Self-contained, so another pipeline (nf-align) can include it: everything comes in through meta and the inputs.
 */

include { GENOTYPE_INTERVAL } from '../../../modules/local/genotype_interval'

workflow GENOTYPE {
    take:
    // [meta, cram, crai]: the CRAM and its index as URIs, read by GATK, not staged. meta needs id, ref (the
    // reference the CRAM was written against, with its .fai and <name>.dict beside it), refDict (that .dict) and
    // selfsm (VerifyBamID's output, for the contamination)
    crams
    intervals        // value: the genotyping interval lists, in genomic order
    sequencing_type  // genome | exome
    downsample_reads // HaplotypeCaller's --max-reads-per-alignment-start, for genomes only
    pairhmm_threads

    main:
    def per_interval = GENOTYPE_INTERVAL(
        crams
            .combine(intervals.flatten())
            .map { meta, cram, crai, interval ->
                [meta, cram, crai, interval, file(meta.ref), file("${meta.ref}.fai"), file(meta.refDict), file(meta.selfsm)]
            },
        sequencing_type == 'genome' ? downsample_reads : 0,
        pairhmm_threads,
    )

    emit:
    interval_gvcfs = per_interval // [meta, gvcf, tbi], one per SG per interval
}
