/*
 * gVCF QC: Picard's variant calling metrics and bcftools stats of each gVCF, in parallel, one task per tool's image,
 * against the gVCF's own reference. An SG is emitted only once both have succeeded, so its gVCF QC is all or
 * nothing.
 *
 * Self-contained, so another pipeline (nf-align) can include it: everything comes in through meta and the inputs.
 */

include { VARIANT_CALLING_METRICS } from '../../../modules/local/variant_calling_metrics'
include { BCFTOOLS_STATS } from '../../../modules/local/bcftools_stats'

workflow GVCF_QC {
    take:
    // [meta, gvcf, tbi]. meta needs id, ref (the reference the gVCF was written against, with its .fai beside it)
    // and refDict (its .dict)
    gvcfs
    evaluation_regions // the sequencing type's evaluation regions, as BED
    dbsnp              // value: [vcf, idx]

    main:
    def metrics = VARIANT_CALLING_METRICS(
        gvcfs.map { meta, gvcf, tbi -> [meta, gvcf, tbi, file(meta.refDict)] },
        evaluation_regions,
        dbsnp,
    )
    // the .fai is resolved from the FASTA's Path, not by appending to a string, so a cloud Path keeps its scheme
    def stats = BCFTOOLS_STATS(
        gvcfs.map { meta, gvcf, tbi ->
            def fasta = file(meta.ref)
            [meta, gvcf, tbi, fasta, fasta.resolveSibling("${fasta.name}.fai")]
        },
    )

    // keyed by SG ID; join drops an SG missing from either
    def complete = metrics.map { meta, summary, detail -> [meta.id, meta, summary, detail] }
        .join(stats.map { meta, bcftools -> [meta.id, bcftools] })
        .map { _id, meta, summary, detail, bcftools -> [meta, summary, detail, bcftools] }

    emit:
    qc = complete // [meta, summary metrics, detail metrics, bcftools stats], for SGs whose every QC task succeeded
}
