/*
 * Picard CollectVariantCallingMetrics of one gVCF, as cpg-flow's vcf_qc, over the sequencing type's evaluation
 * regions. They're stored as BED, so they aren't tied to one reference, and turned into an interval list headed
 * with the gVCF's own reference .dict here, as nf-align's PICARD_METRICS does, so the md5s always match.
 */
process VARIANT_CALLING_METRICS {
    tag "${meta.id}"

    input:
    // the .dict of the gVCF's own reference, staged apart so it can't collide with other inputs
    tuple val(meta), path(gvcf), path(tbi), path(dict, stageAs: 'ref/*')
    path evaluation_regions, stageAs: 'intervals/*'
    tuple path(dbsnp, stageAs: 'dbsnp/*'), path(dbsnp_idx, stageAs: 'dbsnp/*')

    output:
    tuple val(meta), path("${meta.id}.variant_calling_summary_metrics"), path("${meta.id}.variant_calling_detail_metrics")

    script:
    """
    ${metricsCommand(meta, gvcf, dict, evaluation_regions, dbsnp, task.memory)}
    """

    stub:
    // the command, as it would run, to .command.out. It holds no single quotes
    """
    echo '${metricsCommand(meta, gvcf, dict, evaluation_regions, dbsnp, task.memory)}'
    touch ${meta.id}.variant_calling_summary_metrics ${meta.id}.variant_calling_detail_metrics
    """
}

def metricsCommand(meta, gvcf, dict, evaluation_regions, dbsnp, memory) {
    // the heap, less headroom for the JVM's own overhead
    def heap = "-Xmx${Math.max(1, memory.toGiga() - 1)}g"
    """\
picard BedToIntervalList -I ${evaluation_regions} -O evaluation.interval_list -SD ${dict}

picard ${heap} CollectVariantCallingMetrics \\
    -I ${gvcf} \\
    -O ${meta.id} \\
    --DBSNP ${dbsnp} \\
    -SD ${dict} \\
    -TI evaluation.interval_list \\
    --GVCF_INPUT true"""
}
