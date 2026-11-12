/*
 * bcftools stats of one gVCF's variant sites, with bin/bcftools_stats.sh: reference blocks are dropped, and records
 * whose ALT is like C,<NON_REF> are kept. Against the gVCF's own reference.
 */
process BCFTOOLS_STATS {
    tag "${meta.id}"

    input:
    // the gVCF's own reference, staged apart so it can't collide with the gVCF
    tuple val(meta), path(gvcf), path(tbi), path(fasta, stageAs: 'ref/*'), path(fai, stageAs: 'ref/*')

    output:
    tuple val(meta), path("${meta.id}.bcftools-stats")

    script:
    """
    ${statsCommand(meta, gvcf, fasta)}
    """

    stub:
    // the command, as it would run, to .command.out
    """
    echo '${statsCommand(meta, gvcf, fasta)}'
    touch ${meta.id}.bcftools-stats
    """
}

def statsCommand(meta, gvcf, fasta) {
    "bcftools_stats.sh ${gvcf} ${fasta} ${meta.id}.bcftools-stats"
}
