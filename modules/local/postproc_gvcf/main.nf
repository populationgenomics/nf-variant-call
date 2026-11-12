/*
 * One SG's interval gVCFs -> its single-sample, reblocked gVCF, with bin/postproc_gvcf.sh: merged if there's more
 * than one, then cpg-flow's postprocessing. Also emits the GATK version and the FREEMIX the SG was genotyped with,
 * for its gvcf analysis meta.
 */
process POSTPROC_GVCF {
    tag "${meta.id}"

    input:
    // the interval gVCFs in genomic order, and the CRAM's own reference, which ReblockGVCF reads
    tuple val(meta), path(gvcfs, stageAs: 'in/*'), path(tbis, stageAs: 'in/*'), path(fasta, stageAs: 'ref/*'), path(fai, stageAs: 'ref/*'), path(dict, stageAs: 'ref/*'), path(selfsm)
    path noalt_bed, stageAs: 'noalt/*'
    val gq_bands

    output:
    tuple val(meta), path("${meta.id}.g.vcf.gz"), path("${meta.id}.g.vcf.gz.tbi"), path('gatk_version.txt'), path('contamination.txt')

    script:
    """
    ${postprocCommand(meta, gvcfs, fasta, noalt_bed, selfsm, gq_bands, task.memory)}
    """

    stub:
    // the command, as it would run, to .command.out. It holds no single quotes
    """
    echo '${postprocCommand(meta, gvcfs, fasta, noalt_bed, selfsm, gq_bands, task.memory)}'
    touch ${meta.id}.g.vcf.gz ${meta.id}.g.vcf.gz.tbi
    echo stub > gatk_version.txt
    echo 0.01 > contamination.txt
    """
}

def postprocCommand(meta, gvcfs, fasta, noalt_bed, selfsm, gq_bands, memory) {
    // the heap, for MergeVcfs and ReblockGVCF, less headroom for the JVM's own overhead and bcftools
    def heap = Math.max(1, memory.toGiga() - 1)
    // one staged gVCF is a path, more are a list
    def inputs = gvcfs instanceof List ? gvcfs.join(' ') : gvcfs
    "postproc_gvcf.sh ${meta.id} ${fasta} ${noalt_bed} ${selfsm} ${gq_bands.join(',')} ${heap} ${inputs}"
}
