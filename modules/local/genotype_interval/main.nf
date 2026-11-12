/*
 * One SG's CRAM, over one genotyping interval -> that interval's gVCF. Ported from cpg-flow's genotype.py, without
 * --dragen-mode and with contamination from the SG's CRAM QC (ADR 0001).
 *
 * The CRAM is not staged: its URI goes in as a val, and GATK reads it from GCS through its index, so each task
 * fetches only its own interval's reads (RISKS #3). The reference, the interval list and the .selfSM are staged.
 */
process GENOTYPE_INTERVAL {
    tag "${meta.id}:${interval.baseName}"

    input:
    // the CRAM's own reference, staged apart so it can't collide with other inputs. GATK finds the .dict beside
    // the FASTA, by name
    tuple val(meta), val(cram), val(crai), path(interval, stageAs: 'intervals/*'), path(fasta, stageAs: 'ref/*'), path(fai, stageAs: 'ref/*'), path(dict, stageAs: 'ref/*'), path(selfsm)
    val downsample_reads // --max-reads-per-alignment-start, or 0 for none (exomes)
    val pairhmm_threads

    output:
    tuple val(meta), path("${meta.id}.${interval.baseName}.g.vcf.gz"), path("${meta.id}.${interval.baseName}.g.vcf.gz.tbi")

    script:
    """
    # on its own line, so set -e fails the task if the .selfSM can't be read
    FREEMIX=\$(freemix.sh ${selfsm})

    ${haplotypeCaller(meta, cram, crai, interval, fasta, downsample_reads, pairhmm_threads, task.memory)}
    """

    stub:
    // the command, as it would run, to .command.out. It holds no single quotes
    """
    echo '${haplotypeCaller(meta, cram, crai, interval, fasta, downsample_reads, pairhmm_threads, task.memory)}'
    touch ${meta.id}.${interval.baseName}.g.vcf.gz ${meta.id}.${interval.baseName}.g.vcf.gz.tbi
    """
}

// cpg-flow's command, for the script and for the stub to print. -L is passed even when there's one interval
def haplotypeCaller(meta, cram, crai, interval, fasta, downsample_reads, pairhmm_threads, memory) {
    // the heap, less headroom for the native PairHMM and the JVM's own overhead
    def heap = Math.max(1, memory.toGiga() - 2)
    def downsample = downsample_reads ? " \\\n    --max-reads-per-alignment-start ${downsample_reads}" : ''
    """\
gatk --java-options "-Xms${heap}g -Xmx${heap}g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \\
    HaplotypeCaller \\
    -R ${fasta} \\
    -I ${cram} \\
    --read-index ${crai} \\
    -L ${interval} \\
    --disable-spanning-event-genotyping \\
    --contamination-fraction-to-filter \$FREEMIX \\
    --native-pair-hmm-threads ${pairhmm_threads} \\
    -O ${meta.id}.${interval.baseName}.g.vcf.gz \\
    -G AS_StandardAnnotation \\
    -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \\
    -ERC GVCF \\
    --create-output-variant-index${downsample}"""
}
