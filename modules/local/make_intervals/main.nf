/*
 * The sequencing type's calling regions, minus the exclusions, split into N genotyping intervals: once per run, and
 * every SG gets the same lists. Ported from cpg-flow's picard.py::generate_intervals, but the exclusions are
 * subtracted even when N is 1, so the regions genotyped never depend on N.
 * The lists are named 0001.interval_list, 0002.interval_list..., so they sort into genomic order. The task fails
 * unless there are exactly N, since each SG's interval gVCFs are grouped by that count.
 */
process MAKE_INTERVALS {
    tag "${count}"

    input:
    path calling_regions, stageAs: 'regions/*'
    path exclusions, stageAs: 'exclusions/*'
    val count             // N
    val break_bands_at    // IntervalListTools' -BREAK_BANDS_AT_MULTIPLES_OF: 100000 for genomes, 0 for exomes

    output:
    path '*.interval_list'

    script:
    // with a SCATTER_COUNT of 1, IntervalListTools writes OUTPUT as one list, otherwise a directory of
    // temp_<i>_of_<N>/scattered.interval_list, with i zero-padded to 4
    def output = count == 1 ? '0001.interval_list' : 'out'
    """
    picard -Xms1000m -Xmx1500m IntervalListTools \\
        -SCATTER_COUNT ${count} \\
        -SUBDIVISION_MODE INTERVAL_SUBDIVISION \\
        -UNIQUE true \\
        -SORT true \\
        -BREAK_BANDS_AT_MULTIPLES_OF ${break_bands_at} \\
        -I ${calling_regions} \\
        -ACTION SUBTRACT -SI ${exclusions} \\
        -OUTPUT ${output}

    if [ -d out ]; then
        for dir in out/temp_*_of_*; do
            mv \$dir/scattered.interval_list \$(basename \$dir | cut -d_ -f2).interval_list
        done
    fi

    made=\$(ls *.interval_list | wc -l)
    if [ \$made -ne ${count} ]; then
        echo "IntervalListTools made \$made interval lists, not ${count}" >&2
        exit 1
    fi
    """

    stub:
    """
    for i in \$(seq 1 ${count}); do touch \$(printf '%04d' \$i).interval_list; done
    """
}
