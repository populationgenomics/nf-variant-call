// cohorts or SG IDs, from any datasets at this access level -> one genotyping sample sheet for this sequencing
// type; see bin/build_samplesheet.py
process BUILD_SAMPLESHEET {
    tag "${access_level}"

    input:
    val access_level
    val sequencing_type
    val cohorts
    val sgids
    val aligned_ref_fasta // the reference of a CRAM whose analysis doesn't record one
    val genotyped_ref_fasta // the reference of a gVCF whose analysis doesn't record one
    val run_dir // only read by publishDir (conf/base.config)

    output:
    path 'samplesheet.tsv'

    script:
    // exactly one of cohorts and sgids is non-empty
    def ids = cohorts ? "--cohorts ${cohorts.join(' ')}" : "--sgids ${sgids.join(' ')}"
    """
    build_samplesheet.py \\
        --access-level ${access_level} \\
        --sequencing-type ${sequencing_type} \\
        ${ids} \\
        --aligned-ref-fasta ${aligned_ref_fasta} \\
        --genotyped-ref-fasta ${genotyped_ref_fasta} \\
        --output samplesheet.tsv
    """

    stub:
    // stub runs have no Metamist, so a test config can hand over a prepared sheet; otherwise it's empty
    def prepared = task.ext.stub_samplesheet
    if (prepared) {
        "cp ${prepared} samplesheet.tsv"
    }
    else {
        """
        printf 'sg_id\\tdataset\\tcram\\tcram_ref\\tselfsm\\tgvcf\\tgvcf_ref\\tgvcf_qc\\n' > samplesheet.tsv
        """
    }
}
