# Risks

1. **Mixed genotyping modes in one VDS.** Every existing gVCF is DRAGEN-mode, and new ones aren't (ADR 0001).
   nf-combine will add new gVCFs to VDSs built from the old ones, unless it decides otherwise. The `gvcf`
   analysis meta records the mode, so they can be told apart.
2. **Unvalidated genotyping configuration.** Plain HaplotypeCaller without BQSR, on bwa-mem3 CRAMs, has had no
   truth-set benchmark yet (ADR 0001).
3. **Reading CRAMs directly from GCS.** Each genotyping interval's task reads the CRAM from GCS through
   its index, instead of copying it to local disk. This hasn't been tested on Google Batch, or with CRAMs on
   both references (the DRAGEN-masked FASTA and unmasked GATK hg38). If it fails, fall back to fewer, larger
   genotyping intervals with the CRAM copied to local disk.
