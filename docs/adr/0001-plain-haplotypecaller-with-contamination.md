# Genotype with plain HaplotypeCaller, no DRAGEN mode or BQSR, with contamination from CRAM QC

cpg-flow genotyped with `HaplotypeCaller --dragen-mode`, but without the DRAGstr calibration that full
DRAGEN-GATK mode expects. That mode was designed for DRAGMAP alignments, and CRAMs now come from bwa-mem3
(nf-align). We drop `--dragen-mode` and keep every other cpg-flow HaplotypeCaller flag. We don't add BQSR
either: GATK4 can't apply it on the fly, so it would mean writing and storing a second, recalibrated CRAM
per SG. Instead, each SG's VerifyBamID FREEMIX, from its `cram_qc` analysis, is passed as
`--contamination-fraction-to-filter`. Genotyping therefore always follows alignment and CRAM QC.

## Considered Options

- Keep `--dragen-mode` exactly as cpg-flow ran it. This is the only setup with a concordance check (about 98%
  against the original gVCFs), but it's DRAGEN mode, half-configured, on a non-DRAGEN aligner.
- Full DRAGEN-GATK, with `ComposeSTRTableFile` and `CalibrateDragstrModel` per SG.
- Broad best practice with BQSR. Rejected for the cost of a recalibrated CRAM.
- Contamination only when a `cram_qc` analysis happens to exist. Rejected: the same CRAM could then give
  different gVCFs depending on the order things ran in.

## Consequences

- This configuration is unvalidated. Before production use, a GIAB sample is benchmarked against truth with
  hap.py (bwa-mem3 plus this genotyping), and compared with the DRAGEN-mode baseline.
- Every existing gVCF and VDS is DRAGEN-mode. Each new `gvcf` analysis records its caller version, that
  DRAGEN mode is off, and the contamination fraction used, so VDSs that mix the two can be split later. Whether a VDS may mix them is
  nf-combine's decision.
