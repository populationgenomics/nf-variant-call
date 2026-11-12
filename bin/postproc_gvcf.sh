#!/usr/bin/env bash
# One SG's interval gVCFs -> its single-sample, reblocked gVCF, <sg>.g.vcf.gz and .tbi, in the working directory.
#   1. More than one interval gVCF is merged with gatk MergeVcfs, in the order given (genomic order).
#   2. cpg-flow's postprocessing, unchanged (genotype.py, postprocess_gvcf): a contig check, INFO/DP dropped,
#      ReblockGVCF, records outside the no-alt regions dropped, INFO/DS dropped, and the sample renamed to the SG.
#   3. The GATK version goes to gatk_version.txt and FREEMIX (freemix.sh) to contamination.txt, for the analysis
#      meta.
#   usage: postproc_gvcf.sh <sg_id> <ref.fasta> <noalt.bed> <selfSM> <GQ bands, comma-separated> <heap GB> <gvcf>...
set -euo pipefail

if [[ $# -lt 7 ]]; then
    echo 'usage: postproc_gvcf.sh <sg_id> <ref.fasta> <noalt.bed> <selfSM> <GQ bands> <heap GB> <gvcf>...' >&2
    exit 1
fi
SG_ID=$1
REF=$2
NOALT_BED=$3
SELFSM=$4
IFS=, read -ra GQ_BANDS <<< "$5"
HEAP_GB=$6
shift 6

# recorded first, so a bad .selfSM or image fails before the long steps
gatk --version | sed -n 's/^.*(GATK) v//p' > gatk_version.txt
if [[ ! -s gatk_version.txt ]]; then
    echo 'postproc_gvcf.sh: no GATK version in gatk --version' >&2
    exit 1
fi
freemix.sh "${SELFSM}" > contamination.txt

if [[ $# -gt 1 ]]; then
    gatk --java-options "-Xms${HEAP_GB}g -Xmx${HEAP_GB}g" \
        MergeVcfs \
        $(printf -- '-I %s ' "$@") \
        -O merged.g.vcf.gz
    GVCF=merged.g.vcf.gz
else
    GVCF=$1
fi

GVCF_NODP=${SG_ID}-nodp.g.vcf.gz
REBLOCKED=${SG_ID}-reblocked.g.vcf.gz
OUTPUT=${SG_ID}.g.vcf.gz

# From here on, cpg-flow's postprocessing, with its temporary files in the working directory.

# Reindexing just to make sure the index is not corrupted
bcftools index -f --tbi $GVCF

EXPECTED_CONTIGS="$(printf 'chr%s ' {1..22})chrX"
FAIL_FILE=contig_check_failures.txt
rm -f $FAIL_FILE
tabix --list-chroms $GVCF > actual_contigs.txt
for contig in $EXPECTED_CONTIGS; do
    if ! grep -qw "$contig" actual_contigs.txt; then
        echo "ERROR: gVCF missing expected contig $contig for ${SG_ID}" | tee -a $FAIL_FILE
    fi
done
if [ -f $FAIL_FILE ]; then
    echo "FAIL: contig check failed for ${SG_ID}"
    exit 1
fi

# Remove INFO/DP field, which contradicts the FORMAT/DP, in the way that
# it has _all_ reads, not just variant-calling-usable reads. If we keep INFO/DP,
# ReblockGVCF would prioritize it over FORMAT/DP and change FORMAT/DP to INFO/DP
# in the resulting merged blocks. It would pick the highest INFO/DP when merging
# multiple blocks, so a variant in a small homopolymer region (surrounded by
# long DP=0 areas), that attracted piles of low-MQ reads with INFO/DP=1000
# will translate into a long GQ<20 block with the same FORMAT/DP=1000,
# which is wrong, because most of this block has no reads.
bcftools view $GVCF \
| bcftools annotate -x INFO/DP \
| bcftools view -Oz -o $GVCF_NODP
tabix -p vcf $GVCF_NODP

gatk --java-options "-Xms${HEAP_GB}g -Xmx${HEAP_GB}g" \
ReblockGVCF \
--reference ${REF} \
-V $GVCF_NODP \
-do-qual-approx \
--floor-blocks $(printf -- '-GQB %s ' "${GQ_BANDS[@]}") \
-O $REBLOCKED \
--create-output-variant-index true

EXISTING_SN=$(bcftools query -l $GVCF)

bcftools view $REBLOCKED -T ${NOALT_BED} \
| bcftools annotate -x INFO/DS \
| bcftools reheader -s <(echo "$EXISTING_SN ${SG_ID}") \
| bcftools view -Oz -o ${OUTPUT}

tabix -p vcf ${OUTPUT}
