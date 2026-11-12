#!/usr/bin/env bash
# bcftools stats of a gVCF's variant sites, to <output>. A record is dropped only when <NON_REF> is its only ALT
# allele: a reference block. A variant's ALT is like C,<NON_REF>, and bcftools' ALT="<NON_REF>" matches any of a
# record's alleles, so without N_ALT=1 the filter would drop every variant too.
#   usage: bcftools_stats.sh <gvcf> <ref.fasta> <output>
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo 'usage: bcftools_stats.sh <gvcf> <ref.fasta> <output>' >&2
    exit 1
fi
GVCF=$1
REF=$2
OUTPUT=$3

bcftools view -Ou -e 'N_ALT=1 && ALT="<NON_REF>"' "${GVCF}" \
    | bcftools stats -F "${REF}" - > "${OUTPUT}"
