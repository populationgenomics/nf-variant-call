"""bin/bcftools_stats.sh over a fixture gVCF, with the real bcftools: only reference blocks are filtered out."""

import shutil
import subprocess
from pathlib import Path

import pytest

BCFTOOLS_STATS = Path(__file__).parents[1] / 'bin' / 'bcftools_stats.sh'

pytestmark = pytest.mark.skipif(shutil.which('bcftools') is None, reason='bcftools not installed')

HEADER = [
    '##fileformat=VCFv4.2',
    '##ALT=<ID=NON_REF,Description="Any possible alternative allele">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##contig=<ID=chr1,length=1000>',
    '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'CPG1']),
]
# (POS, REF, ALT, INFO, GT): two reference blocks around a SNV, a multi-allelic site and an indel
RECORDS = [
    (1, 'A', '<NON_REF>', 'END=99', '0/0'),
    (100, 'C', 'T,<NON_REF>', '.', '0/1'),
    (200, 'G', 'C,T,<NON_REF>', '.', '1/2'),
    (300, 'AT', 'A,<NON_REF>', '.', '0/1'),
    (302, 'T', '<NON_REF>', 'END=1000', '0/0'),
]


def write_fixture(tmp_path: Path) -> tuple[Path, Path]:
    """The gVCF, bgzipped and indexed, and a reference whose bases agree with its REF alleles."""
    bases = ['A'] * 1000
    for pos, ref, *_ in RECORDS:
        bases[pos - 1 : pos - 1 + len(ref)] = ref
    fasta = tmp_path / 'ref.fasta'
    fasta.write_text('>chr1\n' + ''.join(bases) + '\n')

    plain = tmp_path / 'CPG1.g.vcf'
    lines = [f'chr1\t{pos}\t.\t{ref}\t{alt}\t50\t.\t{info}\tGT\t{gt}' for pos, ref, alt, info, gt in RECORDS]
    plain.write_text('\n'.join(HEADER + lines) + '\n')
    gvcf = tmp_path / 'CPG1.g.vcf.gz'
    subprocess.run(['bcftools', 'view', '-Oz', '-o', str(gvcf), '--write-index=tbi', str(plain)], check=True)
    return gvcf, fasta


def summary(stats: Path) -> dict[str, int]:
    """The SN (summary numbers) section of bcftools stats: 'number of records' and so on, to their counts."""
    rows = [line.split('\t') for line in stats.read_text().splitlines() if line.startswith('SN\t')]
    return {key.removeprefix('number of ').rstrip(':'): int(value) for _sn, _id, key, value in rows}


def test_only_reference_blocks_are_dropped(tmp_path):
    gvcf, fasta = write_fixture(tmp_path)
    stats = tmp_path / 'CPG1.bcftools-stats'
    subprocess.run([str(BCFTOOLS_STATS), str(gvcf), str(fasta), str(stats)], check=True)
    counts = summary(stats)
    # exactly the three variant records: the SNV, the multi-allelic site and the indel
    assert counts['records'] == 3
    assert counts['SNPs'] == 2
    assert counts['indels'] == 1
    # <NON_REF> is an allele to bcftools, so every variant site in a gVCF counts as multi-allelic
    assert counts['multiallelic sites'] == 3


def test_a_missing_gvcf_fails(tmp_path):
    _gvcf, fasta = write_fixture(tmp_path)
    stats = tmp_path / 'CPG1.bcftools-stats'
    result = subprocess.run(
        [str(BCFTOOLS_STATS), str(tmp_path / 'absent.g.vcf.gz'), str(fasta), str(stats)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
