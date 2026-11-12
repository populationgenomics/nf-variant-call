"""
bin/postproc_gvcf.sh over small fixture gVCFs, with the real bcftools and tabix. GATK isn't run: a fake `gatk` on
PATH logs its arguments, and stands in for MergeVcfs (bcftools concat) and ReblockGVCF (a copy).
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

BIN = Path(__file__).parents[1] / 'bin'

pytestmark = pytest.mark.skipif(
    not (shutil.which('bcftools') and shutil.which('tabix')),
    reason='bcftools and tabix not installed',
)

CONTIGS = [f'chr{n}' for n in range(1, 23)] + ['chrX', 'chrY', 'chrM', 'chrUn_KI270302v1']
NOALT = [c for c in CONTIGS if not c.startswith('chrUn')]

FAKE_GATK = """#!/usr/bin/env bash
# logs each call, then does just enough of the tool for the script to carry on
set -euo pipefail
echo "$*" >> "$GATK_LOG"
if [[ $1 == --version ]]; then
    printf 'Using GATK jar /gatk/gatk.jar\\nThe Genome Analysis Toolkit (GATK) v4.6.2.0\\nHTSJDK Version: 4.2.0\\n'
    exit 0
fi
shift 2  # --java-options "..."
tool=$1
shift
inputs=()
while [[ $# -gt 0 ]]; do
    case $1 in
        -I|-V) inputs+=("$2"); shift 2 ;;
        -O) output=$2; shift 2 ;;
        *) shift ;;
    esac
done
if [[ $tool == MergeVcfs ]]; then
    bcftools concat -Oz -o "$output" "${inputs[@]}"
else
    cp "${inputs[0]}" "$output"
fi
tabix -f -p vcf "$output"
"""

SELFSM = '#SEQ_ID\tRG\tCHIP_ID\t#SNPS\t#READS\tAVG_DP\tFREEMIX\nCPG1\tNA\tNA\t9000\t100\t30\t0.0123\n'


def write_gvcf(path: Path, contigs: list[str], sample: str = 'NA12878') -> Path:
    """A bgzipped, indexed gVCF with a reference block and an SNV on each contig, and INFO/DP and INFO/DS."""
    header = [
        '##fileformat=VCFv4.2',
        '##ALT=<ID=NON_REF,Description="Any possible alternative allele">',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">',
        '##INFO=<ID=DS,Number=0,Type=Flag,Description="Downsampled">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        *(f'##contig=<ID={c},length=1000000>' for c in CONTIGS),
        '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample]),
    ]
    records = []
    for contig in contigs:
        records.append(f'{contig}\t100\t.\tA\t<NON_REF>\t.\t.\tEND=199;DP=30\tGT\t0/0')
        records.append(f'{contig}\t200\t.\tC\tT,<NON_REF>\t50\t.\tDP=40;DS\tGT\t0/1')
    plain = path.with_suffix('')
    plain.write_text('\n'.join(header + records) + '\n')
    subprocess.run(['bcftools', 'view', '-Oz', '-o', str(path), str(plain)], check=True)
    subprocess.run(['tabix', '-p', 'vcf', str(path)], check=True)
    return path


@pytest.fixture
def postproc(tmp_path):
    """Run postproc_gvcf.sh in tmp_path/work over the given gVCFs; returns (result, workdir, gatk calls)."""
    fake_bin = tmp_path / 'fake_bin'
    fake_bin.mkdir()
    (fake_bin / 'gatk').write_text(FAKE_GATK)
    (fake_bin / 'gatk').chmod(0o755)
    noalt = tmp_path / 'noalt.bed'
    noalt.write_text(''.join(f'{c}\t0\t1000000\n' for c in NOALT))
    selfsm = tmp_path / 'CPG1.verify-bamid.selfSM'
    selfsm.write_text(SELFSM)
    workdir = tmp_path / 'work'
    workdir.mkdir()
    log = tmp_path / 'gatk.log'
    env = os.environ | {'PATH': f'{fake_bin}:{BIN}:{os.environ["PATH"]}', 'GATK_LOG': str(log)}

    def _run(*gvcfs: Path, bands: str = '13,20,30,40'):
        result = subprocess.run(
            [str(BIN / 'postproc_gvcf.sh'), 'CPG1', 'ref.fasta', str(noalt), str(selfsm), bands, '6',
             *map(str, gvcfs)],
            cwd=workdir, env=env, capture_output=True, text=True, check=False,
        )  # fmt: skip
        calls = log.read_text().splitlines() if log.exists() else []
        return result, workdir, calls

    return _run


def query(gvcf: Path, *args: str) -> str:
    return subprocess.run(['bcftools', *args, str(gvcf)], capture_output=True, text=True, check=True).stdout


def test_one_gvcf_is_postprocessed_without_merging(postproc, tmp_path):
    gvcf = write_gvcf(tmp_path / 'CPG1.0001.g.vcf.gz', CONTIGS)
    result, workdir, calls = postproc(gvcf)
    assert result.returncode == 0, result.stdout + result.stderr
    assert not [c for c in calls if 'MergeVcfs' in c]
    reblock = next(c for c in calls if 'ReblockGVCF' in c)
    assert '--reference ref.fasta -V CPG1-nodp.g.vcf.gz -do-qual-approx --floor-blocks ' in reblock
    assert '-GQB 13 -GQB 20 -GQB 30 -GQB 40 ' in reblock
    assert (workdir / 'CPG1.g.vcf.gz').exists()
    assert (workdir / 'CPG1.g.vcf.gz.tbi').exists()


def test_the_sample_is_renamed_and_info_dp_and_ds_are_dropped(postproc, tmp_path):
    gvcf = write_gvcf(tmp_path / 'CPG1.0001.g.vcf.gz', CONTIGS)
    result, workdir, _calls = postproc(gvcf)
    assert result.returncode == 0, result.stdout + result.stderr
    out = workdir / 'CPG1.g.vcf.gz'
    assert query(out, 'query', '-l') == 'CPG1\n'
    header = query(out, 'view', '-h')
    assert '##INFO=<ID=DP,' not in header
    assert '##INFO=<ID=DS,' not in header
    assert set(query(out, 'query', '-f', '%INFO\n').split()) == {'END=199', '.'}


def test_records_outside_the_noalt_regions_are_dropped(postproc, tmp_path):
    gvcf = write_gvcf(tmp_path / 'CPG1.0001.g.vcf.gz', CONTIGS)
    result, workdir, _calls = postproc(gvcf)
    assert result.returncode == 0, result.stdout + result.stderr
    contigs = query(workdir / 'CPG1.g.vcf.gz', 'query', '-f', '%CHROM\n').split()
    assert sorted(set(contigs)) == sorted(NOALT)


def test_interval_gvcfs_are_merged_in_the_order_given(postproc, tmp_path):
    first = write_gvcf(tmp_path / 'CPG1.0001.g.vcf.gz', CONTIGS[:12])
    second = write_gvcf(tmp_path / 'CPG1.0002.g.vcf.gz', CONTIGS[12:])
    result, workdir, calls = postproc(first, second)
    assert result.returncode == 0, result.stdout + result.stderr
    merge = next(c for c in calls if 'MergeVcfs' in c)
    assert f'-I {first} -I {second} -O merged.g.vcf.gz' in merge
    assert '-V CPG1-nodp.g.vcf.gz' in next(c for c in calls if 'ReblockGVCF' in c)
    contigs = query(workdir / 'CPG1.g.vcf.gz', 'query', '-f', '%CHROM\n').split()
    assert sorted(set(contigs)) == sorted(NOALT)


def test_a_missing_autosome_fails_the_contig_check(postproc, tmp_path):
    gvcf = write_gvcf(tmp_path / 'CPG1.0001.g.vcf.gz', [c for c in CONTIGS if c != 'chr21'])
    result, workdir, calls = postproc(gvcf)
    assert result.returncode != 0
    assert 'ERROR: gVCF missing expected contig chr21 for CPG1' in result.stdout
    assert 'FAIL: contig check failed for CPG1' in result.stdout
    assert not [c for c in calls if 'ReblockGVCF' in c]
    assert not (workdir / 'CPG1.g.vcf.gz').exists()


def test_gatk_version_and_contamination_are_recorded(postproc, tmp_path):
    gvcf = write_gvcf(tmp_path / 'CPG1.0001.g.vcf.gz', CONTIGS)
    result, workdir, _calls = postproc(gvcf)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (workdir / 'gatk_version.txt').read_text() == '4.6.2.0\n'
    assert (workdir / 'contamination.txt').read_text() == '0.0123\n'


def test_custom_gq_bands_are_passed_to_reblock(postproc, tmp_path):
    gvcf = write_gvcf(tmp_path / 'CPG1.0001.g.vcf.gz', CONTIGS)
    result, _workdir, calls = postproc(gvcf, bands='10,20')
    assert result.returncode == 0, result.stdout + result.stderr
    assert '--floor-blocks -GQB 10 -GQB 20 -O ' in next(c for c in calls if 'ReblockGVCF' in c)
