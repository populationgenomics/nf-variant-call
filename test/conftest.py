"""Metamist-shaped fixtures shared by the sample sheet and stub-run tests. Metamist is never contacted."""

from pathlib import Path

import pytest

PROJECT = 'ds'
OTHER = 'other'  # a second dataset, so runs span two


def touch(path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()
    return str(path)


def sg(
    sg_id: str,
    analyses: list[dict] | None = None,
    project: str = PROJECT,
    **fields: object,
) -> dict:
    """A sequencing group, shaped like the SgFields fragment of the sample sheet queries."""
    return {
        'id': sg_id,
        'type': 'genome',
        'technology': 'short-read',
        'sample': {'project': {'name': project}},
        'analyses': analyses or [],
    } | fields


def cohort(cohort_id: str, sgs: list[dict], status: str = 'active') -> dict:
    return {'id': cohort_id, 'status': status, 'sequencingGroups': sgs}


def outputs(path: str | None, **secondary: str) -> dict | None:
    """An analysis's outputs, as cpg_utils.metamist_registration records them: --secondary <key>=<path>."""
    if not path:
        return None
    return {'path': path, 'secondary_files': {key: {'path': value} for key, value in secondary.items()}}


def cram_analysis(analysis_id: int, cram: str | None, **meta: str) -> dict:
    """A CRAM analysis as nf-align registers it: the index as a secondary file. meta may record ref_fasta."""
    return {'id': analysis_id, 'type': 'cram', 'outputs': outputs(cram, index=f'{cram}.crai'), 'meta': meta}


def cram_qc_analysis(analysis_id: int, cram: str, verifybamid: str | None = None) -> dict:
    """A CRAM QC analysis as nf-align registers it: the CRAM is its primary file, VerifyBamID's .selfSM a secondary."""
    secondary = {'verifybamid': verifybamid} if verifybamid else {}
    return {'id': analysis_id, 'type': 'cram_qc', 'outputs': outputs(cram, **secondary), 'meta': {}}


def gvcf_analysis(analysis_id: int, gvcf: str | None, **meta: str) -> dict:
    """A gVCF analysis: the .tbi as a secondary file. cpg-flow's record no ref_fasta, this pipeline's do."""
    return {'id': analysis_id, 'type': 'gvcf', 'outputs': outputs(gvcf, index=f'{gvcf}.tbi'), 'meta': meta}


def gvcf_qc_analysis(analysis_id: int, gvcf: str, bcftools_stats: bool = True) -> dict:
    """
    A gVCF QC analysis: the gVCF is its primary file, the metrics its secondaries. A back-filled one, of cpg-flow's
    QC files, has no bcftools_stats.
    """
    stem = gvcf.removesuffix('.g.vcf.gz')
    secondary = {
        'variant_calling_summary_metrics': f'{stem}.variant_calling_summary_metrics',
        'variant_calling_detail_metrics': f'{stem}.variant_calling_detail_metrics',
    } | ({'bcftools_stats': f'{stem}.bcftools-stats'} if bcftools_stats else {})
    return {'id': analysis_id, 'type': 'gvcf_qc', 'outputs': outputs(gvcf, **secondary), 'meta': {}}


def to_genotype(sg_id: str, data: Path, project: str = PROJECT, analysis_id: int = 1, **meta: str) -> dict:
    """An SG with a CRAM and its CRAM QC, as local placeholder files so a stub run can stage them."""
    cram = touch(data / sg_id / f'{sg_id}.cram')
    touch(data / sg_id / f'{sg_id}.cram.crai')
    selfsm = touch(data / sg_id / f'{sg_id}.verify-bamid.selfSM')
    analyses = [cram_analysis(analysis_id, cram, **meta), cram_qc_analysis(analysis_id + 1, cram, verifybamid=selfsm)]
    return sg(sg_id, analyses, project=project)


def genotyped(sg_id: str, data: Path, project: str = PROJECT, qc: bool = False, analysis_id: int = 1) -> dict:
    """An SG with a gVCF that cpg-flow made (so its analysis records no reference), with gVCF QC if qc."""
    gvcf = touch(data / sg_id / f'{sg_id}.g.vcf.gz')
    touch(data / sg_id / f'{sg_id}.g.vcf.gz.tbi')
    analyses = [gvcf_analysis(analysis_id, gvcf)] + ([gvcf_qc_analysis(analysis_id + 1, gvcf)] if qc else [])
    return sg(sg_id, analyses, project=project)


@pytest.fixture
def cohort_sgs(tmp_path: Path) -> list[dict]:
    """
    SGs across two datasets. CPG1-CPG3 are to genotype: CPG1's CRAM analysis records its reference, CPG2's and
    CPG3's (older CRAMs) don't. CPG4 and CPG5 are genotyped, and only CPG4's gVCF has gVCF QC.
    """
    data = tmp_path / 'data'
    return [
        to_genotype('CPG1', data, ref_fasta=str(tmp_path / 'refs' / 'unmasked.fasta')),
        to_genotype('CPG2', data),
        to_genotype('CPG3', data, project=OTHER),
        genotyped('CPG4', data, project=OTHER, qc=True),
        genotyped('CPG5', data),
    ]
