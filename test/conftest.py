"""Metamist-shaped fixtures shared by the sample sheet and stub-run tests. Metamist is never contacted."""

from pathlib import Path

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
