"""bin/build_samplesheet.py: Metamist responses -> sample sheet rows. The query is validated against the
schema bundled with the metamist client; responses are hand-built dicts in the shape that query returns."""

import csv
import logging
import re
from argparse import Namespace

import build_samplesheet
import pytest
from build_samplesheet import SamplesheetError, build_rows, fetch_cohort_sgs, fetch_sgs, row_for, summary, write_sheet
from conftest import (
    OTHER,
    PROJECT,
    cohort,
    cram_analysis,
    cram_qc_analysis,
    gvcf_analysis,
    gvcf_qc_analysis,
    sg,
)
from metamist.graphql import gql, validate

# the references of a CRAM or gVCF whose analysis doesn't record one (--aligned-ref-fasta, --genotyped-ref-fasta)
ALIGNED = 'gs://refs/aligned_masked.fasta'
GENOTYPED = 'gs://refs/genotyped_masked.fasta'
CRAM = 'gs://b/CPG1.cram'
SELFSM = 'gs://b/CPG1.verify-bamid.selfSM'
GVCF = 'gs://b/CPG1.g.vcf.gz'


def row(*analyses: dict) -> dict:
    return row_for(sg('CPG1', list(analyses)), ALIGNED, GENOTYPED)


def ready(cram: str = CRAM, **meta: str) -> list[dict]:
    """The analyses of an SG that's ready to genotype: a CRAM, and its CRAM QC with a VerifyBamID file."""
    return [cram_analysis(1, cram, **meta), cram_qc_analysis(2, cram, verifybamid=SELFSM)]


def rows(sgs: list[dict], access_level: str = 'standard', sequencing_type: str = 'genome') -> list[dict]:
    return build_rows(sgs, access_level, sequencing_type, ALIGNED, GENOTYPED)


@pytest.mark.parametrize('query', [build_samplesheet.COHORTS_QUERY, build_samplesheet.SGS_QUERY])
def test_queries_match_metamist_schema(query):
    assert validate(gql(query), use_local_schema=True)


def test_rows_from_the_cohort_fixture(cohort_sgs):
    sheet = rows(cohort_sgs)
    assert [(r['sg_id'], r['dataset'], bool(r.get('cram')), bool(r.get('gvcf'))) for r in sheet] == [
        ('CPG1', PROJECT, True, False),
        ('CPG2', PROJECT, True, False),
        ('CPG3', OTHER, True, False),
        ('CPG4', OTHER, False, True),
        ('CPG5', PROJECT, False, True),
    ]
    assert sheet[0]['cram_ref'].endswith('unmasked.fasta')
    assert sheet[1]['cram_ref'] == ALIGNED
    assert sheet[0]['selfsm'].endswith('CPG1.verify-bamid.selfSM')
    assert (sheet[3]['gvcf_qc'], sheet[4]['gvcf_qc']) == ('2', '')


# to genotype


def test_to_genotype_row_carries_the_cram_its_reference_and_selfsm():
    assert row(*ready(ref_fasta='gs://refs/unmasked.fasta')) == {
        'sg_id': 'CPG1',
        'cram': CRAM,
        'cram_ref': 'gs://refs/unmasked.fasta',
        'selfsm': SELFSM,
    }


def test_to_genotype_cram_without_a_recorded_reference_gets_aligned_ref_fasta():
    assert row(*ready())['cram_ref'] == ALIGNED
    # meta can be null too
    analyses = [{'id': 1, 'type': 'cram', 'outputs': {'path': CRAM}, 'meta': None}, ready()[1]]
    assert row(*analyses)['cram_ref'] == ALIGNED


def test_latest_cram_analysis_wins():
    new = 'gs://b/new.cram'
    analyses = [cram_analysis(9, new), cram_analysis(3, CRAM), cram_qc_analysis(10, new, verifybamid=SELFSM)]
    assert row(*analyses)['cram'] == new


def test_to_genotype_without_a_cram_fails():
    with pytest.raises(SamplesheetError, match='CPG1: to genotype, but has no CRAM'):
        row()
    # an analysis without an output doesn't count
    with pytest.raises(SamplesheetError, match='CPG1: to genotype, but has no CRAM'):
        row(cram_analysis(1, None))


def test_to_genotype_without_cram_qc_fails():
    with pytest.raises(SamplesheetError, match=rf'CPG1: to genotype, but its CRAM {re.escape(CRAM)} has no CRAM QC'):
        row(cram_analysis(1, CRAM))


def test_cram_qc_of_a_different_cram_does_not_count():
    analyses = [cram_analysis(9, 'gs://b/new.cram'), cram_analysis(3, CRAM), cram_qc_analysis(4, CRAM, SELFSM)]
    with pytest.raises(SamplesheetError, match=r'its CRAM gs://b/new\.cram has no CRAM QC'):
        row(*analyses)


def test_cram_qc_without_a_verifybamid_file_fails():
    with pytest.raises(SamplesheetError, match=r'CRAM QC of gs://b/CPG1\.cram \(analysis 2\) has no verifybamid'):
        row(cram_analysis(1, CRAM), cram_qc_analysis(2, CRAM))


def test_latest_cram_qc_of_the_cram_gives_the_selfsm():
    analyses = [
        cram_analysis(1, CRAM),
        cram_qc_analysis(2, CRAM, verifybamid='gs://b/old.selfSM'),
        cram_qc_analysis(5, CRAM, verifybamid=SELFSM),
    ]
    assert row(*analyses)['selfsm'] == SELFSM


# genotyped


def test_genotyped_row_carries_the_gvcf_and_needs_no_cram():
    assert row(gvcf_analysis(1, GVCF)) == {'sg_id': 'CPG1', 'gvcf': GVCF, 'gvcf_ref': GENOTYPED, 'gvcf_qc': ''}


def test_genotyped_row_carries_the_reference_its_analysis_records():
    assert row(gvcf_analysis(1, GVCF, ref_fasta='gs://refs/unmasked.fasta'))['gvcf_ref'] == 'gs://refs/unmasked.fasta'


def test_genotyped_sg_is_never_genotyped_again():
    # even with a newer CRAM than its gVCF came from (stale gVCFs are an accepted gap)
    result = row(gvcf_analysis(1, GVCF), *ready(cram='gs://b/newer.cram'))
    assert result['gvcf'] == GVCF
    assert not result.get('cram')


def test_latest_gvcf_analysis_wins():
    new = 'gs://b/new.g.vcf.gz'
    assert row(gvcf_analysis(3, GVCF), gvcf_analysis(9, new))['gvcf'] == new


def test_genotyped_row_records_its_gvcfs_own_gvcf_qc():
    assert row(gvcf_analysis(1, GVCF), gvcf_qc_analysis(2, GVCF), gvcf_qc_analysis(4, GVCF))['gvcf_qc'] == '4'


def test_gvcf_qc_of_a_different_gvcf_does_not_count():
    # two gvcf_qc analyses, each on its own gVCF: only the current gVCF's counts
    new = 'gs://b/new.g.vcf.gz'
    analyses = [gvcf_analysis(9, new), gvcf_analysis(3, GVCF), gvcf_qc_analysis(4, GVCF), gvcf_qc_analysis(10, new)]
    assert row(*analyses)['gvcf_qc'] == '10'
    assert row(*analyses[:3])['gvcf_qc'] == ''


def test_back_filled_gvcf_qc_without_bcftools_stats_counts():
    assert row(gvcf_analysis(1, GVCF), gvcf_qc_analysis(2, GVCF, bcftools_stats=False))['gvcf_qc'] == '2'


def test_gvcf_qc_without_a_gvcf_analysis_is_not_genotyped():
    assert row(gvcf_qc_analysis(5, GVCF), *ready())['cram'] == CRAM


# scope, datasets and errors


def test_sgs_are_deduplicated_and_filtered_to_the_sequencing_type(caplog):
    # CPG1 is in two cohorts
    sgs = [
        sg('CPG1', ready()),
        sg('CPG2', ready(), type='exome'),
        sg('CPG1', ready()),
        sg('CPG3', ready(), technology='long-read'),
    ]
    assert [r['sg_id'] for r in rows(sgs)] == ['CPG1']
    assert [r['sg_id'] for r in rows(sgs, sequencing_type='exome')] == ['CPG2']
    assert 'CPG3: skipping long-read genome' in caplog.text


def test_sgs_from_several_projects_keep_their_own_dataset():
    # one cohort, two projects: the SG's own project decides, not the cohort's
    sgs = [sg('CPG1', ready()), sg('CPG2', [gvcf_analysis(1, GVCF)], project=OTHER)]
    assert [(r['sg_id'], r['dataset']) for r in rows(sgs, 'full')] == [('CPG1', PROJECT), ('CPG2', OTHER)]


def test_test_projects_give_their_dataset_without_the_suffix():
    sgs = [sg('CPG1', ready(), project='ds-test'), sg('CPG2', ready(), project='other-test')]
    assert [r['dataset'] for r in rows(sgs, 'test')] == [PROJECT, OTHER]


@pytest.mark.parametrize('access_level', ['full', 'standard'])
def test_test_projects_fail_a_production_run(access_level):
    sgs = [sg('CPG1', ready()), sg('CPG2', ready(), project='ds-test'), sg('CPG3', ready(), project='other-test')]
    with pytest.raises(SamplesheetError) as error:
        rows(sgs, access_level)
    assert str(error.value).splitlines() == [
        f'CPG2: project ds-test does not match access level {access_level}',
        f'CPG3: project other-test does not match access level {access_level}',
    ]


def test_production_projects_fail_a_test_run():
    sgs = [sg('CPG1', ready(), project='ds-test'), sg('CPG2', ready()), sg('CPG3', ready(), project=OTHER)]
    with pytest.raises(SamplesheetError) as error:
        rows(sgs, 'test')
    assert str(error.value).splitlines() == [
        'CPG2: project ds does not match access level test',
        'CPG3: project other does not match access level test',
    ]


def test_no_in_scope_sgs_fails():
    with pytest.raises(SamplesheetError, match='No short-read exome sequencing groups to genotype'):
        rows([sg('CPG1', ready())], sequencing_type='exome')


def test_bad_sgs_are_reported_together():
    sgs = [
        sg('CPG1'),
        sg('CPG2', [cram_analysis(1, CRAM)]),
        sg('CPG3', [gvcf_analysis(1, GVCF)]),  # fine
        sg('CPG4', [cram_analysis(1, CRAM), cram_qc_analysis(2, CRAM)]),
    ]
    with pytest.raises(SamplesheetError) as error:
        rows(sgs)
    assert [line.split(':')[0] for line in str(error.value).splitlines()] == ['CPG1', 'CPG2', 'CPG4']


def test_summary_counts_each_kind_of_sg(cohort_sgs):
    assert summary(rows(cohort_sgs)) == (
        '5 SGs from 2 dataset(s): 3 to genotype, 2 already genotyped (1 of them without gVCF QC)'
    )


def test_main_logs_the_counts(monkeypatch, tmp_path, cohort_sgs, caplog):
    monkeypatch.setattr(build_samplesheet, 'fetch_cohort_sgs', lambda _ids: cohort_sgs)
    args = Namespace(
        cohorts=['COH1'],
        sgids=None,
        access_level='standard',
        sequencing_type='genome',
        aligned_ref_fasta=ALIGNED,
        genotyped_ref_fasta=GENOTYPED,
        output=str(tmp_path / 'sheet.tsv'),
    )
    with caplog.at_level(logging.INFO):
        build_samplesheet.main(args)
    assert '3 to genotype, 2 already genotyped (1 of them without gVCF QC)' in caplog.text


# fetching


def test_fetch_cohorts_reports_unusable_cohorts(monkeypatch):
    calls = []

    def fake_query(_query, variables):
        calls.append(variables)
        return {
            'cohorts': [
                cohort('COH1', [sg('CPG1')]),
                cohort('COH2', [sg('CPG2')], status='archived'),
                cohort('COH3', []),
            ],
        }

    monkeypatch.setattr(build_samplesheet, 'query', fake_query)
    # COH4 may not exist, or may be in a project this run's credentials can't read
    with pytest.raises(
        SamplesheetError,
        match=r"missing or not visible to this run: \['COH4'\], inactive: \['COH2'\], empty: \['COH3'\]",
    ):
        fetch_cohort_sgs(['COH1', 'COH2', 'COH3', 'COH4'])
    assert calls == [{'cohorts': ['COH1', 'COH2', 'COH3', 'COH4']}]


def test_fetch_cohort_sgs_flattens_cohorts(monkeypatch):
    response = {'cohorts': [cohort('COH1', [sg('CPG1'), sg('CPG2', project=OTHER)]), cohort('COH2', [sg('CPG1')])]}
    monkeypatch.setattr(build_samplesheet, 'query', lambda _query, **_variables: response)
    assert [s['id'] for s in fetch_cohort_sgs(['COH1', 'COH2'])] == ['CPG1', 'CPG2', 'CPG1']


def test_fetch_sgs_reports_missing_sgs(monkeypatch):
    calls = []

    def fake_query(_query, variables):
        calls.append(variables)
        return {'sequencingGroups': [sg('CPG1'), sg('CPG2', project=OTHER)]}

    monkeypatch.setattr(build_samplesheet, 'query', fake_query)
    assert [s['id'] for s in fetch_sgs(['CPG1', 'CPG2'], 'genome')] == ['CPG1', 'CPG2']
    # CPG3 may not exist, or may be in a project this run's credentials can't read
    with pytest.raises(SamplesheetError, match=r"missing or not visible to this run: \['CPG3'\]"):
        fetch_sgs(['CPG1', 'CPG2', 'CPG3'], 'genome')
    assert calls[-1] == {'sgids': ['CPG1', 'CPG2', 'CPG3']}


def test_fetch_sgs_fails_on_sgs_of_the_wrong_type(monkeypatch):
    # named SGs fail rather than being skipped, unlike SGs found through a cohort
    sgs = [sg('CPG1'), sg('CPG2', type='exome'), sg('CPG3', technology='long-read')]
    response = {'sequencingGroups': sgs}
    monkeypatch.setattr(build_samplesheet, 'query', lambda _query, **_variables: response)
    with pytest.raises(
        SamplesheetError,
        match=r"not short-read genome: \['CPG2 \(short-read exome\)', 'CPG3 \(long-read genome\)'\]",
    ):
        fetch_sgs(['CPG1', 'CPG2', 'CPG3'], 'genome')


def test_write_sheet_has_every_column(tmp_path, cohort_sgs):
    path = tmp_path / 'sheet.tsv'
    write_sheet(rows(cohort_sgs), str(path))
    with path.open() as handle:
        written = list(csv.DictReader(handle, delimiter='\t'))
    assert list(written[0]) == build_samplesheet.COLUMNS
    assert written[0]['gvcf'] == ''
    assert written[4]['cram'] == ''
