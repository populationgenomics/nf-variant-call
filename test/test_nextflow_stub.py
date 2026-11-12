"""Offline `nextflow run -stub` checks of main.nf's wiring: params, the run-records check, the published sample
sheet, row parsing and the pre-flight check. The SGs span two datasets (conftest.cohort_sgs): CPG1-CPG3 are to
genotype, CPG4 and CPG5 are genotyped, and only CPG5 has no gVCF QC. Placeholder files stand in for data, and the
tests are skipped when nextflow isn't on PATH."""

import csv
import shutil
import subprocess
from pathlib import Path

import pytest
from build_samplesheet import build_rows, write_sheet
from conftest import OTHER, PROJECT, touch

pytestmark = pytest.mark.skipif(shutil.which('nextflow') is None, reason='nextflow not installed')

MAIN_NF = Path(__file__).parents[1] / 'main.nf'

# every reference param, each passed as a local placeholder so a stub run never reaches GCS
REFERENCES = {
    'aligned_ref_fasta': 'masked.fasta',
    'genotyped_ref_fasta': 'masked.fasta',
    'calling_regions_genome': 'genome_calling.interval_list',
    'calling_regions_exome': 'exome_calling.interval_list',
    'calling_exclusions': 'exclusions.interval_list',
    'noalt_bed': 'noalt.bed.gz',
    'dbsnp_vcf': 'dbsnp.vcf',
    'genome_evaluation_regions': 'genome_evaluation.bed',
    'exome_evaluation_regions': 'exome_evaluation.bed',
}
# files that sit beside the references, and the reference CPG1's CRAM analysis records (conftest.cohort_sgs)
COMPANIONS = [
    'masked.fasta.fai',
    'masked.dict',
    'dbsnp.vcf.idx',
    'unmasked.fasta',
    'unmasked.fasta.fai',
    'unmasked.dict',
]


def read_tsv(path: Path, sep: str = '\t') -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=sep))


# the run records are under the analysis dataset's <output_root>/<ds>
def run_dir(out: Path, ar_guid: str = 'local-test') -> Path:
    return out / PROJECT / 'genotyping_samplesheets' / ar_guid


# each SG's outputs are under its own dataset's <output_root>/<ds>, with exomes under exome/
def gvcf_dir(out: Path, ds: str, prefix: str = '') -> Path:
    return out / ds / prefix / 'gvcf'


def gvcf_qc_dir(out: Path, ds: str, sg_id: str, prefix: str = '') -> Path:
    return out / ds / prefix / 'gvcf_qc' / sg_id


GVCF_QC_SUFFIXES = ['variant_calling_summary_metrics', 'variant_calling_detail_metrics', 'bcftools-stats']
QC_ONLY = 'CPG5'  # with --qc_genotyped: genotyped, without gVCF QC


def completed(tasks: list[dict]) -> set[str]:
    """Names of the tasks that ran, other than building the sample sheet, which every new run does first."""
    return {t['name'] for t in tasks if t['status'] == 'COMPLETED' and t['name'] != 'BUILD_SAMPLESHEET (standard)'}


def edit_sheet(tmp_path: Path, sg_id: str, **values: str) -> None:
    """Set columns of one SG's row in the sheet BUILD_SAMPLESHEET's stub hands over."""
    sheet = tmp_path / 'sheet.tsv'
    rows = read_tsv(sheet)
    for row in rows:
        if row['sg_id'] == sg_id:
            row.update(values)
    write_sheet(rows, str(sheet))


@pytest.fixture
def run(tmp_path, cohort_sgs):
    """
    Run main.nf in stub mode over a sheet built from conftest.cohort_sgs, across two datasets, with PROJECT as the
    analysis dataset; returns (result, output_root, trace rows). BUILD_SAMPLESHEET's stub hands over sheet.tsv,
    which a test may edit before its first run.
    """
    refs = tmp_path / 'refs'
    sheet = tmp_path / 'sheet.tsv'
    masked = str(refs / REFERENCES['aligned_ref_fasta'])
    write_sheet(build_rows(cohort_sgs, 'standard', 'genome', masked, masked), str(sheet))
    sheet_config = tmp_path / 'sheet.config'
    sheet_config.write_text(f"process {{ withName: 'BUILD_SAMPLESHEET' {{ ext.stub_samplesheet = '{sheet}' }} }}\n")
    for name in [*set(REFERENCES.values()), *COMPANIONS]:
        touch(refs / name)

    def _run(
        *args: str,
        sequencing_type: str = 'genome',
        inputs: tuple = ('--cohorts', 'COH1,COH2'),
        ar_guid: str = 'local-test',
        omit: tuple = (),  # reference params left off the command line
    ):
        out = tmp_path / 'results'
        trace = tmp_path / 'trace.txt'
        trace.unlink(missing_ok=True)
        result = subprocess.run(
            [
                'nextflow', 'run', str(MAIN_NF), '-stub', '-profile', 'test',
                '-with-trace', str(trace),
                '-c', str(sheet_config),
                *inputs,
                '--output_root', str(out),
                '--dataset', PROJECT,
                '--arGuid', ar_guid,
                '--sequencing_type', sequencing_type,
                *(a for name, ref in REFERENCES.items() if name not in omit for a in (f'--{name}', str(refs / ref))),
                *args,
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
            check=False,
        )  # fmt: skip
        tasks = read_tsv(trace) if trace.exists() else []
        return result, out, tasks

    return _run


def test_valid_params_run(run):
    result, _out, _tasks = run()
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize('inputs', [(), ('--cohorts', 'COH1', '--sgids', 'CPG1')])
def test_exactly_one_of_cohorts_or_sgids_is_needed(run, inputs):
    result, _out, tasks = run(inputs=inputs)
    assert result.returncode != 0
    assert 'Provide exactly one of --cohorts or --sgids' in result.stdout + result.stderr
    assert not tasks


@pytest.mark.parametrize(
    ('args', 'message'),
    [
        (('--sequencing_type', 'wgs'), "--sequencing_type must be genome or exome, got 'wgs'"),
        (('--accessLevel', 'prod'), "--accessLevel must be full, standard or test, got 'prod'"),
        (('--register', 'maybe'), 'Parameter `register` with type Boolean cannot be assigned to maybe'),
        (('--qc_genotyped', 'maybe'), 'Parameter `qc_genotyped` with type Boolean cannot be assigned to maybe'),
        (('--genotyping_intervals', '-1'), '--genotyping_intervals must be 0 (the default for the sequencing type) or'),
        (('--pairhmm_threads', '0'), '--pairhmm_threads must be at least 1, got 0'),
        (('--downsample_reads', '0'), '--downsample_reads must be at least 1, got 0'),
        (('--register_wait_minutes', '0'), '--register_wait_minutes must be at least 1, got 0'),
    ],
)
def test_invalid_params_fail_early(run, args, message):
    result, _out, tasks = run(*args)
    assert result.returncode != 0
    assert message in result.stdout + result.stderr
    assert not tasks


@pytest.mark.parametrize(
    ('bands', 'message'),
    [
        ('[13, 0, 30]', '--reblock_gq_bands must be a list of positive GQs, got [13, 0, 30]'),
        ('[]', '--reblock_gq_bands must be a list of positive GQs, got []'),
    ],
)
def test_reblock_gq_bands_must_be_positive(run, tmp_path, bands, message):
    # a list, so it's set from config rather than the command line
    config = tmp_path / 'bands.config'
    config.write_text(f'params.reblock_gq_bands = {bands}\n')
    result, _out, tasks = run('-c', str(config))
    assert result.returncode != 0
    assert message in result.stdout + result.stderr
    assert not tasks


@pytest.mark.parametrize('name', list(REFERENCES))
def test_missing_required_params_fail_early(run, tmp_path, name):
    # every reference has a gs:// default, so unset it from config and leave it off the command line
    unset = tmp_path / 'unset.config'
    unset.write_text(f'params.{name} = null\n')
    result, _out, tasks = run('-c', str(unset), omit=(name,))
    assert result.returncode != 0
    assert f'Parameter `{name}` is required' in result.stdout + result.stderr
    assert not tasks


@pytest.mark.parametrize('name', ['gatk_container', 'picard_container', 'cpg_container'])
def test_workflow_needs_every_container(run, tmp_path, name):
    # the test profile sets placeholder containers, and an empty value can't be passed on the command line
    unset = tmp_path / 'unset.config'
    unset.write_text(f'params.{name} = null\n')
    result, _out, tasks = run('-c', str(unset))
    assert result.returncode != 0
    assert f'Parameter `{name}` is required' in result.stdout + result.stderr
    assert not tasks


@pytest.mark.parametrize('inputs', [('--cohorts', 'COH1,COH2'), ('--sgids', 'CPG1,CPG2')])
def test_cohorts_or_sgids_are_built_into_the_published_sheet(run, tmp_path, inputs):
    result, out, tasks = run(inputs=inputs)
    assert result.returncode == 0, result.stdout + result.stderr
    assert {t['name'] for t in tasks} == {'BUILD_SAMPLESHEET (standard)', 'PREFLIGHT'}
    # the run records, just the sample sheet, go to the analysis dataset only
    assert {p.name for p in run_dir(out).iterdir()} == {'samplesheet.tsv'}
    assert not (out / OTHER / 'genotyping_samplesheets').exists()
    assert (run_dir(out) / 'samplesheet.tsv').read_text() == (tmp_path / 'sheet.tsv').read_text()


def test_reusing_an_ar_guid_fails_before_any_task(run):
    result, out, _tasks = run()
    assert result.returncode == 0, result.stdout + result.stderr

    result, _out, tasks = run()
    assert result.returncode != 0
    console = result.stdout + result.stderr
    assert 'Run records for arGuid local-test already exist' in console
    assert f'{run_dir(out)}/samplesheet.tsv' in console
    assert not tasks

    # a new arGuid gets its own run records
    result, out, _tasks = run(ar_guid='second')
    assert result.returncode == 0, result.stdout + result.stderr
    assert (run_dir(out, 'second') / 'samplesheet.tsv').exists()


def test_row_without_a_dataset_fails_early(run, tmp_path):
    edit_sheet(tmp_path, 'CPG4', dataset='')
    result, _out, tasks = run()
    assert result.returncode != 0
    assert 'CPG4: a row needs dataset, the dataset the SG belongs to' in result.stdout + result.stderr
    assert not completed(tasks)


@pytest.mark.parametrize('column', ['cram', 'cram_ref', 'selfsm'])
def test_row_to_genotype_without_its_inputs_fails_early(run, tmp_path, column):
    edit_sheet(tmp_path, 'CPG3', **{column: ''})
    result, _out, tasks = run()
    assert result.returncode != 0
    assert f'CPG3: a row to genotype needs {column}' in result.stdout + result.stderr
    assert not completed(tasks)


@pytest.mark.parametrize(
    ('sg_id', 'column', 'args'),
    [('CPG1', 'cram_ref', ()), (QC_ONLY, 'gvcf_ref', ('--qc_genotyped', 'true'))],
)
def test_reference_without_a_fasta_suffix_fails_early(run, tmp_path, sg_id, column, args):
    edit_sheet(tmp_path, sg_id, **{column: 'gs://bucket/ref.fasta.gz.bak'})
    result, _out, tasks = run(*args)
    assert result.returncode != 0
    console = result.stdout + result.stderr
    assert f"{sg_id}: reference 'gs://bucket/ref.fasta.gz.bak' must end in .fasta, .fa or .fna" in console
    assert not completed(tasks)


def test_reference_of_a_genotyped_sg_is_only_checked_if_it_is_qcd(run, tmp_path):
    edit_sheet(tmp_path, QC_ONLY, gvcf_ref='gs://bucket/ref.fasta.gz.bak')
    result, _out, _tasks = run()
    assert result.returncode == 0, result.stdout + result.stderr


def test_existing_outputs_fail_the_run_and_every_clash_is_listed(run, tmp_path):
    out = tmp_path / 'results'
    clashes = [
        # an unregistered gVCF at an SG's destination is a clash: the run only registers files it wrote
        gvcf_dir(out, PROJECT) / 'CPG1.g.vcf.gz',
        gvcf_dir(out, PROJECT) / 'CPG1.g.vcf.gz.tbi',
        gvcf_qc_dir(out, OTHER, 'CPG3') / 'CPG3.bcftools-stats',
        gvcf_qc_dir(out, PROJECT, QC_ONLY) / f'{QC_ONLY}.variant_calling_summary_metrics',
    ]
    ignored = [
        gvcf_qc_dir(out, OTHER, 'CPG9') / 'CPG9.bcftools-stats',  # not in the run
        gvcf_qc_dir(out, OTHER, 'CPG4') / 'CPG4.bcftools-stats',  # genotyped, with gVCF QC: nothing to publish
    ]
    for path in clashes + ignored:
        touch(path)
    result, _out, tasks = run('--qc_genotyped', 'true')
    assert result.returncode != 0
    # log.warn/log.error go to stderr, or to stdout when nextflow detects an AI agent (e.g. CLAUDECODE is set)
    console = result.stdout + result.stderr
    assert f'{len(clashes)} output(s) already exist' in console
    for path in clashes:
        assert str(path) in console
    assert 'CPG9' not in console
    assert 'CPG4' not in console
    assert not completed(tasks)


def test_exome_destinations_are_under_exome(run, tmp_path):
    out = tmp_path / 'results'
    touch(gvcf_dir(out, PROJECT) / 'CPG1.g.vcf.gz')  # a genome's destination, so no clash for an exome run
    clashes = [
        gvcf_dir(out, PROJECT, 'exome') / 'CPG1.g.vcf.gz',
        gvcf_qc_dir(out, OTHER, 'CPG3', 'exome') / 'CPG3.variant_calling_detail_metrics',
        gvcf_qc_dir(out, PROJECT, QC_ONLY, 'exome') / f'{QC_ONLY}.bcftools-stats',
    ]
    for path in clashes:
        touch(path)
    result, _out, _tasks = run('--qc_genotyped', 'true', sequencing_type='exome')
    assert result.returncode != 0
    console = result.stdout + result.stderr
    assert f'{len(clashes)} output(s) already exist' in console
    for path in clashes:
        assert str(path) in console


def test_cpg_flows_qc_files_are_not_clashes(run, tmp_path):
    # cpg-flow published gVCF QC to qc/, not gvcf_qc/<sg>/
    for suffix in GVCF_QC_SUFFIXES[:2]:
        touch(tmp_path / 'results' / PROJECT / 'qc' / f'{QC_ONLY}.{suffix}')
    result, _out, tasks = run('--qc_genotyped', 'true')
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'PREFLIGHT' in completed(tasks)


def test_genotyped_sgs_publish_nothing_without_qc_genotyped(run, tmp_path):
    # without --qc_genotyped, CPG5's gVCF QC destination is never written, so it can't clash
    touch(gvcf_qc_dir(tmp_path / 'results', PROJECT, QC_ONLY) / f'{QC_ONLY}.bcftools-stats')
    result, _out, _tasks = run()
    assert result.returncode == 0, result.stdout + result.stderr


def test_resume_reuses_the_sheet_and_skips_the_checks_it_passed(run):
    result, out, _tasks = run()
    assert result.returncode == 0, result.stdout + result.stderr

    # the run records exist, but -resume skips their check and takes the sheet from the cache. And as if the
    # first attempt had published CPG1's gVCF: the pre-flight check it passed is cached, so that's no clash
    touch(gvcf_dir(out, PROJECT) / 'CPG1.g.vcf.gz')
    result, _out, tasks = run('-resume')
    assert result.returncode == 0, result.stdout + result.stderr
    assert {(t['name'], t['status']) for t in tasks} == {
        ('BUILD_SAMPLESHEET (standard)', 'CACHED'),
        ('PREFLIGHT', 'CACHED'),
    }


def test_resume_checks_again_after_a_failed_pre_flight_check(run, tmp_path):
    # a failed check isn't cached, so a resumed run still stops on a file it didn't write
    clash = Path(touch(gvcf_dir(tmp_path / 'results', PROJECT) / 'CPG1.g.vcf.gz'))
    for args in [(), ('-resume',)]:
        result, _out, tasks = run(*args)
        assert result.returncode != 0
        assert '1 output(s) already exist' in result.stdout + result.stderr
        assert str(clash) in result.stdout + result.stderr
        assert not completed(tasks)
