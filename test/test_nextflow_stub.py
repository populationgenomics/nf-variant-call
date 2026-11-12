"""Offline `nextflow run -stub` checks of main.nf's wiring: params, the run-records check and the published sample
sheet. The SGs span two datasets (conftest.cohort_sgs). Placeholder files stand in for data, and the tests are
skipped when nextflow isn't on PATH."""

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
    assert [t['name'] for t in tasks] == ['BUILD_SAMPLESHEET (standard)']
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


def test_resume_skips_the_run_records_check_and_reuses_the_sheet(run):
    result, _out, _tasks = run()
    assert result.returncode == 0, result.stdout + result.stderr

    result, _out, tasks = run('-resume')
    assert result.returncode == 0, result.stdout + result.stderr
    assert [(t['name'], t['status']) for t in tasks] == [('BUILD_SAMPLESHEET (standard)', 'CACHED')]
