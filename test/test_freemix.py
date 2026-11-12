"""bin/freemix.sh over fixture .selfSM files: FREEMIX is found by its header name and passed through as is."""

import subprocess
from pathlib import Path

import pytest

FREEMIX = Path(__file__).parents[1] / 'bin' / 'freemix.sh'

# VerifyBamID2's header, as nf-align's VERIFYBAMID writes it
HEADER = [
    '#SEQ_ID', 'RG', 'CHIP_ID', '#SNPS', '#READS', 'AVG_DP', 'FREEMIX', 'FREELK1', 'FREELK0', 'FREE_RH', 'FREE_RA',
    'CHIPMIX', 'CHIPLK1', 'CHIPLK0', 'CHIP_RH', 'CHIP_RA', 'DPREF', 'RDPHET', 'RDPALT',
]  # fmt: skip


def selfsm(tmp_path: Path, freemix: str = '0.00123', header: list[str] = HEADER, rows: int = 1) -> Path:
    values = dict.fromkeys(header, 'NA') | {'#SEQ_ID': 'CPG1', '#SNPS': '9000', 'FREEMIX': freemix}
    path = tmp_path / 'CPG1.verify-bamid.selfSM'
    lines = ['\t'.join(header)] + ['\t'.join(values[name] for name in header)] * rows
    path.write_text('\n'.join(lines) + '\n')
    return path


def freemix(path: Path) -> subprocess.CompletedProcess:
    return subprocess.run([str(FREEMIX), str(path)], capture_output=True, text=True, check=False)


@pytest.mark.parametrize('value', ['0.00123', '0', '1', '0.35', '1e-05', '0.0'])
def test_freemix_is_printed_as_written(tmp_path, value):
    # no cap: a high value such as 0.35 is passed through for HaplotypeCaller to use
    result = freemix(selfsm(tmp_path, value))
    assert result.returncode == 0, result.stderr
    assert result.stdout == f'{value}\n'


def test_freemix_is_found_by_its_header_name(tmp_path):
    reordered = ['FREEMIX', *(name for name in HEADER if name != 'FREEMIX')][::-1]
    result = freemix(selfsm(tmp_path, '0.042', header=reordered))
    assert result.returncode == 0, result.stderr
    assert result.stdout == '0.042\n'


def test_missing_file_fails(tmp_path):
    result = freemix(tmp_path / 'absent.selfSM')
    assert result.returncode != 0
    assert 'is missing or empty' in result.stderr


def test_empty_file_fails(tmp_path):
    path = tmp_path / 'empty.selfSM'
    path.touch()
    result = freemix(path)
    assert result.returncode != 0
    assert 'is missing or empty' in result.stderr


def test_missing_column_fails(tmp_path):
    result = freemix(selfsm(tmp_path, header=[name for name in HEADER if name != 'FREEMIX']))
    assert result.returncode != 0
    assert 'no FREEMIX column' in result.stderr
    assert not result.stdout


@pytest.mark.parametrize('rows', [0, 2])
def test_other_than_one_data_row_fails(tmp_path, rows):
    result = freemix(selfsm(tmp_path, rows=rows))
    assert result.returncode != 0
    assert f'expected one data row, found {rows}' in result.stderr
    assert not result.stdout


@pytest.mark.parametrize('value', ['NA', '', '-0.1', '1.5', 'nan', '0.1.2'])
def test_value_outside_0_to_1_fails(tmp_path, value):
    result = freemix(selfsm(tmp_path, value))
    assert result.returncode != 0
    assert 'FREEMIX is not a number in [0, 1]' in result.stderr
    assert not result.stdout
