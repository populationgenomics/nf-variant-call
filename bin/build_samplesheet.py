#!/usr/bin/env python3

"""
Build the genotyping sample sheet for one sequencing type, from a list of cohorts or a list of sequencing
group IDs. The SGs can come from any number of Metamist projects, and each row records its SG's dataset.

Run by the BUILD_SAMPLESHEET process at the start of main.nf (see docs/DESIGN.md §1). The input rules follow
nf-align's builder:
- cohorts are looked up across projects; missing (or not visible to this run), inactive and empty cohorts
  fail the run
- sequencing groups are de-duplicated across cohorts
- an SG's dataset is its sample's project, without any -test suffix. Its project must match the access
  level: a -test project at access level test, and any other project otherwise. Mismatches fail the run
- SGs from cohorts that aren't short-read, or aren't of the requested sequencing type, are skipped with a warning
- SG IDs are looked up across projects; missing SGs (or ones not visible to this run), and SGs that aren't
  short-read of the requested sequencing type, fail the run, as they were asked for by name
And these genotyping rules:
- an SG with a completed, active 'gvcf' analysis is genotyped (the latest one counts), and is never genotyped
  again. Its row carries the gVCF, its reference (the one the analysis records, or --genotyped-ref-fasta for a
  gVCF cpg-flow made), and the ID of the latest completed, active 'gvcf_qc' analysis whose primary file is that
  same gVCF, or nothing. It needs no CRAM. With --qc_genotyped, main.nf QCs a genotyped SG without gVCF QC
- any other SG is to genotype. It needs a completed, active 'cram' analysis, a 'cram_qc' analysis whose primary
  file is that same CRAM, and a 'verifybamid' secondary file on that cram_qc analysis, for the contamination
  estimate. A gap in any fails the run, in both modes: genotyping always follows alignment and CRAM QC. Its row
  carries the CRAM, its reference (the one the analysis records, or --aligned-ref-fasta for a CRAM made before
  nf-align) and the VerifyBamID .selfSM
- locations are used exactly as recorded in Metamist

    build_samplesheet.py --access-level standard --sequencing-type genome --cohorts COH0001 COH0002 \\
        --aligned-ref-fasta gs://.../masked.fasta --genotyped-ref-fasta gs://.../masked.fasta --output samplesheet.tsv
    build_samplesheet.py --access-level test --sequencing-type genome --sgids CPG0001 CPG0002 \\
        --aligned-ref-fasta gs://.../masked.fasta --genotyped-ref-fasta gs://.../masked.fasta --output samplesheet.tsv
"""

import csv
import logging
from argparse import ArgumentParser, Namespace
from pathlib import Path

from metamist.graphql import gql, query

logger = logging.getLogger(__name__)

COLUMNS = ['sg_id', 'dataset', 'cram', 'cram_ref', 'selfsm', 'gvcf', 'gvcf_ref', 'gvcf_qc']
SEQUENCING_TYPES = {'genome', 'exome'}
ACCESS_LEVELS = {'full', 'standard', 'test'}

# the SG fields both queries fetch: any existing CRAM, CRAM QC, gVCF and gVCF QC analyses, with their secondary files
SG_FIELDS = """
fragment SgFields on GraphQLSequencingGroup {
    id
    type
    technology
    sample {
        project {
            name
        }
    }
    analyses(type: {in_: ["cram", "cram_qc", "gvcf", "gvcf_qc"]}, status: {eq: COMPLETED}, active: {eq: true}) {
        id
        type
        outputs
        meta
    }
}
"""

# one round trip: cohort status, and each cohort's SGs. Top-level, so the cohorts can be in any project
COHORTS_QUERY = (
    """
query BuildSamplesheetFromCohorts($cohorts: [String!]!) {
    cohorts(id: {in_: $cohorts}) {
        id
        status
        sequencingGroups {
            ...SgFields
        }
    }
}
"""
    + SG_FIELDS
)

# one round trip: the requested SGs, queried directly, in any project
SGS_QUERY = (
    """
query BuildSamplesheetFromSgs($sgids: [String!]!) {
    sequencingGroups(id: {in_: $sgids}) {
        ...SgFields
    }
}
"""
    + SG_FIELDS
)


class SamplesheetError(Exception):
    """Invalid cohorts, sequencing group IDs, or sequencing group metadata."""


def fetch_cohort_sgs(cohort_ids: list[str]) -> list[dict]:
    """
    The SGs of these cohorts, in cohort order and possibly repeated. Fails on any cohort that is missing
    (or not visible to this run's credentials), inactive, or empty.
    """
    cohorts = query(gql(COHORTS_QUERY), variables={'cohorts': cohort_ids})['cohorts']

    missing = sorted(set(cohort_ids) - {cohort['id'] for cohort in cohorts})
    inactive = sorted(cohort['id'] for cohort in cohorts if cohort['status'] != 'active')
    empty = sorted(cohort['id'] for cohort in cohorts if not cohort['sequencingGroups'])
    if missing or inactive or empty:
        raise SamplesheetError(
            f'Unusable cohorts - missing or not visible to this run: {missing}, inactive: {inactive}, empty: {empty}',
        )
    return [sg for cohort in cohorts for sg in cohort['sequencingGroups']]


def fetch_sgs(sg_ids: list[str], sequencing_type: str) -> list[dict]:
    """
    The SGs with these IDs. Each was asked for by name, so one that is missing (or not visible to this run's
    credentials), or isn't short-read of this sequencing type, fails the run rather than being skipped.
    """
    sgs = query(gql(SGS_QUERY), variables={'sgids': sg_ids})['sequencingGroups']

    missing = sorted(set(sg_ids) - {sg['id'] for sg in sgs})
    wrong_type = sorted(
        f'{sg["id"]} ({sg["technology"]} {sg["type"]})'
        for sg in sgs
        if sg['technology'] != 'short-read' or sg['type'] != sequencing_type
    )
    if missing or wrong_type:
        raise SamplesheetError(
            f'Unusable SGs - missing or not visible to this run: {missing}, '
            f'not short-read {sequencing_type}: {wrong_type}',
        )
    return sgs


def get_primary_path(analysis: dict) -> str | None:
    """Get the primary path for an analysis - outputs.path."""
    # null for an analysis without a file
    return (analysis.get('outputs') or {}).get('path')


def latest_analysis(analyses: list[dict], analysis_type: str, primary: str | None = None) -> dict | None:
    """
    The latest analysis of this type that has an output, if any. With primary, only one whose primary file is
    that path counts, so QC of an earlier CRAM or gVCF doesn't.
    """
    found = sorted(
        (
            a
            for a in analyses
            if a['type'] == analysis_type and get_primary_path(a) and primary in {None, get_primary_path(a)}
        ),
        key=lambda a: a['id'],
    )
    return found[-1] if found else None


def secondary_file(analysis: dict, secondary_type: str) -> str:
    """
    The path of one of an analysis's secondary files, or ''. cpg_utils.metamist_registration's
    --secondary <type>=<path> puts it at outputs.secondary_files.<type>.path.
    """
    secondary = ((analysis.get('outputs') or {}).get('secondary_files') or {}).get(secondary_type) or {}
    return secondary.get('path') or ''


def recorded_ref(analysis: dict, default: str) -> str:
    """The reference an analysis records in its meta, or default for one registered before it was recorded."""
    return (analysis.get('meta') or {}).get('ref_fasta') or default


def row_for(sg: dict, aligned_ref_fasta: str, genotyped_ref_fasta: str) -> dict:
    """
    The sample sheet row for one SG: a genotyped SG's gVCF, or a to-genotype SG's CRAM and contamination estimate.
    aligned_ref_fasta and genotyped_ref_fasta are the references of a CRAM or gVCF whose analysis records none.
    """
    sg_id, analyses = sg['id'], sg['analyses']
    if gvcf := latest_analysis(analyses, 'gvcf'):
        path = get_primary_path(gvcf)
        gvcf_qc = latest_analysis(analyses, 'gvcf_qc', primary=path)
        return {
            'sg_id': sg_id,
            'gvcf': path,
            'gvcf_ref': recorded_ref(gvcf, genotyped_ref_fasta),
            # empty if this gVCF has no gVCF QC registered
            'gvcf_qc': str(gvcf_qc['id']) if gvcf_qc else '',
        }

    cram = latest_analysis(analyses, 'cram')
    if not cram:
        raise SamplesheetError(f'{sg_id}: to genotype, but has no CRAM. Align it with nf-align first')
    path = get_primary_path(cram)
    cram_qc = latest_analysis(analyses, 'cram_qc', primary=path)
    if not cram_qc:
        raise SamplesheetError(f'{sg_id}: to genotype, but its CRAM {path} has no CRAM QC. Run nf-align --qc_aligned')
    selfsm = secondary_file(cram_qc, 'verifybamid')
    if not selfsm:
        raise SamplesheetError(
            f'{sg_id}: to genotype, but the CRAM QC of {path} (analysis {cram_qc["id"]}) has no verifybamid file',
        )
    return {
        'sg_id': sg_id,
        'cram': path,
        'cram_ref': recorded_ref(cram, aligned_ref_fasta),
        'selfsm': selfsm,
    }


def build_rows(
    sgs: list[dict],
    access_level: str,
    sequencing_type: str,
    aligned_ref_fasta: str,
    genotyped_ref_fasta: str,
) -> list[dict]:
    """
    One row for every in-scope SG, each SG once, with the SG's dataset. Reports every bad SG at once, not one
    per run. aligned_ref_fasta and genotyped_ref_fasta are the references of a CRAM or gVCF whose analysis
    records none.
    """
    rows: list[dict] = []
    seen: set[str] = set()
    errors: list[str] = []
    for sg in sgs:
        if sg['id'] in seen:
            continue
        seen.add(sg['id'])
        # a test run must never publish to a -main bucket, nor a full run to a -test one
        project = sg['sample']['project']['name']
        if project.endswith('-test') != (access_level == 'test'):
            errors.append(f'{sg["id"]}: project {project} does not match access level {access_level}')
            continue
        if sg['technology'] != 'short-read' or sg['type'] != sequencing_type:
            logger.warning(
                f'{sg["id"]}: skipping {sg["technology"]} {sg["type"]}, this run is short-read {sequencing_type}'
            )
            continue
        try:
            row = row_for(sg, aligned_ref_fasta, genotyped_ref_fasta)
        except SamplesheetError as error:
            errors.append(str(error))
            continue
        rows.append({'dataset': project.removesuffix('-test')} | row)
    if errors:
        raise SamplesheetError('\n'.join(errors))
    if not rows:
        raise SamplesheetError(f'No short-read {sequencing_type} sequencing groups to genotype')
    return rows


def write_sheet(rows: list[dict], path: str) -> None:
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, COLUMNS, delimiter='\t', restval='')
        writer.writeheader()
        writer.writerows(rows)


def summary(rows: list[dict]) -> str:
    """How many SGs are to genotype, genotyped, and genotyped without gVCF QC."""
    genotyped = [row for row in rows if row.get('gvcf')]
    unqcd = sum(not row.get('gvcf_qc') for row in genotyped)
    datasets = len({row['dataset'] for row in rows})
    return (
        f'{len(rows)} SGs from {datasets} dataset(s): {len(rows) - len(genotyped)} to genotype, '
        f'{len(genotyped)} already genotyped ({unqcd} of them without gVCF QC)'
    )


def main(args: Namespace) -> None:
    if args.cohorts:
        sgs = fetch_cohort_sgs(sorted(set(args.cohorts)))
    else:
        sgs = fetch_sgs(sorted(set(args.sgids)), args.sequencing_type)
    rows = build_rows(sgs, args.access_level, args.sequencing_type, args.aligned_ref_fasta, args.genotyped_ref_fasta)
    write_sheet(rows, args.output)
    logger.info(f'{args.sequencing_type} at access level {args.access_level}: {summary(rows)} -> {args.output}')


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        '--access-level',
        required=True,
        choices=sorted(ACCESS_LEVELS),
        help='test: every SG must be in a <dataset>-test project; full or standard: none may be',
    )
    parser.add_argument('--sequencing-type', required=True, choices=sorted(SEQUENCING_TYPES))
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--cohorts', nargs='+', help='Metamist cohort IDs')
    inputs.add_argument('--sgids', nargs='+', help='Metamist sequencing group IDs, queried directly')
    parser.add_argument(
        '--aligned-ref-fasta',
        required=True,
        help="the reference of a CRAM whose cram analysis doesn't record one",
    )
    parser.add_argument(
        '--genotyped-ref-fasta',
        required=True,
        help="the reference of a gVCF whose gvcf analysis doesn't record one",
    )
    parser.add_argument('--output', required=True, help='Sample sheet TSV to write')
    main(parser.parse_args())
