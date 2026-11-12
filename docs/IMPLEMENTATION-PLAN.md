# nf-variant-call: implementation plan

This plan builds the pipeline in [DESIGN.md](DESIGN.md) step by step. Each step ends in a state that's committed,
lint-clean and tested. Most steps can be tested offline: pytest for the scripts, and `nextflow run -stub` for
the wiring, using nf-align's harness. Real compute only comes in phases 0 and 9.

Sources to copy from:
- **nf-align** (`../nf-align`): run shape, tests, sample sheet builder, pre-flight, registration, config conventions.
- **cpg-flow's untested port** (`../cpg-flow-align-genotype/nextflow/`): the genotyping modules.
- **cpg-flow** (`../cpg-flow-align-genotype/src/align_genotype/jobs/genotype.py`, `jobs/picard.py`): the reference
  for the exact commands.

Conventions in every step:
- Each new process has a `stub:` block from the start.
- Each step adds its own tests. Don't leave them all for the end.
- `nextflow lint main.nf`, `uv run ruff check .` and `uv run pytest -rs` pass before each commit.

## Phase 0: de-risking and prerequisites (in parallel with phases 1–3)

These don't block the offline build, but 0.1 can change the shape of `GENOTYPE_INTERVAL`, so start it first.

**0.1 Test reading CRAMs from GCS (RISKS #3).** On Google Batch, in the test access level, run one bare
HaplotypeCaller task by hand, with the image and command from DESIGN §4, `-I gs://…cram --read-index gs://…crai`
and a single interval list. Do it for two CRAMs: one from nf-align (unmasked reference) and one older CRAM
(DRAGEN-masked reference, interval list built against the unmasked `.dict`). For one of them, also run with the
CRAM copied to local disk.
- *Done when:* both runs read from GCS succeed, the local-disk run's gVCF matches the GCS run's (ignoring
  header lines), and the runtime and network costs are noted.
- *If it fails:* change DESIGN §4 to copy the CRAM to local disk and use fewer genotyping intervals, update
  RISKS #3, and write a short ADR.
- *Also check:* whether GATK accepts an interval list whose `.dict` md5s differ from the reference's. If it doesn't,
  `MAKE_INTERVALS` needs one set of lists per reference (see 4.2).

**0.2 External prerequisites.** None of these is code in this repo.
- Create the `gvcf_qc` analysis type in Metamist, as `cram_qc` was.
- Convert cpg-flow's `genome_evaluation_interval_lists` to BED and store it next to nf-align's
  `exome_evaluation_regions`.
- Check that `<name>.dict` and `.fai` exist beside the DRAGEN-masked FASTA, and create them if not.
- Confirm that the cpg-flow `gatk:4.6.2.0` image contains `bcftools`, `tabix` and `gatk MergeVcfs`, and pin its
  exact tag.

## Phase 1: scaffold

**1.1 Repo tooling.** Copy these from nf-align and rename them:
- `pyproject.toml`: name `nf-variant-call`, the same dependencies, ruff and pytest settings.
- `.github/workflows/{lint,nf_lint,test}.yaml`.
- `.gitignore` entries for `work/`, `results/` and `.nextflow*`.
- Run `uv lock`.

The test workflow also installs `bcftools`, which phase 5's filter test needs.

**1.2 Nextflow skeleton.**
- `nextflow.config`: the manifest, every param from DESIGN's two tables with its default, and the `standard`
  and `test` profiles (adapted from nf-align: Google Batch, spot, `maxSpotAttempts = 3`, and the `test` profile
  setting `output_root = 'results'`, `register = false` and stub containers).
- `conf/base.config`: the error strategy (137/140 retried twice, otherwise ignore; registration retried once),
  per-attempt `machineType` closures, and `dataset` resource labels. No processes yet.
- `main.nf`: the typed `params {}` block, `validateParams()`, and an empty workflow.
- `nextflow_schema.json`.

Copy these helpers from nf-align's `main.nf` unchanged: `typePrefix`, `metamistProject`, `datasetRoot`,
`runDir`, `uriOf`, `idList`, `fail`, `refDict`.
- *Tests:* port nf-align's `test/conftest.py` (the `run` fixture and the Metamist-shaped builders), and its
  param tests (`test_invalid_params_fail_early`, `test_missing_required_params_fail_early`,
  `test_exactly_one_of_cohorts_or_sgids_is_needed`, `test_workflow_needs_every_container`).
- *Done when:* `nextflow lint` is clean and the param tests pass.

## Phase 2: sample sheet

**2.1 `bin/build_samplesheet.py`.** Start from nf-align's builder: the GraphQL query, cohort and SG ID modes,
access-level checks, and errors collected into one report. Change the query to fetch each SG's `cram`,
`cram_qc`, `gvcf` and `gvcf_qc` analyses, with `outputs` and `meta`. Then add these rules from DESIGN §1:
- Each SG is either **genotyped** (it has an active `gvcf`; take the latest) or **to genotype**.
- **To genotype** requires a CRAM, a `cram_qc` whose primary file is that CRAM, and a `verifybamid` secondary
  file on it. Any gap fails the run, with `--cohorts` and with `--sgids`.
- **Genotyped** requires nothing else. `gvcf_qc` is the ID of an active `gvcf_qc` analysis whose primary file is
  this gVCF, or empty. Back-filled analyses without `bcftools_stats` count.
- `cram_ref` and `gvcf_ref` come from `ref_fasta` in each analysis's meta. Otherwise they're
  `--aligned-ref-fasta` and `--genotyped-ref-fasta`.
- Out-of-scope SGs are skipped with a warning (`--cohorts`) or fail the run (`--sgids`).
- The run fails if no SG is in scope.
- The counts are logged.
- The columns are `sg_id, dataset, cram, cram_ref, selfsm, gvcf, gvcf_ref, gvcf_qc`.

*Tests* (`test/test_build_samplesheet.py`, with fake GraphQL responses):
- one test per rule above
- an SG whose `cram_qc` points at a different CRAM
- an SG with two `gvcf_qc` analyses on different gVCFs
- SGs from two datasets
- the `-test` project mapping

**2.2 `BUILD_SAMPLESHEET` module and run-records check.** Copy nf-align's module, including its
`ext.stub_samplesheet` stub, and `checkNoRunRecords()`, pointed at `genotyping_samplesheets/<arGuid>/`. Publish
the sheet there.
- *Tests:* stub runs check the published sheet, and that reusing an `arGuid` fails before any task.

## Phase 3: row parsing, destinations, pre-flight

**3.1 `parseRow` and destinations in `main.nf`.**
- `parseRow` validates each row:
  - a `dataset` is present
  - a to-genotype row has `cram`, `cram_ref` and `selfsm`
  - a reference ends in `.fasta`, `.fa` or `.fna`
- It builds `meta` with these fields: `id`, `dataset`, `project`, `cram`, `ref`, `selfsm`, `gvcf`, `gvcfDir`,
  `qcDir`.
- Add `gvcfFiles(meta)`, `gvcfQcFiles(meta)` and `outputs(meta)`, keyed by registration key, as the single source
  of the layout in DESIGN §5.
- Branch the rows into three sets: genotyped, QC-only (`--qc_genotyped` and empty `gvcf_qc`), and to genotype.
- Add `publishDir` rules to `conf/base.config` (`mode: 'copy'`, `overwrite: false`, destinations from `meta`).

*Tests:*
- rows with no dataset, a to-genotype row missing `cram`/`selfsm`/`cram_ref`, and a bad reference suffix all fail
  early
- exome destinations go under `exome/`
- `accessLevel test` resolves to `-test` buckets

**3.2 `PREFLIGHT`.** Copy nf-align's native `exec:` module. Feed it every to-genotype and QC-only row, with
`outputs(meta)`. Genotyped SGs that aren't QC-only are left out.
- *Tests:* ported from nf-align:
  - existing outputs fail the run and every clash is listed
  - a resume skips a pre-flight check that passed
  - a resume checks again after a failed one
  - an unregistered gVCF is a clash
  - cpg-flow's `qc/<sg>.variant_calling_*` files are not

## Phase 4: GENOTYPE subworkflow

**4.1 `bin/freemix.sh`.** Read the `FREEMIX` column, by header name not position, from a VerifyBamID2 `.selfSM`,
and print it. Exit non-zero if the file is missing or empty, if the column is absent, or if the value isn't a
number in [0, 1].
- *Tests:* pytest over fixture `.selfSM` files, including one with columns reordered, and an uncapped high value
  such as 0.35, which is passed through.

**4.2 `MAKE_INTERVALS`.** Port `nextflow/modules/local/intervals` from cpg-flow:
- picard image, the command in DESIGN §3
- N is resolved in `main.nf` (0 means 30 for genomes, 1 for exomes)
- exclusions are always subtracted
- outputs are named `0001.interval_list`, `0002.interval_list` and so on, so they sort into genomic order

If 0.1 found that the `.dict` md5s matter, this runs once per distinct reference, keyed by reference.
- *Tests:* a stub run makes exactly N lists, for N = 1, 3 and the defaults.

**4.3 `GENOTYPE_INTERVAL`.** Port `nextflow/modules/local/haplotypecaller` from cpg-flow:
- Input is `[meta, cram_uri (val), crai_uri (val), interval]`. The reference, the `.selfSM` and
  `bin/freemix.sh` are staged.
- Drop `--dragen-mode`.
- Add `--contamination-fraction-to-filter $(freemix.sh <selfsm>)`.
- Downsampling applies to genomes only.
- The Java heap is the task's memory minus headroom.
- The output is `<sg>.<i>.g.vcf.gz` and `.tbi`.

*Tests:*
- the stub's `.command.sh` has no `--dragen-mode`, and has the contamination flag and `-L`
- the CRAM is not staged into the task directory

**4.4 `bin/postproc_gvcf.sh`.**
- If more than one interval gVCF is given, run `gatk MergeVcfs` over them in the order given.
- Then run cpg-flow's postprocessing steps verbatim (from `genotype.py:303-353`), with the GQ bands as arguments.
- Print the GATK version to `gatk_version.txt`, and FREEMIX (via `freemix.sh`) to `contamination.txt`.

*Tests:* where `bcftools` is available, a small fixture gVCF:
- missing chr21 fails the contig check
- `INFO/DP` and `INFO/DS` are removed
- the sample is renamed to the SG ID

GATK steps aren't run in tests.

**4.5 `POSTPROC_GVCF` and the subworkflow.**
- The module runs `postproc_gvcf.sh` in the gatk image and emits `[meta, gvcf, tbi]`. The version and
  contamination are read into `meta` (`meta + [gatk_version: …, contamination: …]`).
- `subworkflows/local/genotype/main.nf` combines the CRAMs with the intervals, groups them with
  `groupKey(meta, N)` and `groupTuple(sort: …)` by interval index, and runs `POSTPROC_GVCF`.
- It takes `[meta, cram, crai]` plus the intervals, as in DESIGN §4.

*Tests:*
- with N = 3, each SG gives one published gVCF and `.tbi` under `[exome/]gvcf/`
- an injected failure in one SG's interval drops only that SG
- `-resume` then finishes it

## Phase 5: gVCF QC

**5.1 `BCFTOOLS_STATS`.** The variant-sites filter drops a record only when `<NON_REF>` is its only ALT allele,
for example `bcftools view -e 'N_ALT=1 && ALT="<NON_REF>"'`. Confirm the exact expression against the fixture
before trusting it. Then run `bcftools stats -F <ref> -` into `<sg>.bcftools-stats`. Use the gatk image.

*Tests:* a fixture gVCF with:
- a reference block (`<NON_REF>`)
- a SNV (`C,<NON_REF>`)
- a multi-allelic site (`C,T,<NON_REF>`)
- an indel

Run the real filter command from the module script on it (pytest, with bcftools needed). Exactly the three
variant records survive.

**5.2 `VARIANT_CALLING_METRICS`.** Model it on nf-align's `picard_metrics`: `BedToIntervalList` with the gVCF's
own reference `.dict`, then `CollectVariantCallingMetrics` with dbSNP, `-SD` and `--GVCF_INPUT true`. The
outputs are `<sg>.variant_calling_{summary,detail}_metrics`.
- *Tests:* the stub's `.command.sh` uses the sequencing type's evaluation BED, and the `.dict` of the gVCF's
  reference, not the alignment reference.

**5.3 `GVCF_QC` subworkflow.** Join the two outputs by SG, so an SG comes out only if both succeeded, and publish
to `[exome/]gvcf_qc/<sg>/`. Wire it to new gVCFs (with `meta.ref`) and to QC-only SGs (with the row's
`gvcf_ref`).

*Tests:*
- a QC-only SG is QC'd without being genotyped
- a run of only genotyped SGs with `--qc_genotyped` does QC only
- an injected QC failure publishes and registers nothing for that SG's QC

## Phase 6: registration

**6.1 `REGISTER_ANALYSIS`.** Copy nf-align's module unchanged: it polls for files, then calls
`cpg_utils.metamist_registration`. Include it twice, as `REGISTER_GVCF` and `REGISTER_GVCF_QC`.

**6.2 Wiring.**
- **`REGISTER_GVCF`** takes the gVCF, with `index` as a secondary file. Its meta is `sequencing_type`,
  `ar-guid`, `ref_fasta`, `cram`, `caller=gatk-<version>`, `dragen_mode=false` and `contamination`.
- **`REGISTER_GVCF_QC`** is joined on both the QC output and `REGISTER_GVCF`'s completion. For a QC-only SG it
  needs only the QC output. Its secondary keys are `variant_calling_summary_metrics`,
  `variant_calling_detail_metrics` and `bcftools_stats`.
- Neither runs with `--register false`.

*Tests:* ported from nf-align's registration tests:
- the published path, and the SG's own project (including `-test`)
- one `gvcf_qc` per gVCF
- `gvcf_qc` is never registered when `gvcf` registration failed
- nothing is registered with `--register false`
- per-SG tasks carry the SG's `dataset` label

## Phase 7: resources and docs

**7.1 Resources.** Fill in `conf/base.config` from DESIGN's resources table: exact `machineType` per attempt,
memory scaled by attempt, and disk by sequencing type. The table is missing a row for
`VARIANT_CALLING_METRICS`, so add one: 2 CPUs, 8 GB, 50/20 GB disk and `e2-highmem-2/4/8`, to be confirmed.
- *Test:* port nf-align's label test, and add a check that every process has a `machineType`.

**7.2 Docs.**
- Update DESIGN's Status section to "done, stub-tested", and list the judgement calls made while building.
- Add `docs/AMBIGUOUS-CHOICES.md` if any came up.
- Expand the README's Running section.

## Phase 8: combining hook (optional now)

Check that `GENOTYPE` and `GVCF_QC` can be included from another pipeline without depending on `main.nf`'s
helper functions: everything they need comes in through `meta` and their inputs. Fix any leaks. Don't start
merging into nf-align yet.

## Phase 9: real runs and validation

**9.1 Test-level run.** In the test access level, run a small cohort with SGs from two datasets: genome and
exome, one older CRAM, and one already-genotyped SG with `--qc_genotyped`. Check:
- the output paths and the registrations (meta, secondary files, projects)
- that `-resume` after a deliberately failed task works
- that `nf-combine`'s sample sheet builder accepts the new gVCFs
- the cost per SG, from Seqera

**9.2 GIAB benchmark (ADR 0001).** Align a GIAB sample with nf-align and genotype it here. Run hap.py against
truth and against its DRAGEN-mode cpg-flow gVCF. Record the results in DESIGN and update RISKS #2.

**9.3 Tuning.** Adjust the default number of genotyping intervals and the resources from 9.1's numbers.

**9.4 Production.** Run the first production cohort.

## Follow-ups (other repos and tasks)

- Back-fill cpg-flow's gVCF QC results as `gvcf_qc` analyses, using the same secondary keys and no
  `bcftools_stats`.
- nf-combine decides whether a VDS may mix genotyping modes (RISKS #1).
- Commit the shared glossary changes in nf-align, and consider retitling it now that three repos use it.
