# nf-variant-call: design brief

This pipeline genotypes the short-read sequencing groups of a set of Metamist cohorts, or of a list of
sequencing groups, which can come from several datasets. For each SG it produces a single-sample, reblocked
gVCF and a set of gVCF QC metrics. Both are published to the SG's own dataset and registered in it. On request,
SGs that are already genotyped get gVCF QC alone. It's a Nextflow (≥ 26.04, strict syntax) pipeline that
runs on Google Batch against data in GCS, and it can be tested offline with stub runs.

It follows nf-align's design wherever the two overlap. Its commands come from the genotyping part of
`cpg-flow-align-genotype`, and from the untested Nextflow port of it in that repo's `nextflow/`. It changes them
in a few places: no DRAGEN mode, contamination taken from CRAM QC, CRAMs read straight from GCS, and new gVCF QC.
For now the pipeline is standalone, but its `GENOTYPE` subworkflow is shaped so nf-align can include it and run
alignment and genotyping as one workflow (see [Combining with nf-align](#combining-with-nf-align)).

Terms are defined in [GLOSSARY.md](../GLOSSARY.md), which this pipeline shares with nf-align and nf-combine.
[RISKS.md](RISKS.md) lists what may bite later, and the decisions with real trade-offs are recorded in
[adr/](adr/).

**Status: design only.** Nothing is implemented yet. See [Status and next steps](#status-and-next-steps).

## Scope of a run

The scope is the same as nf-align's. A run is launched under one **analysis dataset** and covers **one
sequencing type**. Its SGs can come from any number of datasets. Each SG is published to its own dataset's
bucket and registered in that dataset's Metamist project.

| param                                                 | meaning                                                                                                                                                                           |
|-------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `dataset`                                             | the analysis dataset (no `-test` suffix). The run is launched under it, it pays for the run, and it holds the run records                                                         |
| `sequencing_type`                                     | `genome` or `exome`                                                                                                                                                               |
| `accessLevel`                                         | `full` / `standard` (equivalent), or `test`. Applies to the whole run                                                                                                             |
| `arGuid`                                              | the analysis-runner GUID. Mandatory                                                                                                                                               |
| `cohorts` **or** `sgids`                              | exactly one: Metamist cohort IDs, or Metamist SG IDs                                                                                                                              |
| `output_root`                                         | null by default. If set (the `test` profile sets `results`), outputs go to `<output_root>/<ds>/` instead of the dataset bucket                                                    |
| `aligned_ref_fasta`                                   | the reference of a CRAM whose `cram` analysis doesn't record one: the DRAGEN-masked FASTA, as in nf-align                                                                         |
| `genotyped_ref_fasta`                                 | the reference of a gVCF whose `gvcf` analysis doesn't record one (every gVCF cpg-flow made): the same DRAGEN-masked FASTA. Only used for QC-only SGs                              |
| `calling_regions_{genome,exome}`                      | the regions genotyped, as interval lists: cpg-flow's `{genome,exome}_calling_interval_lists`                                                                                      |
| `calling_exclusions`                                  | regions removed from the calling regions: cpg-flow's `hg38_telomeres_and_centromeres_intervals_list`                                                                              |
| `noalt_bed`                                           | the primary contigs plus chrM. Records outside it are dropped from the final gVCF                                                                                                 |
| `dbsnp_vcf`                                           | dbSNP 138, with its `.idx` alongside, for `CollectVariantCallingMetrics`                                                                                                          |
| `{genome,exome}_evaluation_regions`                   | the gVCF QC evaluation regions, as BED. The exome one is nf-align's `exome_evaluation_regions`. The genome one is cpg-flow's `genome_evaluation_interval_lists`, converted to BED |
| `gatk_container`, `picard_container`, `cpg_container` | the images. `gatk` is cpg-flow's `gatk:4.6.2.0` image, which also holds bcftools and tabix. `picard` is the image nf-align uses. `cpg` holds `cpg_utils` and `metamist`           |

A reference's `.fai` and `<name>.dict` must sit beside it. That applies to `aligned_ref_fasta` too: an SG whose
CRAM predates nf-align is genotyped against it.

Switches and tuning, all with defaults in `nextflow.config`:

| param                   | default            | meaning                                                                                      |
|-------------------------|--------------------|----------------------------------------------------------------------------------------------|
| `register`              | `true`             | `false` skips all Metamist registration (the `test` profile sets it)                         |
| `qc_genotyped`          | `false`            | also run gVCF QC on genotyped SGs whose gVCF has no `gvcf_qc` analysis                       |
| `register_wait_minutes` | 30                 | how long registration polls for the published files                                          |
| `genotyping_intervals`  | 0                  | how many genotyping intervals each SG is split into. 0 means 30 for genomes and 1 for exomes |
| `downsample_reads`      | 20                 | HaplotypeCaller's `--max-reads-per-alignment-start`, for genomes only                        |
| `pairhmm_threads`       | 4                  | HaplotypeCaller's `--native-pair-hmm-threads`                                                |
| `reblock_gq_bands`      | `[13, 20, 30, 40]` | `ReblockGVCF`'s GQ bands                                                                     |

As in nf-align, `main.nf`'s `params {}` block declares the types, `validateParams()` checks what the types can't,
there is no `outputDir`, and the launcher supplies `-work-dir` and `google.project`.

**Two references, two genotyping modes.** Each CRAM is genotyped against the reference it was written with.
That's unmasked GATK hg38 for nf-align CRAMs, and the DRAGEN-masked FASTA for older ones. The reference is
recorded in the `gvcf` analysis. Every gVCF made before this pipeline was genotyped in DRAGEN mode, and this
pipeline's aren't (ADR 0001). Both facts are in the analysis meta, so a VDS that mixes them can be split.
Whether a VDS may mix them is nf-combine's decision.

## Flow

```
run-records check ─ --cohorts | --sgids ─ BUILD_SAMPLESHEET ─┐   (cached on -resume)
                                                             ▼
sample sheet ─ parse rows, per-SG destinations into meta ─ PREFLIGHT (every row) ─┬─ genotyped: nothing more
                                                                                  ├─ QC-only (--qc_genotyped)
MAKE_INTERVALS (once per run) ─────────────────────────────────┐                  └─ to genotype
                                                               ▼
   to genotype: GENOTYPE: GENOTYPE_INTERVAL × N ─ groupKey ─ POSTPROC_GVCF ─ publishDir: gVCF, .tbi
                  ├─ REGISTER_GVCF (one gvcf analysis: gVCF + .tbi) ──────────────────┐
                  └─ GVCF_QC: VARIANT_CALLING_METRICS, BCFTOOLS_STATS ─ publishDir ───┴─ REGISTER_GVCF_QC

   QC-only:  GVCF_QC (the registered gVCF, its own reference) ─ publishDir ─ REGISTER_GVCF_QC
```

### 1. Sample sheet (`BUILD_SAMPLESHEET`, `bin/build_samplesheet.py`)

This works like nf-align's builder: one GraphQL query through the top-level `cohorts(id:)` or
`sequencingGroups(id:)` lookup, taking `--access-level` rather than a project. It fails the run on missing,
inactive or empty cohorts, on access-level mismatches, and when no SG is in scope. Every offender is reported
in one error. For each SG it fetches the project, the latest completed, active `cram`, `cram_qc` and `gvcf`
analyses, and the `gvcf_qc` analyses.

- **Out of scope:** with `--cohorts`, an SG that isn't short-read or isn't `sequencing_type` is skipped with a
  warning. With `--sgids`, it fails the run.
- **To genotype:** an SG with no active `gvcf` analysis. It needs:
  - an active `cram` analysis
  - a `cram_qc` analysis whose primary file is that CRAM
  - a `verifybamid` secondary file on that `cram_qc` analysis

  Missing any of these fails the run, in both modes. Genotyping always follows alignment and CRAM QC. The fix
  is to run nf-align, or nf-align with `--qc_aligned`.
- **Genotyped:** an SG with an active `gvcf` analysis. It's never genotyped again, whichever CRAM its gVCF
  came from (stale gVCFs are an accepted gap). Its `gvcf_qc` column holds the ID of a completed, active
  `gvcf_qc` analysis whose primary file is that gVCF, or is empty. Back-filled analyses of cpg-flow's QC files
  count. With `--qc_genotyped`, an empty `gvcf_qc` makes the SG QC-only. A genotyped SG needs no CRAM or
  `cram_qc`.
- Each row records the SG's `dataset`: its project, with any `-test` removed.
- It logs how many SGs are to genotype, how many are genotyped, and how many genotyped SGs have no gVCF QC.
- The sheet is published to `genotyping_samplesheets/<arGuid>/samplesheet.tsv` in the analysis dataset bucket.
  As in nf-align, every new run builds its own sheet, and a resumed run takes it from the cache.

The columns are `sg_id, dataset, cram, cram_ref, selfsm, gvcf, gvcf_ref, gvcf_qc`, one row per SG. A row with
`gvcf` set is a genotyped SG. Otherwise the SG is to genotype, and `cram`, `cram_ref` and `selfsm` are required.
- `cram_ref` is the `ref_fasta` in the CRAM analysis's meta, or else `aligned_ref_fasta`.
- `gvcf_ref` is the `ref_fasta` in the gVCF analysis's meta, or else `genotyped_ref_fasta`.
- The CRAM index is `<cram>.crai`, as nf-align publishes it.

### 2. Pre-flight

This is the same as nf-align. The pipeline never overwrites a file: `publishDir` uses `overwrite: false`, so
clashes are caught up front.

- **Run records:** the run fails if `genotyping_samplesheets/<arGuid>/samplesheet.tsv` already exists. This
  check is skipped on `-resume`.
- **Outputs (`PREFLIGHT`, a cached native `exec:` process):** before anything is genotyped, every URI the run
  would publish is checked: each to-genotype SG's gVCF, `.tbi` and gVCF QC files, and each QC-only SG's gVCF QC
  files. If any exist, the run fails and lists them all. An unregistered gVCF at an SG's destination is a clash:
  the run only registers files it wrote. cpg-flow's QC files live at `[exome/]qc/<sg>.variant_calling_*`, so
  they never clash with `[exome/]gvcf_qc/<sg>/`.

### 3. Genotyping intervals (`MAKE_INTERVALS`)

This runs once per run, in the picard image. It splits the sequencing type's calling regions, minus
`calling_exclusions`, into N lists. N is `genotyping_intervals`, or 30 for genomes and 1 for exomes:

```
picard IntervalListTools -SCATTER_COUNT <N> -SUBDIVISION_MODE INTERVAL_SUBDIVISION -UNIQUE true -SORT true \
    -BREAK_BANDS_AT_MULTIPLES_OF <100000 genome | 0 exome> -I <calling regions> -ACTION SUBTRACT -SI <exclusions> \
    -OUTPUT out
```

Exclusions are subtracted even when N is 1, unlike cpg-flow, so the regions genotyped never depend on N. Every SG
in the run gets the same lists. The task is cached, so `-resume` keeps them.

### 4. GENOTYPE subworkflow

It **takes** `[meta, cram, crai]` for each SG to genotype, plus the interval lists. Its meta carries `ref`,
the CRAM's reference, and `selfsm`, the VerifyBamID output. It **emits** `gvcf`: `[meta, gvcf, tbi]`.

**`GENOTYPE_INTERVAL`** (gatk image) runs once per SG per genotyping interval. It's cpg-flow's command without
`--dragen-mode`, with contamination added (ADR 0001):

```
gatk --java-options "<heap> -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" HaplotypeCaller \
  -R <ref> -I gs://.../<sg>.cram --read-index gs://.../<sg>.cram.crai -L <interval list> \
  --disable-spanning-event-genotyping --contamination-fraction-to-filter <FREEMIX> \
  --native-pair-hmm-threads <pairhmm_threads> \
  -O <sg>.<i>.g.vcf.gz -G AS_StandardAnnotation \
  -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
  -ERC GVCF --create-output-variant-index [--max-reads-per-alignment-start <downsample_reads>]   # genomes only
```

- **The CRAM is not staged.** Its `gs://` URI goes in as a `val`, and GATK reads it from GCS through its index,
  so each task fetches only its own region. The reference (FASTA, `.fai`, `.dict`), the interval list and the
  `.selfSM` are staged. Reading from GCS is untested on Batch (RISKS #3). If it fails, the fallback is to stage
  the CRAM and use fewer genotyping intervals.
- **FREEMIX** is read from the staged `.selfSM` by `bin/freemix.sh`, the only parser of that file. It's passed
  as is, with no cap. A `.selfSM` that's missing or unreadable fails the SG.
- `-L` is always passed, even when N is 1.

**`POSTPROC_GVCF`** (gatk image) runs once per SG, after its N interval gVCFs are grouped with
`groupKey(meta, N)`. If any interval failed, the group is incomplete and is dropped, so the SG fails. The steps
are:

1. If N > 1, `gatk MergeVcfs` the interval gVCFs, in interval order.
2. Run cpg-flow's postprocessing script unchanged:
   - `bcftools index --tbi`
   - fail unless chr1–22 and chrX are present (exomes too)
   - drop `INFO/DP`
   - `ReblockGVCF --reference <ref> -do-qual-approx --floor-blocks -GQB <each reblock_gq_bands>`
   - keep only `noalt_bed`
   - drop `INFO/DS`
   - rename the sample to the SG ID
   - `tabix`

The output is `<sg>.g.vcf.gz` and `.tbi`. The task also records the GATK version (`gatk --version`) for the
analysis meta.

### 5. Outputs

Everything is published with `publishDir`, `mode: 'copy'` and `overwrite: false`. `main.nf` works out each SG's
destinations into `meta` and lists every published file in `outputs(meta)`, which the pre-flight check and both
registrations read, as in nf-align.

| target  | path under the SG's output root (`exome/` prefix for exomes)                                                                  |
|---------|-------------------------------------------------------------------------------------------------------------------------------|
| gVCF    | `[exome/]gvcf/<sg>.g.vcf.gz`, plus `.tbi`: cpg-flow's path                                                                    |
| gVCF QC | `[exome/]gvcf_qc/<sg>/`: `<sg>.variant_calling_summary_metrics`, `<sg>.variant_calling_detail_metrics`, `<sg>.bcftools-stats` |

Only SGs genotyped by this run publish a gVCF. QC-only SGs publish only their gVCF QC files. Per-SG tasks are
labelled with `meta.dataset` for billing, and `BUILD_SAMPLESHEET`, `PREFLIGHT` and `MAKE_INTERVALS` with the
analysis dataset.

### 6. gVCF QC (`GVCF_QC`)

Every new gVCF, and every QC-only SG's registered gVCF, goes through two tasks against the gVCF's own
reference. That's `meta.ref` for a new gVCF, and the row's `gvcf_ref` for a QC-only SG.

- **`VARIANT_CALLING_METRICS`** (picard image): `BedToIntervalList` turns the sequencing type's evaluation BED
  into an interval list headed with the reference's own `.dict`, the same way nf-align's Picard QC avoids
  md5 mismatches. Then:
  ```
  picard CollectVariantCallingMetrics -I <gvcf> -O <sg> --DBSNP <dbsnp_vcf> -SD <ref .dict> \
      -TI <evaluation interval list> --GVCF_INPUT true
  ```
- **`BCFTOOLS_STATS`** (gatk image): runs on variant sites only. A record is dropped only when `<NON_REF>` is its
  only ALT allele (a reference block). A real variant's ALT is like `C,<NON_REF>`, so matching any `<NON_REF>`
  would drop every variant too. The filter is tested against a fixture that includes such a record. Output:
  `bcftools stats -F <ref> - > <sg>.bcftools-stats`.

An SG comes out of `GVCF_QC` only if both tasks succeeded, so gVCF QC is all or nothing. A QC failure never
affects the gVCF's registration: that SG can get QC later through `--qc_genotyped`.

### 7. Registration (`REGISTER_GVCF`, `REGISTER_GVCF_QC`)

Both use the `REGISTER_ANALYSIS` process taken from nf-align, included twice. Neither runs with
`--register false`. Both poll for up to `register_wait_minutes` for their published files.

`REGISTER_GVCF` runs once per newly genotyped SG, in the SG's own project:

```
python3 -m cpg_utils.metamist_registration --project <SG's dataset>[-test] --output <gvcf> --type gvcf \
    --secondary index=<gvcf>.tbi --sgs <sg> \
    --meta sequencing_type=<type> ar-guid=<arGuid> ref_fasta=<ref> cram=<source CRAM> \
           caller=gatk-<version> dragen_mode=false contamination=<FREEMIX>
```

`REGISTER_GVCF_QC` runs once the SG's gVCF QC is complete *and* its `gvcf` analysis exists (registered by
`REGISTER_GVCF`, or already, for a QC-only SG). So a `gvcf_qc` analysis never points at an unregistered gVCF:

```
python3 -m cpg_utils.metamist_registration --project <SG's dataset>[-test] --output <gvcf> --type gvcf_qc \
    --secondary variant_calling_summary_metrics=<dir>/<sg>.variant_calling_summary_metrics \
        variant_calling_detail_metrics=<dir>/<sg>.variant_calling_detail_metrics \
        bcftools_stats=<dir>/<sg>.bcftools-stats \
    --sgs <sg> --meta sequencing_type=<type> ar-guid=<arGuid>
```

`gvcf_qc` is a new Metamist analysis type, and it has to exist before a registering run. The back-fill of
cpg-flow's QC results is a separate job. It should use the same two `variant_calling_*` keys, and has no
`bcftools_stats`, so consumers must treat `bcftools_stats` as optional.

### 8. Failures

These are the same as nf-align:
- **Retries:** per-SG tasks are retried twice on exit codes 137 or 140, with memory and machine type scaled by
  attempt. Any other failure drops the SG and the rest carry on. Registration is retried once and then ignored.
- **Spot preemption** is retried by `maxSpotAttempts = 3`.
- **Run-level failures:** `BUILD_SAMPLESHEET`, `PREFLIGHT`, `MAKE_INTERVALS` and the run-records check fail the
  whole run.
- **No failure report:** the run exits 0, and Seqera shows the failed tasks.
- **Recovery is `-resume`.** A new run is different: an SG whose genotyping finished but whose registration
  failed has a published, unregistered gVCF, and a new run's pre-flight check stops on it until the gVCF is
  registered or deleted by hand.

## Combining with nf-align

`GENOTYPE` takes `[meta, cram, crai]` with `meta.ref` and `meta.selfsm`. In a combined workflow:
- **Newly aligned SGs:** these come from `ALIGN.out.cram`, joined with the `.selfSM` from `CRAM_QC`'s
  `VERIFYBAMID`. A CRAM whose QC failed never reaches genotyping, so it's just a failed SG.
- **Already-aligned SGs:** these come from the sample sheet, as here.

A freshly aligned CRAM lives in the `gs://` work directory, so reading it from GCS still works. The two sample
sheet builders and pre-flight checks would need merging. That's left until the workflows are combined.

## Accepted gaps

- **Stale gVCFs:** a genotyped SG is never genotyped again, even if its CRAM has since been replaced. Each
  new gVCF records its source CRAM, so stale ones can be found later.
- **Mixed genotyping modes:** VDSs may hold both DRAGEN-mode and non-DRAGEN gVCFs (RISKS #1).
- **Unvalidated configuration:** until the hap.py benchmark, this genotyping setup is unproven (ADR 0001,
  RISKS #2).
- **Back-filled gVCF QC lacks `bcftools stats`:** `bcftools_stats` is an optional key.
- **Inherited from nf-align:**
  - concurrent runs can both pass the pre-flight check
  - one stray file holds up a new run
  - no failure report
- **QC-only SGs** need a `.dict` beside `genotyped_ref_fasta`. **Older CRAMs** need one beside `aligned_ref_fasta`.

## Resources (`conf/base.config`)

These are starting points, from cpg-flow, to be tuned on real runs. As in nf-align, every process names one
exact `machineType` per attempt.

| process                           | CPUs           | memory          | disk (genome / exome)              | machine type, by attempt                          |
|-----------------------------------|----------------|-----------------|------------------------------------|---------------------------------------------------|
| `GENOTYPE_INTERVAL`               | 4              | 26 GB × attempt | 20 / 20 GB (the CRAM isn't staged) | `n2-highmem-4`, `n2-highmem-8`, `n2-highmem-16`   |
| `POSTPROC_GVCF`                   | 2              | 7 GB × attempt  | N × 1.5 + 2 GB / 20 GB             | `e2-highmem-2`, `e2-highmem-4`, `e2-highmem-8`    |
| `BCFTOOLS_STATS`                  | 1              | 2 GB × attempt  | 50 / 20 GB                         | `e2-standard-2`, `e2-standard-2`, `e2-standard-4` |
| `MAKE_INTERVALS`                  | 1              | 2 GB            | default                            | `e2-standard-2`                                   |
| `BUILD_SAMPLESHEET`, `REGISTER_*` | 1              | 4 / 2 GB        | default                            | `e2-standard-2`                                   |
| `PREFLIGHT`                       | native `exec:` |                 |                                    |                                                   |

## Files (planned)

```
main.nf                          param types and validation, run-records check, sample sheet, destinations and outputs, pre-flight, GENOTYPE, GVCF_QC, registration
nextflow.config                  param defaults, profiles (standard: Google Batch, docker, test)
nextflow_schema.json             launch-form schema
conf/base.config                 containers, resources, error strategy, labels, publishDir
subworkflows/local/genotype/     GENOTYPE
subworkflows/local/gvcf_qc/      GVCF_QC
modules/local/<process>/         build_samplesheet, preflight, make_intervals, genotype_interval, postproc_gvcf, variant_calling_metrics, bcftools_stats, register_analysis
bin/build_samplesheet.py         cohorts or SG IDs -> sample sheet (needs only `metamist`)
bin/freemix.sh                   .selfSM -> FREEMIX
bin/postproc_gvcf.sh             merge + cpg-flow's postprocessing
test/                            pytest: sample sheet unit tests, stub runs of main.nf, the <NON_REF> filter fixture
```

## Running

```bash
analysis-runner seqera ...   # or, directly:
nextflow run main.nf --dataset <analysis ds> --sequencing_type genome --arGuid <guid> --cohorts COH1,COH2

nextflow run main.nf -resume ...                       # recovery
nextflow run main.nf -profile test -stub -c sheet.config ...   # offline wiring check
uv run pytest
```

## Status and next steps

**Not started.** The design above was settled in a grilling session on 2026-11-12. Before the first production
run:

1. Create the `gvcf_qc` analysis type in Metamist.
2. Convert cpg-flow's genome evaluation interval list to BED, and check there's a `.dict` beside the
   DRAGEN-masked FASTA.
3. Test reading CRAMs from GCS on Google Batch, with one CRAM on each reference (RISKS #3).
4. Benchmark a GIAB sample with hap.py: bwa-mem3 plus this genotyping, against truth and against the
   DRAGEN-mode baseline (ADR 0001).
5. Separately: back-fill cpg-flow's gVCF QC results as `gvcf_qc` analyses.

**Judgement calls made while writing this, for review:**
- Exclusions are subtracted even when N is 1, unlike cpg-flow.
- A genotyped SG needs no CRAM or `cram_qc`, so `--qc_genotyped` works on SGs whose CRAMs were archived.
- There is no switch to turn gVCF QC off, unlike nf-align's `--cram_qc`.
- The GATK version in the meta is taken at runtime from the image, not from a param.
- The resource numbers.
