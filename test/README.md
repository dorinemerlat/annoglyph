# NCBI integration datasets

The three test groups cover 11 assemblies: Hanseniaspora (5 species), Pieris
(3 species) and Utricularia (3 species). Accessions and taxon IDs are recorded
in [assemblies.json](assemblies.json). See the [assembly statistics](assembly_stats.md)
([TSV](assembly_stats.tsv)) computed from the installed FASTAs.

## Installation and preparation

For a fresh installation:

```bash
module load ncbi-datasets-cli/18.14.0
bash test/create_test_data.sh
```

**The download script recreates each dataset directory**, including any results
and cache below it. To regenerate evidence, samplesheets and statistics from
already downloaded assemblies while retaining results/cache:

```bash
python3 test/prepare_test_sets.py
```

Preparation creates, for each group:

- `samplesheet_<group>.csv`: all assemblies, used in full in both modes.
- `test_<group>/protein_evidence_2000.fasta`: deterministic round-robin selection
  of 2,000 distinct proteins of 90–5,000 amino acids.
- `test_<group>/protein_evidence_full.fasta`: all distinct usable downloaded
  proteins without a length cutoff.
- Adjacent protein TSV files: assembly, original header, length and SHA256.
- `assembly_stats.md` and `assembly_stats.tsv`: sizes, sequence counts, N50,
  L50, longest sequence, GC (ACGT denominator), ambiguous bases and accessions.

Both evidence modes remove a trailing stop and deduplicate by sequence.
Empty sequences and sequences containing unsupported characters or internal
stops are excluded. U. reniformis has no protein annotation in its package;
its evidence therefore comes from the other Utricularia assemblies. Evidence
is derived from the tested genomes: these are integration fixtures rather
than independent annotation benchmarks.

## Minimal versus full

| Group | Minimal profile | Full profile |
| --- | --- | --- |
| Hanseniaspora | `test_hanseniaspora_minimal` | `test_hanseniaspora_full` |
| Pieris | `test_pieris_minimal` | `test_pieris_full` |
| Utricularia | `test_utricularia_minimal` | `test_utricularia_full` |

Minimal mode uses 2,000 proteins and STAR with at most three eligible SRA runs
per species, excluding runs above 1 billion bases. Each selected run is capped
at its first 500,000 spots (one pair counts as one spot). Evaluation is disabled
for minimal profiles. Full profiles enable `evaluation_mode` and download
OrthoDB/UniProtKB evidence from a shared clade above genus, then exclude the
target genus separately for each genome. The assembly-derived full protein
fixtures are preserved but are not used by these profiles. Full mode uses all
reads from all eligible SRA runs; technical eligibility filters still apply.

Evaluation requires an ETE NCBI taxonomy SQLite database configured through
`protein_evidence_taxonomy_db` (or the existing `omark_taxonomy_db`). Unknown
protein taxonomy is excluded and counted in per-genome JSON provenance reports.
An empty filtered set is an error. Downloads can be substantial. See
[genus-exclusion evaluation](../README.md#genus-exclusion-evaluation) for details.
If you override `protein_set`, that user-supplied bank is used unchanged,
including in full profiles.
No eligible SRA runs means no downloaded RNA-seq evidence for that species.

BRAKER3 uses a task-local copy of GeneMark-ETP with an exon-count compatibility
fix. The bundled `super.pl` writes intron counts in training summaries, while
`gmes/parse_set.pl` expects exon counts. In intron-poor genomes this can produce
negative transition probabilities and abort masking-penalty estimation. The
fix reports the number of CDS exons, without changing the training genes or
prediction thresholds. `test_patch_genemark_super.py` checks single-exon and
multi-exon summaries; an unrecognized upstream summary format fails explicitly.
The original tools are documented in the
[GeneMark-ETP source](https://github.com/gatech-genemark/GeneMark-ETP).

The supplied samplesheets contain no local RNA-seq; local FASTQs added later
are not subsampled by these parameters. Neither mode invokes VARUS.

Minimal profiles require more than 100,000 mapped alignments per RNA-seq BAM
to enable RNA-supported prediction (`rnaseq_min_mapped_reads`). The default
threshold for full profiles and production runs remains 1,000,000.

All three annotation categories remain enabled and genome sequences are never
truncated. The larger plant/insect genomes can still take substantial time.
SRA searches depend on current public availability, with metadata recorded by
the workflow; full mode does not mean all protein databases or every public
RNA-seq dataset. The legacy group profile names select minimal mode.

## Launch

See the [main README](../README.md#integration-tests) for environment setup and
all six Nextflow commands. For example, from the project root:

```bash
nextflow run main.nf -profile singularity,test_hanseniaspora_minimal
nextflow run main.nf -profile singularity,test_hanseniaspora_full
```

Add your site configuration with `-c conf/local.config` when needed.
Machine-specific configurations and batch launchers are kept locally and are not
versioned.

Results and persistent caches are separated by mode under
`test/test_<group>/<minimal|full>/{results,cache}`. Existing older output
directories are retained. Add `-resume` to resume the same mode.

General program options are defined in `nextflow.config`. Put local database
paths and execution settings in the ignored `conf/local.config`, which is loaded
automatically when present.
No complete pipeline run is implied by these instructions; execution is left
to the user.

## RNA-seq preprocessing and per-run alignment checks

Portable synthetic inputs live in `test/data/rnaseq`. Run these checks from
separate launch directories, using absolute paths to the source tree; all
generated work, caches and results stay under the launch directory. Load
Nextflow 25.10.2+ and Singularity first and configure an image cache, for example:

```bash
export NXF_SINGULARITY_CACHEDIR=/path/to/annoglyph/.singularity
nextflow run /path/to/annoglyph/test_rnaseq_utils.nf
nextflow run /path/to/annoglyph/test_rnaseq_runs.nf -profile singularity \
    -c /path/to/annoglyph/test/validation/rnaseq_runs.config -stub-run
```

The channel check covers two genome IDs, a shared SRA accession, single/paired
libraries and colliding local filenames. It expects eight per-library BAMs and
two merged per-genome BAMs with the original metadata restored. SRA acquisition
is stubbed, so it makes no NCBI requests. Stub outputs are not real annotations.

The following execute real tools on tiny inputs:

```bash
nextflow run /path/to/annoglyph/test_rnaseq_tools.nf -profile singularity \
    -c /path/to/annoglyph/test/validation/rnaseq_runs.config
nextflow run /path/to/annoglyph/test_rnaseq_qc.nf -profile singularity \
    -c /path/to/annoglyph/test/validation/rnaseq_runs.config
```

The first verifies fastp, STAR, separate sorting and merging, then checks that
the resulting BAM contains 90 mapped reads in coordinate order. The QC check
verifies default filtering and mate synchronization. Repeat the QC check with
`--rnaseq_clean_reads false --test_qc_expected_mixed 4` to retain every read, or
with `--rnaseq_min_read_length 80 --test_qc_expected_good 0
--test_qc_expected_mixed 0` to test valid empty outputs. Changing QC settings
uses a separate persistent cache namespace; existing raw SRA downloads remain
reusable. These tests do not establish annotation accuracy on full genomes.
For additional channel checks, run `test_rnaseq_runs.nf -stub-run` with
`--test_read_count 5000000000` (64-bit mapped counts),
`--test_clean_read_count 0` (no alignment for empty libraries), or
`--test_local_libraries false` (one run per genome, no redundant merge).
Use a separate launch/cache directory for each stub scenario.

## Disk-space failure recovery

The default error strategy retries explicit `No space left on device` failures,
including scratch allocation errors recorded only in `.command.log` before
the tool starts. It also checks `.command.err` and `.command.out`. Other error
handling remains unchanged; the existing `maxRetries = 4` allows at most five
attempts per task. Persistent storage failures still fail the run.

Slurm schedules each retry normally: a retry may use another node, but does not
automatically exclude the previous node. Keep large downloads on shared storage
and limit concurrent downloads in the site configuration when local `/tmp`
capacity is insufficient. Configuration changes apply to subsequent launches;
they do not change an already running Nextflow session.

To validate recovery locally, run from a separate launch directory without a
Slurm/site profile:

```bash
nextflow -c /path/to/annoglyph/conf/modules.config run /path/to/annoglyph/test/test_disk_space_retry.nf
```

This simulates wrapper, stdout and stderr failures, each succeeding on its
second attempt. Add `--failure_mode ordinary` to verify that an unrelated tool
error is not retried, or `--failure_mode persistent` to verify that retries stop
at the configured limit. Both negative checks must exit unsuccessfully.

## Execution statistics

Use `--run_group` to give an experiment and its resumes the same identifier. Test profiles set it automatically: `test_minimal` for `test`,
`test_<genus>_minimal` or `test_<genus>_full` for the three genus datasets,
and `test_full1` / `test_full2` for the legacy profiles. Short genus aliases use
the same group as their minimal profile. Override with `--run_group my_experiment`
to distinguish independent experiments. These defaults group reports and traces;
the Nextflow execution name remains automatic unless supplied with `-name`.
After the run, collect its traces with `bin/summarize_execution.py`; see
[execution statistics](../README.md#execution-statistics-across-resumed-runs)
for the tables, Slurm accounting, plot selection modes and requested-resource
comparison. The earlier minimal validation traces use `minimal_validation`.
