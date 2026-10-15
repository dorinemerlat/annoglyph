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
per species, excluding runs above 1 billion bases. These are whole-run filters,
not read subsampling. Full mode uses the complete prepared protein pool and
STAR with both SRA limits disabled; technical eligibility filters still apply.
No eligible SRA runs means no downloaded RNA-seq evidence for that species.
The supplied samplesheets contain no local RNA-seq; local FASTQs added later
are not subsampled by these parameters. Neither mode invokes VARUS.

All three annotation categories remain enabled and genome sequences are never
truncated. The larger plant/insect genomes can still take substantial time.
SRA searches depend on current public availability, with metadata recorded by
the workflow; full mode does not mean all protein databases or every public
RNA-seq dataset. The legacy group profile names select minimal mode.

## Launch

See the [main README](../README.md#integration-tests) for environment setup and
all six Nextflow commands. For example, from the project root on IFB:

```bash
nextflow run main.nf -profile singularity,test_hanseniaspora_minimal -c conf/ifb.config -c test/validation/slurm_runtime.config
nextflow run main.nf -profile singularity,test_hanseniaspora_full -c conf/ifb.config -c test/validation/slurm_runtime.config
```

Or submit the launcher, which loads the modules/environment:

```bash
sbatch nextflow.sbatch test_hanseniaspora_minimal
sbatch nextflow.sbatch test_pieris_full
sbatch nextflow.sbatch test_utricularia_minimal -resume
```

The launcher defaults to `test_hanseniaspora_minimal`. It uses the IFB account
configured in `nextflow.sbatch`; `slurm_runtime.config` excludes four nodes
previously found without Singularity. Adapt execution settings to your site.
Results and persistent caches are separated by mode under
`test/test_<group>/<minimal|full>/{results,cache}`. Existing older output
directories are retained. Add `-resume` to resume the same mode.

Bank paths and main program options are defined in `nextflow.config`.
No complete pipeline run is implied by these instructions; execution is left
to the user.
