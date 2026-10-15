# Preprocessing modules

Custom modules live in `modules/local/`, preserving their original directory hierarchy.
Preprocessing and repeat annotation channels have been migrated. Gene annotation
modules retain their existing tuple signatures; their include paths have been updated.

Preprocessing genome channels now emit `[meta, fasta]`, with `meta.id` holding the
identifier derived from the samplesheet name. FASTA and mapping outputs from
`REFORMAT_GENOME` retain this metadata. Joins temporarily use `meta.id` as a key,
then remove that separate key before passing data to a module.

The following official modules are vendored under `modules/nf-core/`:

| Module | Configuration |
| --- | --- |
| `seqkit/stats` | `--tabular` is built in; `ext.args` supplies `--seq-type dna --all`. |
| `busco/download` | Downloads the selected shared lineage into `busco_downloads/lineages/`. |
| `busco/busco` | Mode `genome`, selected shared lineage, local download root, `--offline --force`, 10 CPUs, no intermediate cleanup. |

`modules.json` records the upstream revision. Local patches replace `eval()` version outputs with `versions.yml` files to support
`storeDir`; patch paths are recorded in `modules.json`.
Process configuration overrides their containers to retain AnnoGlyph's existing SeqKit
2.8.2 and BUSCO 6.0.0 versions.

`NORMALIZE_SEQKIT_STATS` restores AnnoGlyph's `genome_id` column before assembly-size
parsing and table merging. The merged table remains `all_seqkit.stats`.

BUSCO uses the upstream module's batch layout and file names. Assembly summaries
are emitted through `busco_assembly_results` as `[meta, json]`. The shared dataset
directory channel emits `[[id: lineage], dataset_directory]`. Summary files may be
absent when BUSCO finds no genes; batch failure reports are available from the
upstream module. Their names and log contents differ from the old custom module.

The three patched modules use AnnoGlyph's global `storeDir` cache, including their
`versions.yml` files. A run without `-resume` can reuse these stored outputs.
BUSCO uses `stageOutMode = 'copy'` because its directory and nested file outputs
overlap; moving the directory first would remove files still awaiting staging.
The supplied lineage root is staged as `input_lineages` so it cannot be mistaken
for a newly downloaded lineage output and copied into every per-genome cache.
The BUSCO patch also handles an unset Augustus configuration variable in AnnoGlyph's
existing container, which uses the default Miniprot predictor for genome mode.

The main workflow runs preprocessing and repeat annotation. Gene annotation calls
remain commented out pending their separate channel migration.

Validation: the real preprocessing run completed on both `test_full2` fungal
assemblies with the existing Singularity images and the local executor. A subsequent
run without `-resume` reused all three nf-core modules from `storeDir`. Slurm testing
also revealed that cpu-node-128, cpu-node-129 and cpu-node-130 lack the container runtime; the
repeat annotation test excludes those three nodes with a temporary test configuration.
