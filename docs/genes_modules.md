# Protein-coding gene modules

The active gene workflow uses `[meta, files…]`. `meta.id` is the sample identifier,
`meta.rnaseq_aligner` identifies STAR/VARUS evidence, and
`meta.annotation_method` distinguishes `braker2`, `galba`, `braker-with-star` and
`braker-with-varus`. Metadata is copied when adding these fields. Genome/BAM and
genome/GFF associations use the sample identifier; functional joins use the full
sample/method metadata. A shared evidence channel is broadcast with `.first()`.

## nf-core replacements

All modules are pinned in `modules.json` to commit
`b2bcd8312599aacd849e69721e20ac8138018f83` of nf-core/modules. Adjacent `.diff`
files record the changes needed by this pipeline.

| Module | Preserved behavior | Adaptation |
| --- | --- | --- |
| `samtools/merge` | Original VARUS image; merge without indexing, 4 threads | File-based versions; retain original thread count |
| `agat/spkeeplongestisoform` | AGAT 1.6.1; longest isoform per gene | Original canonical GFF filename; omit the upstream additional `--cpu`; file-based versions |
| `agat/spextractsequences` | AGAT 1.6.1; `--protein`; original protein FASTA filename | GFF and corresponding genome travel in one tuple; file-based versions |
| `interproscan` | 5.59-91.0; 10 CPUs, 8 GB; `-f TSV -dp -goterms -iprlookup -pa -t p` | Original prefix; remove stop symbols from an input copy; complete IFB image data; file-based versions |
| `eggnogmapper` | 2.1.12; 8 CPUs, 16 GB; `--report_orthologs`; `/tmp` scratch/temp; bank `/shared/bank/emapperdb/5.0.2/` | Preserve output subdirectory and all four outputs; disable automatic upstream `--dbmem`; file-based versions |

The eggNOG database is staged as a symlink, avoiding a copy of the shared bank.
The explicit `-m diamond` matches eggNOG-mapper's default search mode.
On IFB, InterProScan uses the installed
`/shared/software/singularity/images/interproscan-5.59-91.0.sif`, containing the
complete 5.59-91.0 database (46 GB). The former BioContainers image contains
only test models and unprepared HMM indices. On another environment, set
`params.interproscan_db` to a complete matching data directory. The BioContainers
image requires external data; the wrapper rejects its built-in test database
for real annotation. The original 18 analyses are selected explicitly with `-appl`, because the
complete IFB image also provides optional Phobius, SignalP and TMHMM analyses
that were absent from the original image. Those additional analyses are not
automatically enabled.

The InterProScan image includes BusyBox `ps`, which lacks the `-p` and column
syntax used by its Java process monitor. A small adapter is enabled only when
that command fails; it preserves the ten configured CPUs and all analyses.

## Local modules retained

- BRAKER2 and BRAKER3: species name, fungal mode, explicit tool paths, writable
  AUGUSTUS configuration, and five-attempt fallback to an empty GFF are retained.
  The nf-core BRAKER3 wrapper does not implement that fallback and requires
  additional outputs even when prediction fails.
- GALBA: no corresponding nf-core module at the pinned commit.
- STAR single/paired: retain the combined genome indexing/alignment commands,
  coordinate-sorted BAM, `intronMotif`, gzip detection and memory limits.
  The nf-core STAR submodules add outputs/default options and need a distinct
  index channel. Multiple paired libraries now use every R1/R2 pair; the previous
  shell code considered only `input1` and `input2`.
- SRA/runlist and VARUS: preserve acquisition and selection logic. VARUS creates
  its empty fallback BAM from a SAM dictionary; `/dev/null` was not a valid SAM.
- Gene-length filtering: retain `agat_sp_filter_gene_by_length.pl --test '>='
  -s 90`. nf-core's `spfilterbyorfsize` invokes a different script and is not an
  exact replacement.
- BLASTP: preserve outfmt 11 archive and the exact 24-column `blast_formatter`
  TSV, `1e-6`, word size 3 and max HSPs 20. The nf-core module does not emit the
  archive. Add the column header required by the existing annotation-transfer
  parser (the former module wrote a headerless TSV). Correct the image registry to quay.io, retaining version 2.17.0.
- GFF rewriting, evidence preparation, annotation transfer, ToxProt and SignalP
  retain the pipeline-specific behavior. SignalP's summary filename now matches
  its declared sample/method output. BLAST now honors `params.blast_db`; the IFB
  configuration forwards the installed Swiss-Prot bank through `BLASTDB`.

## Persistent cache

The existing global `storeDir` rule remains active. Annotation modules include
`annotation_method/id` in their tags; BRAKER3 and mapped-read counts include
`rnaseq_aligner/id`, preventing collisions. SRA downloads retain their accession
cache key, allowing existing downloads to be reused. nf-core version outputs are
`versions.yml` files rather than `eval` values, so they can be restored from
`storeDir` without `-resume`. Existing local caches are reused where names match.
Taxonomy JSON is parsed with `JsonSlurperClassic`, avoiding mutable lazy maps
shared across concurrent joins.

The main workflow now activates the gene block when `params.module_genes` is
true. The already-commented QC/contamination, comparison/selection and
complementary-analysis blocks remain disabled and have not been migrated in this
step. Active outputs expose all method-specific functional GFFs and canonical
proteins; no best-method selection is performed.

## Validation

The root-level `test_genes_channels.nf` harness exercises two samples, local
single/paired reads, SRA stubs, BAM merging, STAR/VARUS evidence, BRAKER2/BRAKER3/
GALBA and the functional joins. It asserts output cardinality and sample/method
filenames. `test/validation/genes_channels.config` uses a separate cache; stub
outputs never populate the full dataset cache.

```bash
nextflow run test_genes_channels.nf -profile singularity \
  -c test/validation/genes_channels.config -stub-run
```

The full dataset test is run with the user-provided IFB environment and
`-profile singularity,test_full2 -c conf/ifb.config`, plus
`test/validation/slurm_runtime.config` to exclude the three nodes without a
Singularity executable. Logs reside in `test/validation/genes_*.log`.

Channel tests passed for STAR, VARUS and both aligners, including the exactly
one-million mapped-read rejection boundary and SignalP/ToxProt stub branches.
A repeat run without `-resume` restored all outputs through `storeDir`.
The real SAMtools merge smoke test retained both input reads; real AGAT outputs
for the two full-test BRAKER2 annotations matched both the legacy canonical
GFFs and protein FASTAs byte for byte.

A real InterProScan smoke test with explicitly prepared image test data
completed successfully, producing an empty TSV (no matching model in that
small fixture bank). This validates execution only, not biological annotation
with a full database. The IFB complete-image functional test is separate.

The complete IFB-image functional smoke test passed with the original 18
InterProScan analyses, real BLAST/eggNOG and functional annotation transfer.

The yeast test profiles now share `conf/test_proteins.config`: 2,000 reviewed,
nonredundant Saccharomyces proteins (~1 MB), with provenance and regeneration
in `test/data/README.md`. Real miniprot runs produced 1,994 and 1,963 alignment
models for the two full2 genomes. Gene predictions use the separate
`genes_annotation/test_proteins_2000` cache; production defaults are unchanged.
The previous 13 GB evidence run was stopped before relaunching this test.

`test/validation/continue_genes_full2.sh` runs the corrected full pipeline,
then repeats it without `-resume` to validate persistent cache reuse. Its logs
are `genes_full2_fixed.log` / `genes_full2_store.log`; `.exit` files record the
exit status. The bounded full test is still running: GALBA's pygustus home
write failure required binding its task directory as container home, while
BRAKER3 currently retries GeneMark-ETP failures. Empty GALBA outputs from the
home failure were preserved outside the active cache to allow a corrected
rerun. Full integration success is not yet established.

## Automatic selection and final export

QC/OMArk filtering and automatic comparison/selection have now been enabled
and migrated to metadata-first tuples. Gene workflow `gff`/`proteins` outputs
now contain the selected method; `method_gff`/`method_proteins` retain all
method-specific outputs. See [ncrna_and_export.md](ncrna_and_export.md) for
the selection rules, noncoding priority filter, final bundle and submission
limits. Earlier validation results above concern the preceding migration;
these newly enabled stages have only been reviewed statically, not executed.
