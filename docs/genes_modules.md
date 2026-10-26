# Protein-coding gene modules

The active gene workflow uses `[meta, files…]`. `meta.id` is the sample identifier,
`meta.rnaseq_aligner` identifies STAR/VARUS evidence, and
`meta.annotation_method` distinguishes `braker2`, `galba`, `braker-with-star` and
`braker-with-varus`. Metadata is copied when adding these fields. Genome/BAM and
genome/GFF associations use the sample identifier; functional joins use the full
sample/method metadata. A shared evidence channel is broadcast with `.first()`.
Protein BUSCO task tags include the annotation method, keeping persistent caches
separate for each predictor; genome-only BUSCO tags retain the sample identifier.
The InterProScan container override is resolved at task execution, after local
database configuration has been loaded.
GALBA GTF identifiers reused on different contigs or strands are disambiguated
before GFF3 conversion. This preserves every CDS coordinate and phase while
preventing unrelated genes from being merged by AGAT. Existing GALBA-derived
caches made before this correction must be regenerated from the corrected GFF3.
Gene-structure and assembly-correlation figures use the selected method for each
genome, including only its introns; method-comparison figures retain all candidates.

## RNA-seq runs, cleaning and BAM sorting

STAR processes each SRA run or local FASTQ library independently. A paired-end
library stays an ordered R1/R2 pair. All reads passing light fastp cleaning are
retained; there is no new subsampling or tissue-selection step. The optional
existing SRA limits still apply when explicitly configured. Local FASTQs receive
the same cleaning as downloaded reads. VARUS acquisition is unchanged.

Cleaning parameters are defined in `nextflow.config` and validated by the schema:

| Parameter | Default | Meaning |
| --- | --- | --- |
| `rnaseq_clean_reads` | `true` | Enable light cleaning; false produces QC reports while retaining every read |
| `rnaseq_min_read_length` | `30` | Minimum retained length, in bases |
| `rnaseq_qualified_quality_phred` | `15` | Bases below this Phred score are unqualified |
| `rnaseq_unqualified_percent_limit` | `40` | Maximum percentage of unqualified bases per read |
| `rnaseq_max_n_bases` | `5` | Maximum N bases per read |
| `rnaseq_trim_adapters` | `true` | Detect and remove adapters |
| `rnaseq_trim_poly_g` | `false` | Optional poly-G trimming |

No deduplication, pair merging, sliding-window quality trimming or poly-X
trimming is enabled. A pair is retained only when both mates pass. Empty
libraries are reported and skipped; other libraries of that genome can proceed.
fastp 1.3.6 is pinned to the image used by the nf-core modules revision recorded
below. JSON/HTML reports are published under
`results/rnaseq/qc/<genome>/<run>/<input-fingerprint>/`.

STAR writes unsorted BAM. A separate `SORT_RNASEQ_BAM` task uses `samtools sort`
with memory bounded per thread (initially 1024 MiB, configurable with
`rnaseq_sort_memory_per_thread_mb`) and 20% task-memory headroom. Temporary sort
files stay in the shared task directory. Sorting and merging stage inputs as
symlinks rather than copying all FASTQs, BAMs and indices. Failed sorting can be
retried independently of alignment; standard Nextflow retries increase the
requested memory on exit 137/140 and walltime on exit 143.
SRA downloads and conversions also run on shared storage, with at most four
concurrent tasks per stage, avoiding large archive/conversion scratch on `/tmp`.

STAR starts with a scheduling estimate of index bytes times
`rnaseq_star_index_memory_factor` (1.5), plus 4 GiB and
`rnaseq_star_memory_per_million_reads_gb` (0.25) GiB per million individual
post-cleaning reads, with a floor of `rnaseq_star_min_memory_gb` (8). This is a
conservative initial estimate, not an exact peak-memory prediction. Existing
`program_resources.star` overrides remain available.

Sorted run BAMs are merged by genome in a deterministic run order, and the
original genome metadata is restored before gene prediction. One input BAM is
passed through without a redundant merge. Counts after merging are alignment
counts, as before; they are not a base-quality score.

Persistent QC, alignment and sort caches are separated by **genome/run/input
fingerprint**. Local-library IDs derive from their paths, and the input
fingerprint includes paths, sizes and modification times. A `rnaseq_runs_v1`
namespace incorporating cleaning settings separates these products and derived
coding annotations from historical pooled results. Downloaded SRA archives,
FASTQs remain reusable. Genome indices are shared across runs and cleaning
settings, with an ID-based cache key (historical filename-based indices may be
rebuilt once). The merge key includes the retained
run inventory. Changing cleaning parameters invalidates derived caches.

Validate channel associations with `test_rnaseq_runs.nf` and
`test/validation/rnaseq_runs.config` in a separate launch/work/cache directory.

## nf-core replacements

All modules are pinned in `modules.json` to commit
`b2bcd8312599aacd849e69721e20ac8138018f83` of nf-core/modules. Adjacent `.diff`
files record the changes needed by this pipeline.

| Module | Preserved behavior | Adaptation |
| --- | --- | --- |
| `samtools/merge` | Original VARUS image; merge without indexing, 4 threads | File-based versions; retain original thread count |
| `agat/spkeeplongestisoform` | AGAT 1.6.1; longest isoform per gene | Original canonical GFF filename; omit the upstream additional `--cpu`; file-based versions |
| `agat/spextractsequences` | AGAT 1.6.1; `--protein`; original protein FASTA filename | GFF and corresponding genome travel in one tuple; file-based versions |
| `interproscan` | 5.59-91.0; 10 CPUs, 16 GB; `-f TSV -dp -goterms -iprlookup -pa -t p` | Original prefix; remove stop symbols from an input copy; complete IFB image data; file-based versions |
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
- STAR single/paired: use nf-core indexing/alignment per run, retaining
  `intronMotif` and gzip detection. Coordinate sorting is a separate local
  samtools task; each library receives the correct genome index.
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
