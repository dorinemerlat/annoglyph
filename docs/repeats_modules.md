# Repeat annotation modules

Repeat annotation now consumes `[meta, fasta]` and emits `[meta, masked_fasta]`
and `[meta, repeat_gff]`. Local modules retain their existing biological commands,
using `meta.id` instead of a separate identifier in tuples.

| Previous process | nf-core module | Preserved options |
| --- | --- | --- |
| REPEATMODELER (database step) | repeatmodeler/builddatabase | `-name <id> -engine ncbi` |
| REPEATMODELER (analysis step) | repeatmodeler/repeatmodeler | `-engine ncbi -LTRStruct -threads 30` |
| CD_HIT_EST | cdhit/cdhitest | `-c 0.8 -n 10 -l 80 -aS 0.8 -g 1 -sc 1 -sf 1 -d 0 -M 0 -T 10` |
| REPEATMASKER | repeatmasker/repeatmasker | `-e ncbi -lib <library> -a -gccalc -norna -excln -s -xsmall -gff -pa 30` |

Container overrides retain RepeatModeler 2.0.5, CD-HIT 4.8.1 and RepeatMasker 4.2.1.
The custom unmasking, HMMER database preparation/scans, MCHelper, classification,
GFF conversion and repeat landscape steps remain local.

The local patches recorded in `modules.json` provide:

- BuildDatabase arguments and a captured command log.
- File-based version outputs for CD-HIT and RepeatMasker (`storeDir` compatibility).
- CD-HIT's explicit `-M 0`, instead of the upstream memory-derived limit.
- RepeatMasker's genome and library in the same input tuple, so per-genome libraries
  stay associated with the correct genome even if upstream tasks finish out of order.
- RepeatMasker `.align`, `.cat*` and log outputs, and preservation of the input FASTA
  extension when constructing the tool's output names.

RepeatModeler uses its native log rather than the old captured stdout log; its families
FASTA keeps the original `<id>-families.fa` name. RepeatMasker names its GFF
`genome_<id>.fa.gff`; EXOGAP's classification-aware GFF conversion still emits
`<id>_repeats.gff`. `.cat*` is optional because availability depends on the tool run.

Both pooled and per-genome library modes have a stub integration test. These tests
validate channel association, output cardinality and metadata and use separate caches
so simulated outputs cannot populate the real analysis cache:

```bash
nextflow run test_repeats_channels.nf -profile singularity \
    -c test/validation/repeats_channels.config -stub-run
nextflow run test_repeats_channels.nf -profile singularity \
    -c test/validation/repeats_channels.config -stub-run \
    --group_consensus_sequences false --cacheDir test/validation/cache_individual
```

The test expects the `test_full2` assemblies and downloaded BUSCO dataset to exist.
Run the first command again without `-resume` to verify stored-output reuse.

Real module smoke test: CD-HIT-EST and RepeatMasker completed with the preserved
options and existing containers; RepeatMasker produced non-empty alignments. A second
run without `-resume` reused both module outputs and their version files.

For the IFB full fungal test, `test/validation/slurm_runtime.config` excludes the three
nodes found to lack Singularity:

```bash
nextflow run main.nf -profile singularity,test_full2 -c conf/ifb.config \
    -c test/validation/slurm_runtime.config
```

The full `test_full2` preprocessing + repeat annotation run completed successfully
on Slurm with this exclusion configuration. The four nf-core repeat processes ran
on the real assemblies/libraries. Existing results were reused for the custom local
steps. Both repeat alignment files are non-empty and all seven nf-core modules
(preprocessing plus repeats) have version files in their `storeDir` caches.
A final full run without `-resume` completed successfully, reusing every process
from `storeDir` without submitting new Slurm jobs.
