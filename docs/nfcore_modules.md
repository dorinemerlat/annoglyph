# nf-core modules in AnnoGlyph

The active workflow uses nf-core modules for SeqKit statistics, BUSCO, RepeatModeler,
RepeatMasker, CD-HIT-EST, SAMtools merge, AGAT canonical isoform/sequence extraction,
InterProScan, eggNOG-mapper, and the tools below. Exact upstream revisions and
compatibility patches are recorded in `modules.json`. Runtime containers and
scientific options remain configurable through `conf/modules.config`,
`params.program_args` and `params.program_resources`.

| Workflow operation | nf-core module | Preserved behavior |
| --- | --- | --- |
| Ribosomal RNA prediction | `barrnap` | Separate `euk` and `mito` calls; `--quiet`; metadata and GFF outputs. |
| Transfer RNA prediction | `trnascanse` | Eukaryotic `-E`, progress, GFF, sequences, statistics and secondary structures. |
| Rfam searches | `infernal/cmscan` | Assembly-specific `-Z`, `--cut_ga --rfam --nohmmonly --oskip`, clan handling and format 2. |
| Protein completeness | `busco/busco` | Protein mode, local shared lineage, `--offline --force`, JSON/text summaries before/after contamination filtering. |
| RNA-seq alignment | `star/genomegenerate`, `star/align` | One index per genome, single/paired reads, compressed inputs, adaptive index size, `intronMotif`, coordinate-sorted BAMs and logs. |
| SRA retrieval | `custom/sratoolsncbisettings`, `sratools/prefetch`, `sratools/fasterqdump` | Explicit NCBI settings, `--max-size 100G`, `--split-files`, gzip compression, genome/run associations. |

BRAKER2, BRAKER3 and GALBA retain their local implementations and existing failure
handling. BLAST and AGAT statistics retain their custom output handling.

## Adapters and compatibility patches

- Infernal's nf-core output is compressed. The small local
  `NORMALIZE_INFERNAL_OUTPUTS` adapter restores the uncompressed tables consumed
  by the existing ncRNA reconciliation, without changing hit filtering.
- The tRNAscan-SE patch exports `.struct` files, correctly declares compressed
  FASTA output, and allows the log to be absent with `--progress`. Temporary
  files use the task directory. GFF output retains the prefix expected downstream.
- New modules record `versions.yml` for the pipeline's existing provenance and
  persistent-cache handling, matching the convention of modules already installed.
  Template-based SRA modules also record versions after successful execution.
- Previously installed BUSCO includes its existing offline/staging compatibility
  patch. Protein-mode calls stage the shared download root, preserving lineage names
  such as `saccharomycetes_odb12`. A failed protein-mode batch is an error,
  rather than an apparently successful task without JSON summaries.
- SRA accessions are deduplicated for retrieval; run layout is taken from the
  runlist. Reads are associated with each genome before optional subsampling.
  The minimal-test cap remains 500,000 spots per run. Full profiles retain all reads.
- STAR indexes are associated by genome ID and reused for single/paired alignment.
  Read layout metadata is removed before joining BAMs back to the genome metadata.
- Migrated tools use a dedicated `nfcore_tools_v2` persistent-cache namespace,
  keeping old module outputs intact. This isolates the implementation migration;
  it does not implement general cache invalidation based on program options.

Legacy local tool files remain available for reference but are no longer imported
by the active workflow for the operations listed in the table.
