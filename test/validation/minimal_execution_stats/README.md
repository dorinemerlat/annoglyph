# Execution statistics: minimal_validation

Includes actual executions across selected traces, including failed attempts and recomputations.

| Metric | Value |
| --- | ---: |
| Trace files | 14 |
| Task executions | 118 |
| Completed tasks | 82 |
| Failed tasks | 24 |
| Aborted tasks | 12 |
| Consumed CPU hours | 33.47 |
| Failed/aborted CPU hours | 4.00 |
| Observed span, hours (includes resume gaps) | 6.59 |
| Observed active hours (union of task intervals) | 4.91 |
| Maximum single-task RSS, GiB | 11.40 |

CPU coverage: 118/118 executions; 118 from Slurm, 0 estimated. Runtime coverage: 106/118.

[Global table](execution_summary.tsv) · [Process table](process_summary.tsv) · [Individual tasks](tasks.tsv)

These timings are not a reconstruction of an uninterrupted fresh run. Processes skipped through storeDir require their original traces. Controller resources are excluded. See [measurement notes](README.txt).

## Processes

| Process | Executions | CPU hours | Max runtime, min | Peak RSS, GiB | Max requested time, h | Max requested CPUs | Max requested memory, GiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| AnnoGlyph:GENES_ANNOTATION:FUNCTIONAL_ANNOTATION:EGGNOGMAPPER | 3 | 19.44 | 166.90 | 6.60 | 20.00 | 8.00 | 16.00 |
| AnnoGlyph:GENES_ANNOTATION:FUNCTIONAL_ANNOTATION:INTERPROSCAN | 4 | 5.40 | 39.17 | 6.10 | 4.00 | 10.00 | 8.00 |
| AnnoGlyph:GENES_ANNOTATION:RNASEQ_ASSEMBLY:STAR_SINGLE | 3 | 1.89 | 14.30 | 1.20 | 4.00 | 8.00 | 16.00 |
| AnnoGlyph:GENES_ANNOTATION:RNASEQ_ASSEMBLY:STAR_PAIRED | 1 | 1.78 | 17.43 | 0.94 | 4.00 | 8.00 | 16.00 |
| AnnoGlyph:GENES_ANNOTATION:FUNCTIONAL_ANNOTATION:BLASTP | 3 | 1.09 | 4.05 | 0.50 | 4.00 | 10.00 | 2.00 |
| AnnoGlyph:GENES_ANNOTATION:STRUCTURAL_ANNOTATION:BRAKER2 | 7 | 0.75 | 3.25 | 0.40 | 4.00 | 8.00 | 16.00 |
| AnnoGlyph:GENES_ANNOTATION:STRUCTURAL_ANNOTATION:GALBA | 6 | 0.70 | 12.88 | 0.92 | 4.00 | 8.00 | 16.00 |
| AnnoGlyph:PREPROCESSING:BUSCO_BUSCO | 1 | 0.59 | 4.13 | 2.20 | 4.00 | 10.00 | 10.00 |
| AnnoGlyph:NCGENES_ANNOTATION:INFERNAL | 1 | 0.50 | 1.68 | 2.50 | 120.00 | 30.00 | 16.00 |
| AnnoGlyph:GENES_ANNOTATION:ANNOTATION_QC_AND_CONTAMINATION_FILTERING:OMARK | 4 | 0.45 | 12.22 | 11.40 | 4.00 | 8.00 | 32.00 |
| AnnoGlyph:GENES_ANNOTATION:RNASEQ_ASSEMBLY:DOWNLOAD_SRA | 3 | 0.19 | 1.78 | 1.10 | 4.00 | 10.00 | 10.00 |
| AnnoGlyph:GENES_ANNOTATION:STRUCTURAL_ANNOTATION:BRAKER3 | 5 | 0.16 | 2.12 | 0.12 | 4.00 | 8.00 | 16.00 |
| AnnoGlyph:REPEATS_ANNOTATION:REPEATMODELER_REPEATMODELER | 1 | 0.15 | 2.95 | 0.20 | 4.00 | 8.00 | 16.00 |
| AnnoGlyph:GENES_ANNOTATION:ANNOTATION_QC_AND_CONTAMINATION_FILTERING:BUSCO_PROTEINS | 7 | 0.10 | 1.22 | 1.00 | 4.00 | 10.00 | 50.00 |
| AnnoGlyph:GENES_ANNOTATION:RNASEQ_ASSEMBLY:SUBSAMPLE_SRA_READS | 5 | 0.09 | 1.57 | 0.04 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:GENES_ANNOTATION:ANNOTATION_QC_AND_CONTAMINATION_FILTERING:DOWNLOAD_OMAMER_DB | 2 | 0.04 | 5.25 | 0.04 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:REPEATS_ANNOTATION:REPEATMASKER_REPEATMASKER | 1 | 0.02 | 0.57 | 0.23 | 4.00 | 8.00 | 16.00 |
| AnnoGlyph:REPEATS_ANNOTATION:MCHELPER | 1 | 0.02 | 1.02 | 1.10 | 4.00 | 8.00 | 16.00 |
| AnnoGlyph:COLLECT_PROVENANCE | 2 | 0.01 | 2.12 | 0.33 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:GENES_ANNOTATION:ANNOTATION_COMPARISON_AND_SELECTION:PLOT_ANNOTATION_COMPARISON | 4 | 0.01 | 0.45 | 0.25 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:POSTPROCESSING:EXPORT_ANNOTATION | 2 | 0.01 | 0.53 | 0.20 | 4.00 | 1.00 | 8.00 |
| AnnoGlyph:POSTPROCESSING:COMPARATIVE_REPORT | 2 | 0.01 | 0.39 | 0.23 | 4.00 | 1.00 | 8.00 |
| AnnoGlyph:NCGENES_ANNOTATION:TRNASCANSE | 1 | 0.01 | 0.17 | 0.22 | 4.00 | 8.00 | 2.00 |
| AnnoGlyph:REPEATS_ANNOTATION:HMMSCAN_BUSCO_GENES | 1 | 0.01 | 2.58 | 0.01 | 2.00 | 5.00 | 10.00 |
| AnnoGlyph:PREPROCESSING:BUSCO_DOWNLOAD | 2 | 0.01 | 0.47 | 0.09 | 4.00 | 1.00 | 10.00 |
| AnnoGlyph:GENES_ANNOTATION:ANNOTATION_COMPARISON_AND_SELECTION:PLOT_ANNOTATION_ANALYSIS | 2 | 0.01 | 0.46 | 0.24 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:NCGENES_ANNOTATION:NCRNA_STATS | 2 | 0.01 | 0.32 | 0.20 | 4.00 | 1.00 | 8.00 |
| AnnoGlyph:REPEATS_ANNOTATION:PLOT_REPEAT_LANDSCAPE | 3 | 0.00 | 0.42 | 0.16 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:PREPROCESSING:SELECT_BUSCO_DATASET | 1 | 0.00 | 0.26 | 0.16 | 4.00 | 10.00 | 10.00 |
| AnnoGlyph:REPEATS_ANNOTATION:HMMSCAN_RRNA | 1 | 0.00 | 0.18 | 0.01 | 2.00 | 5.00 | 10.00 |
| AnnoGlyph:GENES_ANNOTATION:ANNOTATION_COMPARISON_AND_SELECTION:GENERATE_ANNOTATION_TABLES | 2 | 0.00 | 0.17 | 0.08 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:GENES_ANNOTATION:STRUCTURAL_ANNOTATION:AGAT_FILTER_GENE_BY_LENGTH | 2 | 0.00 | 0.15 | 0.07 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:GENES_ANNOTATION:STRUCTURAL_ANNOTATION:AGAT_SPKEEPLONGESTISOFORM | 2 | 0.00 | 0.07 | 0.07 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:GENES_ANNOTATION:ANNOTATION_QC_AND_CONTAMINATION_FILTERING:AGAT_STATISTICS | 1 | 0.00 | 0.25 | 0.07 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:NCGENES_ANNOTATION:RFAM_CMPRESS | 1 | 0.00 | 0.15 | 0.01 | 120.00 | 1.00 | 16.00 |
| AnnoGlyph:PREPROCESSING:GENERATE_TAXONOMY_METADATA | 1 | 0.00 | 0.12 | 0.07 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:GENES_ANNOTATION:STRUCTURAL_ANNOTATION:AGAT_SPEXTRACTSEQUENCES | 1 | 0.00 | 0.08 | 0.13 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:REPEATS_ANNOTATION:REPEATMODELER_BUILDDATABASE | 1 | 0.00 | 0.33 | 0.03 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:NCGENES_ANNOTATION:DOWNLOAD_RFAM | 1 | 0.00 | 0.15 | 0.03 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:GENES_ANNOTATION:RNASEQ_ASSEMBLY:GENERATE_RUNLIST | 3 | 0.00 | 0.07 | 0.03 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:GENES_ANNOTATION:FUNCTIONAL_ANNOTATION:TRANSFER_FUNCTIONAL_ANNOTATIONS | 1 | 0.00 | 0.09 | 0.06 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:NCGENES_ANNOTATION:BARRNAP_MITO | 1 | 0.00 | 0.04 | 0.01 | 4.00 | 30.00 | 2.00 |
| AnnoGlyph:GENES_ANNOTATION:RNASEQ_ASSEMBLY:SAMTOOLS_MERGE | 1 | 0.00 | 0.02 | 0.01 | 6.00 | 4.00 | 10.00 |
| AnnoGlyph:GENES_ANNOTATION:ANNOTATION_COMPARISON_AND_SELECTION:SELECT_BEST_ANNOTATION | 1 | 0.00 | 0.10 | 0.07 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:REPEATS_ANNOTATION:POSTPROCESS_MCHELPER | 1 | 0.00 | 0.09 | 0.00 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:REPEATS_ANNOTATION:REFORMAT_REPEATLANDSCAPE | 1 | 0.00 | 0.07 | 0.04 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:NCGENES_ANNOTATION:BARRNAP_EUK | 1 | 0.00 | 0.04 | 0.01 | 4.00 | 30.00 | 2.00 |
| AnnoGlyph:REPEATS_ANNOTATION:CLEAN_REPEAT_FAMILIES | 1 | 0.00 | 0.09 | 0.03 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:PREPROCESSING:REFORMAT_GENOME | 1 | 0.00 | 0.05 | 0.04 | 4.00 | 1.00 | 10.00 |
| AnnoGlyph:GENES_ANNOTATION:STRUCTURAL_ANNOTATION:REFORMAT_GFF3 | 2 | 0.00 | 0.04 | 0.02 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:REPEATS_ANNOTATION:REPEATMASKER_OUT_TO_GFF | 1 | 0.00 | 0.07 | 0.04 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:REPEATS_ANNOTATION:REPEATLANDSCAPE | 1 | 0.00 | 0.08 | 0.02 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:NCGENES_ANNOTATION:MERGE_NCRNA | 1 | 0.00 | 0.02 | 0.02 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:GENES_ANNOTATION:RNASEQ_ASSEMBLY:SAMTOOLS_READ_COUNT | 1 | 0.00 | 0.01 | 0.01 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:GENES_ANNOTATION:ANNOTATION_COMPARISON_AND_SELECTION:CONCAT_ANNOTATION_TABLES | 2 | 0.00 | 0.01 | 0.00 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:PREPROCESSING:SEQKIT_STATS | 1 | 0.00 | 0.00 | 0.00 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:REPEATS_ANNOTATION:UNMASK_GENOME | 1 | 0.00 | 0.00 | 0.00 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:PREPROCESSING:NORMALIZE_SEQKIT_STATS | 1 | 0.00 | 0.00 | 0.00 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:REPEATS_ANNOTATION:MERGE_REPEATLANDSCAPE | 1 | 0.00 | 0.00 | 0.00 | 4.00 | 1.00 | 2.00 |
| AnnoGlyph:PREPROCESSING:MERGE_SEQKIT | 1 | 0.00 | 0.00 | 0.00 | 4.00 | 1.00 | 2.00 |

## Figures

When enabled: [runtime](plots/runtime.png), [consumed CPU](plots/cpu.png), [memory](plots/memory.png), [mean CPU cores](plots/cores.png). A metric with no selected measured modules has no figure; see plots_index.tsv.
