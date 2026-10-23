# Automatic selection, noncoding genes and final annotation export

All active channels start with a metadata map (`meta.id`); selected coding
annotations also carry `meta.annotation_method` and the contamination filtering
status. No Nextflow execution, scientific tool run or test was performed for
this change, at the user's request. Only static checks/review were performed.

## Coding annotation selection

The gene workflow now runs AGAT statistics, protein-mode BUSCO and OMArk.
OMArk contamination candidates trigger the existing chromosome-stretch filter
using `omark_contamination_threshold`. Filtered annotations are reassessed;
clean annotations reuse their original metrics. Empty filtered proteomes are
excluded from selection. Original scoring thresholds (including the expected
15,000-gene criterion) remain unchanged. Missing QC metrics raise an error.
The existing score, ranking and BUSCO tie rules are preserved; any remaining
perfect tie is resolved by method name so exactly one method is selected.
The selection table is published in `results/selection/best_annotations.tsv`.
The selected functional GFF, rather than the raw structural GFF, feeds export.
BRAKER/GALBA still create an empty GFF after five unsuccessful attempts, as
requested. If no method survives for a genome, its export contains the other
available annotations and records the missing coding annotation in validation.tsv.

`--omamer_db /path/to/LUCA.h5` can supply a compatible existing OMA database;
otherwise the workflow downloads it. This database must match the OMAmer
version bundled with the configured OMArk container.
`--omark_taxonomy_db` supplies an ETE3 NCBI taxonomy SQLite database; an
existing `~/.etetoolkit/taxa.sqlite` is detected and staged into the container.

## Noncoding annotations

`module_ncgenes` enables tRNAscan-SE (-E), Barrnap (euk and mito searches), and
Infernal cmscan. All run on preprocessed genomes, independently of repeat
masking. Infernal retains the original GA/Rfam/nohmmonly settings and assembly
size correction. `--oclan --oskip` additionally removes lower-scoring overlaps
within the same Rfam clan, following the Rfam annotation guidance.

Rfam is downloaded and indexed once with cmpress, then broadcast to samples.
Use `--rfam_release VERSION` to pin a release (the default is CURRENT), or
`--rfam_db DIRECTORY` containing Rfam.cm, Rfam.clanin, family.txt and the four
Rfam.cm.i1* indexes. All three metadata/model files must come from the same
release. The metadata identifies tRNA/rRNA hits without guessing from model
names. A hit is suppressed only when it is of the same RNA type and overlaps
at least one nucleotide of a priority call on the same sequence and strand:

- tRNAscan-SE takes priority over Infernal for tRNA.
- Barrnap takes priority over Infernal for rRNA.
- Other Infernal RNA annotations and opposite-strand calls are retained.
- Exact duplicate Barrnap calls with the same name/coordinates/strand are merged.

The combined GFF3 and filtering audit are published under `results/ncgenes/ID/`.

## Final bundle

`results/final/ID/` contains:

- annotation.gff3: selected coding genes + repeats + noncoding annotations.
- genome.masked.fna: softmasked genome from RepeatMasker (preprocessed sequence
  when the repeat branch is disabled).
- cds.fna: concatenated CDS per transcript in biological orientation.
- proteins.faa: translated CDS, named consistently with the final GFF3.
- ncrna.fna: all ncRNA types in one FASTA, with the type in each header;
  exon/intron information is used for splicing where available.
- genome.fsa and annotation.tbl: uppercase assembly and five-column feature
  table for preparing NCBI table2asn input.
- metadata.json, validation.tsv, id_mapping.tsv and checksums.sha256.

GFF3 IDs and parent links are namespaced across sources; a locus_tag is assigned
to each gene. The default prefix is the sample mnemonic, suitable for local use;
`--locus_tag_prefix PREFIX` must be the registered prefix for an actual GenBank
submission. Use a separate run/configuration when different genomes require
different registered prefixes. CDS translation defaults to genetic code 1
(`--translation_table`), with feature-level transl_table overriding that default.
Incomplete CDS, inconsistent phases and internal stops are reported for review;
these are not silently converted into validated annotations.

These are export/preparation files, not a validated database submission.
NCBI still requires BioProject/BioSample, registered locus tags, source and
submitter metadata, a submission template, and table2asn validation. Ensembl
integration is project-specific; ENA annotated submissions require its supported
assembly/annotation format and validation. UniProt normally receives translated
proteins through nucleotide archive submissions; local FASTA IDs are not UniProt
accessions. No deposit is performed by this pipeline.

References:
- https://docs.rfam.org/en/latest/genome-annotation.html
- https://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation/
- https://www.ncbi.nlm.nih.gov/genbank/table2asn/
- https://www.uniprot.org/help/submissions

## Cache

New QC, comparison, ncRNA and postprocessing tasks use the `annotation_export_v1`
cache namespace, avoiding old incompatible output caches. Existing preprocessing,
repeat and structural/functional prediction caches are retained. As with the
existing storeDir scheme, changing input files/options while reusing a cache path
requires selecting a fresh cacheDir or deliberately invalidating affected outputs.

## Noncoding statistics and figures

The `NCRNA_STATS` Python module receives the reconciled GFF3 files, assembly
sizes from genome metadata, taxonomy.tsv and the preprocessing Newick tree.
Outputs are published in `results/ncgenes/statistics/`:

- `tables/ncrna_counts.tsv`: ncRNA types plus a Total row, genomes in columns.
- `tables/trna_isotypes.tsv`: isotype counts, genomes in columns. Unknown
  isotypes are explicit; gene/exon rows are not counted again.
- `tables/trna_anticodons.tsv`: isotype:anticodon counts, genomes in columns,
  for loci with an annotated anticodon only (U normalized to T).
- `tables/rrna_types.tsv`: 5S, 5.8S, 12S, 16S, 18S, 23S, 25S, 26S, 28S,
  SSU/LSU and Unknown counts, genomes in columns.
- `tables/genome_metrics.tsv`: assembly length, ncRNA/tRNA/rRNA counts, tRNA
  density per kb, anticodon availability and taxonomy, with one genome per row.
- `tables/ncrna_loci.tsv`: per-locus counting/classification audit.
- `tables/correlations.tsv` and `tables/group_membership.tsv`: descriptive
  Pearson/Spearman results, group sizes/status and group membership.
- Heatmaps for isotypes, isotype/anticodon and rRNA types, in PNG/PDF/SVG.
  The anticodon heatmap is omitted when no annotated anticodon is available.
- Correlation panels (PNG/PDF/SVG) for all genomes, taxonomic groups and
  distinct tree clades: tRNA count vs genome size and tRNA/kb vs genome size.
  Every point represents a genome. Individual genomes are labeled for groups
  of up to 20 genomes.

The matrix columns and heatmaps use the supplied tree's leaf order. Scientific
names are matched to taxonomy IDs/genome names through taxonomy.tsv; multiple
genomes sharing a species leaf remain adjacent. Unmatched genomes are appended
in ID order and recorded in README.txt. The preprocessing tree currently
represents taxonomy, not an inferred molecular phylogeny; heatmap branch lengths
are unit hierarchy depths, not evolutionary distances.

`--ncrna_stats_group_rank family` selects the taxonomy grouping column (e.g.
order, family or genus). `--ncrna_stats_min_group_size 3` sets the minimum number
of genomes for a correlation; groups below the minimum or with constant values
are plotted but receive NA coefficients and an explicit status. Additional tree
clades need this minimum size. Correlations are ordinary descriptive analyses,
not phylogenetically corrected models. P-values are unadjusted, and the density
correlation is mathematically coupled to genome size through its denominator.
The locus tables permit independent review of all counts and unknown labels.

Only static source review/syntax checks were performed for this module; the
workflow and scientific tools were not run, as requested.

## Public results layout (current)

The directory names described in earlier sections have been superseded by
the consolidated layout in `assets/results_README.md`, which is published as
the sole `results/README.md`. Species exports are assembled centrally under
`species/<genome_id>/`, with genome, repeats, protein_coding, noncoding and
submission subdirectories. Individual statistics and plots are generated at
export; existing assembly BUSCO/SeqKit results, selected protein QC, candidate
method tables, repeat libraries, processed family tables and landscape tables/
plots are included as staged assets. The pooled repeat library is copied into
each affected species directory, since that is the library actually used.

`comparative/` centralizes taxonomy/tree, genome and repeat summary tables/plots,
protein-method comparisons and the ncRNA matrices/heatmaps/correlations.
`pipeline_info/` now lives under the configured output directory. No per-species
README is published. Final exports use a fresh results_layout_v2 cache namespace
to avoid restoring flat-layout bundles, while costly biological analysis caches
are preserved. Legacy public directories are not deleted or moved automatically.
This reorganisation was reviewed statically only, without Nextflow or tests.

## Configurable selection thresholds

Set `params.selection_thresholds` in a Nextflow config file. Omitted keys retain
historical defaults. For example (illustrative expected gene count):

```groovy
params.selection_thresholds = [
    expected_gene_number: 6000,
    gene_number_good_deviation: 2000,
    gene_number_neutral_deviation: 4000,
    busco_complete_min: 80,
    omark_consistent_min: 70
]
```

Setting `expected_gene_number: null` disables the expected-count criterion in
both initial scoring and rank tie-breaking. The independent quality flag
`orthofinder_gene_number_min` remains configurable separately. All score
thresholds and OrthoInspector/OrthoFinder quality-flag thresholds are exposed
in nextflow_schema.json. Threshold names ending in min/max use inclusive score
boundaries; gene-count deviation at the neutral boundary receives -1, as before.
Resolved thresholds and all candidates' initial criterion scores plus selection
flags are published in `comparative/protein_coding/stats/selection_thresholds.json`
and `selection_scores.tsv`. Existing rank/BUSCO/method-name tie rules remain.
Threshold changes use distinct selection and final-export caches, without
invalidating the biological prediction/QC caches.

## Submission, provenance and tool options (current)

See `conf/programs.example.config` and `assets/submission_metadata.example.json`
for editable examples. All example accessions, prefixes and sequence names must
be replaced with real metadata. Program database/image paths are centralized in the `params` block of
`nextflow.config`; `conf/ifb.config` holds cluster execution settings. Other sites must configure their
versioned eggnog_db/BLAST/InterProScan resources explicitly. The sbatch launcher
uses the project/submission directory and environment overrides.

The single published results README now describes the local table2asn option,
per-genome/per-sequence genetic codes, Rfam/OMA defaults, UniProt release checks,
provenance files and their exact limits. Full database hashes are optional;
large-bank metadata inventories and cached-container probes are not represented
as full content/execution evidence. No cache invalidation on program option or
bank changes was added, as requested. The existing selection-specific cache
behavior from the preceding change remains unchanged. No validation tools,
Nextflow executions or tests were run while implementing these changes.
