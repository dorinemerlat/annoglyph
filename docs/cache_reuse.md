# Resume the Myriapoda experiment from cache_save

`nextflow_myriapoda.config` uses `samplesheet_myriapoda_cache_save.csv`,
`cache_save`, and its validated `cache_reuse_v1/manifest.json` import manifest.
The cache retains 37 historical assemblies and their sequence mnemonics. The
active samplesheet contains 36: `strigamia-acuminata-3` is deliberately excluded
and its cached files are left untouched. The samplesheet uses
the cached reference FASTAs as inputs, so a newer NCBI assembly cannot silently
replace an assembly associated with historical annotations. The original
`samplesheet_myriapoda.csv` remains available.

From the repository directory, submit the driver with:

```bash
sbatch finalize_cache_save.sbatch
```

This finalization entry point loads Nextflow 25.10.2 and Singularity and uses
only existing annotation results. It does not run protein-coding predictors,
functional annotation, BUSCO, OMArk, Infernal, tRNAscan-SE or Barrnap. It validates
and translates the selected coding GFF sequence IDs, merges the cached ncRNAs,
regenerates figures, and builds species and comparative exports in `results/`.
Finalization products use `cache_save/finalize_cache_save_v1` separately from the
original results. The updated ncRNA merges use `cache_save/finalize_cache_save_v2`;
coding imports and exports with audited model exclusions use
`cache_save/finalize_cache_save_v3`. Global statistics and reports are rebuilt
from this run.

Official Rfam 15.1 family types cover the imported Infernal hits. Before merging,
`RESOLVE_RFAM_METADATA` checks every accession and automatically looks up missing
types in the official Rfam CURRENT metadata, cached for 24 hours. A network
failure or a genuinely unknown accession does not abort merging: the hit is
retained as `Unclassified ncRNA`, with its original accession, coordinates and
score. Cis-regulatory families are labeled `Regulatory RNA`; other official
types, including future subtypes, are preserved without an accession-by-accession
mapping. These English categories are used in GFF attributes, filtering audits,
counts and figures. `comparative/noncoding/stats/rfam_metadata_resolution.json`
records automatic resolutions, unresolved accessions, metadata provenance and
checksums. The fallback does not assign a regulatory function to unknown hits.

`outputs/cache_reuse_20261023/final_coding_sources.json` records the 36 selected
coding GFFs, historical quality flags, and source methods. Thirty-three come
from `publish_best_annotations`; the other three use the exact method already
chosen in the historical selection table and its cached functional GFF. Those
three retain the historical `no` quality flags; exporting them does not imply
that they pass current quality thresholds.

Selected historical coding GFFs contain a small number of models with fragments
on different contigs, incompatible strands, inconsistent parent spans or bounds
outside the assembly. With the user's approval these models are excluded from
the validated export, including their descendant features, instead of clipping
coordinates or predicting replacement genes. Valid models keep their original
coordinates and functional annotations. Shared features retain valid Parent
references. Empty gene containers are removed. Original cache files are never
edited. Excluded original GFF rows and the reasons appear in each species'
`protein_coding/selection/excluded_models.gff3` and `excluded_models.tsv`;
`cache_import.json` gives counts and source identities. Comparative
`protein_coding/stats/coding_import_summary.tsv` and `excluded_models.tsv`
provide the combined report. An unknown contig still fails the import because
it can indicate an incorrect assembly or ID mapping.

`nextflow_myriapoda.sbatch` runs the full pipeline instead and can calculate
missing predictions and downstream analyses as described below. Use that
entry point only when those additional analyses are intended.

## What is reused

The explicit `legacy_cache_manifest` parameter enables reuse of historical
RepeatMasker annotations, the associated masked genomes, consensus library,
family tables and landscapes, and valid raw BRAKER/GALBA predictions. Their
historical scientific results are retained; this is not a rerun with current
RepeatModeler, MCHelper or RepeatMasker arguments. The manifest is recorded as
an input in the final provenance report; unavailable software versions are not
invented.

Masked genomes are compared base by base via per-sequence SHA-256 hashes with
the cached reference genomes, ignoring softmask case and line wrapping. Only
sequence identifiers are translated. All imported GFF rows are checked for
known sequence names and coordinates within their cached assembly. Historical
RepeatMasker entries with zero length (`start = end + 1`) are excluded and
retained unchanged in `excluded_zero_length_repeats.gff` beside each import.

Complete tRNAscan-SE and Barrnap results are adapted to the current ncRNA
cache paths. Infernal tables were imported separately and use the same cached
sequence identifiers. Original files under the historical cache paths are
preserved. Corrected artifacts live under `cache_save/cache_reuse_v1/imported`;
unchanged artifacts use links.

## What continues to run

This section describes the full pipeline, not `finalize_cache_save.sbatch`.

Empty, absent or invalid raw gene predictions are not considered completed.
BRAKER/GALBA run for missing predictions. RNA-seq work is skipped for genomes
whose requested BRAKER3 aligners are already imported; each aligner is checked
individually when only one prediction is available. Existing BAM files with
old sequence headers are not reused through the new gene-analysis namespace.

Current gene filtering, protein extraction, functional annotation, QC and
selection use a separate `cache_reuse_v1` cache namespace. Current ncRNA merge,
statistics and final export run normally. Historical derived QC/selection
outputs are not relabeled as current results. Assembly BUSCO and taxonomy
outputs whose current contracts differ are generated by the current pipeline.

The standard pipeline is unchanged when `legacy_cache_manifest` is null.
To start a completely new experiment, use a separate cache and the standard
samplesheet rather than the reuse config.

## Taxonomy availability

Taxonomy generation first uses an exact sample/taxid lineage in
`cacheDir/preprocessing/get_taxonomy`. Otherwise, validated remote taxonomy and
mnemonic responses are saved immediately under `cacheDir/taxonomy_metadata/v1`,
so a failed request later in the run does not discard successful downloads.
Requests to NCBI are spaced by one second across runs sharing that cache.
HTTP 429, 408 and temporary 5xx errors, timeouts and network errors are retried
up to eight attempts with exponential delays and respect for `Retry-After`.
Permanent errors such as invalid requests still fail explicitly. Aggregate
taxonomy outputs are rebuilt for the current samplesheet instead of reusing a
table for a different set of genomes. With complete local lineages, no network
request is required. The finalization workflow uses the validated offline
taxonomy table for the 36 active genomes.

## Preparation and validation records

The preparation command is:

```bash
python3 bin/prepare_cache_reuse.py \
  --cache-dir cache_save \
  --samplesheet samplesheet_myriapoda.csv \
  --resume-samplesheet samplesheet_myriapoda_cache_save.csv --workers 4
```

Preparation refuses to replace an existing completed import manifest. It never
edits the original assemblies or predictions. The manifest records source file
identities, validation counts, available predictions and excluded repeat rows.
At workflow initialization, changed or missing source assemblies/mappings and
missing imported artifacts cause an error instead of silent reuse.

The preparation log, per-genome summary and validation outputs are retained in
`outputs/cache_reuse_20261023`. Configuration/samplesheet originals are copied
there as well. Unit tests cover sequence identity, simultaneous renaming,
coordinate bounds and audited zero-length exclusions. The Nextflow channel test
uses a small synthetic fixture with local stub commands; it verifies reuse
routing without running annotation software.
