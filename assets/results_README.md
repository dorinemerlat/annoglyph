# AnnoGlyph results

## Organisation

`species/<genome_id>/` contains the individual genome results. The identifier is
unique to an assembly; multiple assemblies from the same species have separate
directories. `comparative/` contains tables and figures combining genomes.
There is one README for the entire results directory.

```text
species/<genome_id>/
  genome/
    sequences/genome.fna                 Uppercase, reformatted assembly
    sequences/genome.softmasked.fna      RepeatMasker softmask, if enabled
    annotation.gff3                      Combined final annotation
    stats/                              Assembly statistics, mapping, BUSCO
    plots/                              Contig lengths and BUSCO summary
  repeats/
    sequences/                          Consensus library used for masking
    annotation.gff3
    stats/                              Classification, coverage, landscape
    plots/                              Composition and repeat landscape
  protein_coding/
    sequences/cds.fna
    sequences/proteins.faa
    annotation.gff3                      Automatically selected method
    stats/                              Protein lengths, feature counts, QC
    plots/                              Protein lengths, feature counts, BUSCO
    selection/selected_method.tsv
    selection/candidates/<method>/<status>/
  noncoding/
    sequences/ncrna.fna                  RNA type in each header
    annotation.gff3
    stats/                              RNA types, isotypes, anticodons, rRNA
    plots/                              Individual composition figures
  submission/
    genome.fsa
    annotation.tbl
    metadata.json
    id_mapping.tsv
    validation.tsv
    checksums.sha256
comparative/
  taxonomy.tsv
  taxonomy_tree.nwk
  genome/{stats,plots}/
  repeats/{stats,plots}/
  protein_coding/{stats,plots}/
  noncoding/{stats,plots}/
pipeline_info/
  execution_report_<run>.html
  execution_trace_<run>.tsv
  execution_timeline_<run>.html
  pipeline_dag_<run>.html
```

## Interpretation

- Genome size, contig lengths, N50 and coverage lengths are in base pairs.
  GC percentage uses A/C/G/T positions as denominator. Softmask counts include
  all lowercase bases in the analysis genome. The uppercase assembly retains
  the same coordinates and reformatted contig names.
- Repeat coverage is the union of genomic intervals, regardless of strand.
  Counts are annotation rows. Coverage percentages for different repeat groups
  can overlap and must not be summed as a unique genome-wide coverage.
  In pooled-library mode, each species receives the shared consensus library
  actually used by RepeatMasker, not an independently inferred species library.
  Repeat landscape bins represent estimated divergence, not calibrated ages.
- Coding feature counts count GFF rows; the protein count counts selected
  transcripts. CDS FASTA preserves all annotated CDS bases; phase/codon_start
  determines translation. Protein lengths are amino-acid counts.
  QC and candidate tables preserve the before/after contamination-filtering
  comparison. The existing score/rank/BUSCO selection rules are retained;
  remaining exact ties are resolved by method name.
- Noncoding counts count retained RNA loci, not parent genes or exon rows.
  tRNAscan-SE takes priority over same-strand overlapping Infernal tRNA calls;
  Barrnap takes priority over Infernal rRNA calls. Other Infernal RNA calls
  remain. Unknown isotypes/types are explicit. Anticodons are included only
  when annotated, in 5-prime to 3-prime DNA notation (U normalized to T).
- Comparative matrices place genomes in columns. Tree order is used for
  overview plots and ncRNA heatmaps. The preprocessing tree is a taxonomy
  hierarchy, not a molecular phylogeny; displayed branch lengths are hierarchy
  depths. ncRNA correlations are descriptive Pearson/Spearman analyses, with
  no phylogenetic correction or multiple-testing correction. tRNA density is
  count / (assembly size / 1,000); its relationship with genome size is
  mathematically coupled through this denominator. Groups with too few genomes
  or constant values have NA coefficients and an explicit status.

## Optional and missing results

Outputs depend on the enabled modules and successful methods. After five
unsuccessful attempts BRAKER/GALBA still produce an empty GFF that is excluded
from method selection. When no coding annotation survives, the other available
categories are exported and `submission/validation.tsv` records the absence.
Empty categories retain empty/header-only sequence and statistics files; a
placeholder plot may state that no annotation was available. Anticodon plots
are omitted if no annotated anticodon exists. Detailed intermediate tool
outputs remain in the configured cache and work directories.

## Submission preparation

The combined GFF3, FASTA files and five-column feature table prepare an
annotation submission; they are not a validated submission. Check
`submission/validation.tsv`, provide the registered BioProject/BioSample and
locus-tag prefix, source/submitter metadata and submission template, then run
the target archive's validator (e.g. NCBI table2asn). Internal stop codons,
incomplete CDS and phase problems require review. FASTA protein identifiers are
local identifiers, not assigned UniProt accessions. ENA/Ensembl integration
requires the destination's own metadata and format requirements.

`submission/checksums.sha256` contains paths relative to the species directory.
This layout is produced by new runs; pre-existing legacy output directories are
not deleted or moved automatically.

## Reproducibility and configuration

`pipeline_info/provenance/` records resolved parameters, Nextflow/run metadata,
SHA-256 of input assemblies/local reads and relevant metadata files, pipeline
source fingerprints, database inventory/lock and available software version
records. NF-core versions.yml files are task output evidence. Additional cached
image version probes and the image inventory are labeled separately and do not
prove which cached image executed a past task. Runtime image paths are also
recorded by the Nextflow execution trace.

Rfam defaults to 15.0; the automatic OMA download uses an archived LUCA-v2.0.0
file with its published checksum. Program bank/image paths are defined in the main nextflow.config params block. Automatic UniProt evidence records the response release;
`uniprot_evidence_release` can reject a different release. To reproduce an old
REST snapshot, retain/reuse the exported protein evidence or supply an archived
protein_set; the REST service does not provide arbitrary past releases.

Input assemblies and protein evidence are fully hashed. Full database hashing
is optional (`provenance_hash_databases=true`), because banks can be very large.
The default database lock inventories large files by name, size and modification
time and hashes small release metadata. It is not a full content lock. Supplying
`database_lock_file` checks the recorded identities at provenance collection;
this does not invalidate or guard the analysis caches before computation.

`program_args` configures extra options for main tools; `program_resources`
sets CPU, memory and time overrides. Existing defaults apply when keys are
omitted. Required input/output arguments and wrapper-controlled formats should
not be overridden. Cache behavior is unchanged: changing options or bank
settings does not guarantee recomputation of already stored outputs.

`submission_metadata` is JSON indexed by genome ID, with source qualifiers,
BioProject/BioSample, registered locus-tag prefix and optional per-sequence
organelle/topology/genetic-code declarations. Original contig names are mapped
through the preprocessing identifier table. Organelles require an explicit
translation table. Codes are applied to exported CDS translations; nuclear
prediction tools receive the per-genome code when supported by their versions.
This does not implement a dedicated organellar gene predictor, origin-spanning
circular gene reconstruction or RNA-editing inference; those annotations need
specific review and are flagged before local submission validation.

Local official validation is opt-in (`validate_submission=true`) and requires
real metadata, an SBT template, a table2asn executable and its SHA-256. Reports
and generated SQN files are retained in `submission/ncbi_validation/`; errors
can stop the export. No remote submission is performed. BioProject/BioSample
association must be completed in the archive submission portal. ENA and Ensembl
still have their own import and validation procedures; the NCBI status does not
certify acceptance by those services.
