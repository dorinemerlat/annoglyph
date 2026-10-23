# AnnoGlyph

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

AnnoGlyph is a Nextflow DSL2 workflow for comprehensive annotation of genomes from non-model species. It integrates tools for repetitive element identification, protein-coding gene prediction and non-coding gene prediction. The pipeline can annotate multiple related genomes in a single run and formats outputs for deposition in public databases.

Main modules:
- Repetitive element annotation: RepeatModeler, RepeatMasker (uses FamDB).
- Protein-coding gene annotation: BRAKER2, BRAKER3, GALBA, STAR, BLASTP, InterProScan, eggNOG-mapper; BUSCO and OMArk quality control with automatic prediction selection.
- Non-coding gene annotation: tRNAscan-SE, Barrnap and Infernal; overlapping tRNA/rRNA predictions are reconciled before export.

## Installation

Prerequisites:
- Nextflow (recommended >= 24.04)
- Singularity / Apptainer

Clone the repository:
```bash
git clone https://github.com/dorinemerlat/annoglyph.git
cd annoglyph
```

### Python reporting container

The statistics and export scripts require SciPy in addition to the other
Python dependencies. Build the local Singularity/Apptainer image once:

```bash
bash containers/annoglyph_python/build_image.sh
```

The image path is configured by `params.annoglyph_python_container` in
`nextflow.config`. The build uses the published base image and adds the pinned
SciPy dependency; run it with the resources required by your computing site.

## Usage

### Samplesheet CSV file

Prepare a samplesheet CSV describing genomes to annotate.
## Pipeline inputs & options

- **Samplesheet (required):** a CSV describing genomes. Required columns: `name`, `taxid`, `fasta`. Optional columns supported: `RNASeq-dir` (path to directory with FASTQ files), `SRA` (SRA accessions), `repeats_set` (FASTA library).
- **Main config:** `nextflow.config` (or `example/nextflow_example.config`) controls runtime options. Key params:
	- `input`: path to samplesheet CSV (or pass `--input` on CLI).
	- `module_repeats`, `module_genes`, `module_ncgenes`: enable/disable modules (booleans).
	- `group_consensus_sequences`: (repeats) boolean to share consensus across genomes.
	- `reference_library`: optional repeats library FASTA for repeat annotation.
	- `protein_set`: (genes) path to protein FASTA used for gene annotation.
	- `omark_taxonomy_db`: ETE3 NCBI taxonomy SQLite database for OMArk; an existing `~/.etetoolkit/taxa.sqlite` is detected automatically.
	- `outdir`: output directory (default `out`).
	- `publish_dir_mode`: how results are published (default `copy`).
	- `max_memory`, `max_cpus`, `max_time`: global resource caps for processes.
	- `env.NXF_SINGULARITY_CACHEDIR`: where Singularity/Apptainer caches images (default `$projectDir/.singularity/`).

### Example runs

Run with your config and samplesheet:
```bash
nextflow run . -c nextflow.config --input samplesheet.csv -profile singularity
```

Use the example config:
```bash
nextflow run . -c example/nextflow_example.config --input example/samplesheet_example.csv -profile singularity
```

To continue after fixing an error, add `-resume`:
```bash
nextflow run . -c nextflow.config --input samplesheet.csv -profile singularity -resume
```

Minimal example samplesheet (samplesheet.csv):
```
name,taxid,fasta
Felis catus,9685,/path/to/felis-catus.fa
```

Column description:

- name (required): display name used in reports (alphanumeric, spaces, '-' and '_' allowed).
- taxid (required): NCBI taxon id (use nearest parent taxid if exact id not available).
- fasta (required): path to fasta assembly file.
- RNASeq-dir (optional): path to directory containing RNA-Seq FASTQ files for this genome (single-end or paired-end, .fastq or .fastq.gz).
- SRA (optional): comma-separated list of SRA accessions for RNA-Seq data

### Config file

Fill the `nextflow.config` file.

### Launch pipeline

Run the pipeline (example):
```
nextflow run main.nf -profile <singularity|docker|...> --config nextflow.config
```

## Integration tests

Three biological groups are provided: **Hanseniaspora** (5 yeast species),
**Pieris** (3 butterfly species) and **Utricularia** (3 plant species).
See the [assembly statistics table](test/assembly_stats.md) for species,
assembly accessions, sizes, sequence counts, N50 and GC content
([TSV](test/assembly_stats.tsv)). These are three groups, covering 11 assemblies.

| Group | Assemblies | Minimal proteins | Full proteins |
| --- | ---: | ---: | ---: |
| Hanseniaspora | 5 | 2,000 | 21,786 |
| Pieris | 3 | 2,000 | 56,761 |
| Utricularia | 3 | 2,000 | 74,162 |

Counts refer to the currently installed packages and can change after downloading
new annotations.

Each group has two modes. Both annotate the **complete assemblies** and enable
repeats, protein-coding genes, non-coding genes and comparative reports.

| Evidence | `minimal` | `full` |
| --- | --- | --- |
| Proteins | 2,000 distinct proteins, 90–5,000 amino acids | All distinct usable proteins from the downloaded annotations, without a length limit |
| RNA-seq | STAR: at most 3 eligible SRA runs per species, each at most 1 billion bases, capped at 500,000 spots per run | STAR: all eligible runs discovered for each species, without these limits |
| Results/cache | Separate `minimal/` directories | Separate `full/` directories |

The RNA-seq size filter **excludes larger runs**. Each selected run is then
limited to its first 500,000 spots; a paired read pair counts as one spot.
Full profiles retain all reads. No run satisfying the filters means no SRA RNA-seq evidence
for that species. The provided samplesheets have no local RNA-seq files;
limits do not subsample local FASTQs if you add them. Both modes use STAR so
that the minimal RNA-seq limits apply throughout the test (VARUS uses a
separate adaptive sampling strategy).

Minimal profiles require more than 100,000 mapped alignments per RNA-seq BAM
to enable RNA-supported prediction (`rnaseq_min_mapped_reads`). The default
threshold for full profiles and production runs remains 1,000,000.

Prepare evidence and regenerate statistics from already installed assemblies:

```bash
python3 test/prepare_test_sets.py
```

For a fresh download, see the [installation instructions](test/README.md).
The download script recreates dataset directories; use the preparation command
above to retain existing results and caches.

From the project root, load Nextflow and Singularity/Apptainer using the
instructions for your computing environment.

Run one command for the desired group and mode:

```bash
# Hanseniaspora
nextflow run main.nf -profile singularity,test_hanseniaspora_minimal
nextflow run main.nf -profile singularity,test_hanseniaspora_full

# Pieris
nextflow run main.nf -profile singularity,test_pieris_minimal
nextflow run main.nf -profile singularity,test_pieris_full

# Utricularia
nextflow run main.nf -profile singularity,test_utricularia_minimal
nextflow run main.nf -profile singularity,test_utricularia_full
```

Add your site configuration with `-c conf/local.config` when needed.
Machine-specific configurations and batch launchers are kept locally and are not
versioned.

Add `-resume` to resume the same mode. Outputs are written to
`test/test_<group>/<minimal|full>/results`, with persistent caches alongside
in `cache`. Existing outputs at the older locations are retained, but these
new profiles use separate directories. The legacy `test_hanseniaspora`,
`test_pieris` and `test_utricularia` names now select the minimal mode.
General program options are configured in `nextflow.config`. Local bank paths
and execution settings belong in the ignored `conf/local.config`, loaded
automatically when present.

Full profiles enable genus-exclusion evaluation: protein evidence is downloaded
from OrthoDB and UniProtKB using a shared clade above genus, then filtered for
each genome to remove its entire genus. All eligible RNA-seq reads are used.
Minimal profiles retain the prepared 2,000-protein fixtures and leave evaluation
disabled; these are integration tests, not independent accuracy benchmarks.
Full plant/insect runs and automatic evidence downloads can be substantial.
Neither mode reduces genome size. These commands describe the configured tests,
not a guarantee that a complete run has already passed.

### Genus-exclusion evaluation

Set `evaluation_mode = true` to exclude the annotated genome's genus from
**automatically downloaded** protein evidence before BRAKER2, BRAKER3 and GALBA.
The default is `false`. With multiple genera, each genome receives evidence
filtered against its own genus; other target genera are not removed from its set.
The downloaded source FASTA remains unchanged.

```groovy
params {
    evaluation_mode = true
    protein_set = null
    protein_evidence_taxonomy_db = '/path/to/taxa.sqlite'
}
```

The taxonomy database must be an ETE NCBI taxonomy SQLite database. When
`protein_evidence_taxonomy_db` is null, `omark_taxonomy_db` is reused. No database
is downloaded or updated implicitly by this filter. UniProtKB `OX=` taxids and
partitioned OrthoDB taxid-prefixed identifiers are resolved through this database.
Proteins with unresolved taxonomy or no genus are excluded and counted, rather
than retained with uncertain provenance. Missing target genus information or an
empty filtered set stops the run with an explicit error.

A JSON report per genome records the excluded genus and counts of retained,
same-genus and unresolved proteins. Reports are published in
`pipeline_info/provenance/database_metadata/protein_evidence_filter_<genome_id>.json`.
Evaluation predictions and exports use a separate persistent cache namespace
from ordinary annotations. Keep the same mode when resuming.

**User-supplied `protein_set` files are always accepted unchanged**, even when
`evaluation_mode` is enabled: AnnoGlyph makes no assumptions about their header
format or taxonomy. Filter such files beforehand if needed. This option filters
structural protein evidence only; it does not filter functional annotation banks,
BUSCO lineages, RNA-seq or other reference resources, and does not by itself make
an evaluation independent of every source of prior annotation.

## Execution statistics across resumed runs

Give an experiment an identifier and keep it when resuming:

```bash
nextflow run main.nf -profile singularity --input samplesheet.csv --run_group experiment_01
nextflow run main.nf -profile singularity --input samplesheet.csv --run_group experiment_01 -resume
```

Add your site's execution configuration as usual. `run_group` prefixes the trace,
report, timeline and DAG filenames without replacing Nextflow's individual run
name. It accepts letters, digits, underscores, dots and hyphens, starting with a
letter or digit. Use a new group for a different experiment or changed inputs.
Reports retain a timestamp, so resumes do not overwrite previous traces.

Generate unified tables and figures after a run or failure:

```bash
python bin/summarize_execution.py \
    --trace-dir results/pipeline_info --run-group experiment_01 \
    --output results/execution_stats --plots
```

On Slurm, add `--slurm-accounting` to retrieve consumed CPU from `sacct` using
trace job IDs. Otherwise CPU hours are estimated as task runtime × CPU usage / 100.
If an accounting record is unavailable, the script falls back to that estimate
and records its source. Plotting requires Matplotlib; tables use only Python's
standard library.

Outputs are `execution_summary.tsv` (global statistics), `process_summary.tsv`
(one row per fully qualified module), `tasks.tsv` (individual executions),
`summary_metadata.json` (input traces and limitations), `README.md` (readable
summary), and `plots/` (PNG/SVG).
CPU and memory measurement counts show coverage. Cached trace records are counted
separately; repeated execution records are deduplicated, while real retries and
recomputations remain included, with failed/aborted CPU reported separately.
Modules skipped through `storeDir` can be absent from a trace: include their
original traces to recover those measurements.

Plots use **Seaborn box plots with individual task points**, arranged in four
vertically stacked panels in `plots/execution_metrics.png` (also SVG). All panels
share the same modules and workflow groups on the X axis; each has its own
logarithmic Y axis and bold vertical metric label. Requested limits appear as
purple dashed segments without markers. Crosses distinguish failed/aborted tasks from completed tasks.
Choose one selection mode:

- `--plot-selection tools` (default): annotation tools such as RepeatMasker,
  RepeatModeler, MCHelper, GALBA, BRAKER2/3, BLAST, InterProScan, eggNOG-mapper
  and noncoding RNA predictors, without resource thresholds. Utilities such as
  AGAT are excluded. Repeat `--plot-tool 'BRAKER*' --plot-tool 'BLAST*'` to
  replace the default tool list with module names or glob patterns.
- `--plot-selection thresholds`: consider every module, retaining those with at
  least one task above any metric's threshold; retained modules appear on all
  four panels, with all their measured observations.
- `--plot-selection all`: show every module with measurements, without thresholds.

Thresholds apply only in `thresholds` mode:

| Option | Metric | Default threshold |
| --- | --- | ---: |
| `--runtime-threshold` | Task runtime, minutes | 30 |
| `--cpu-threshold` | Consumed CPU, CPU hours | 10 |
| `--memory-threshold` | Peak resident memory, GiB | 8 |
| `--cores-threshold` | Mean CPU use, equivalent cores | 4 |

A threshold of `0` includes all measured modules for that metric. Missing values
are never plotted as zeros. Mean cores and CPU efficiency use Slurm elapsed time with Slurm CPU measurements,
or Nextflow runtime with Nextflow estimates. Requested CPU/memory resources and
consumed CPU are kept distinct in the tables. All modules remain in the tables in every plot mode.

Plot labels use readable tool names (e.g. RepeatMasker, InterProScan), with a
`module_labels.tsv` mapping to the unchanged full process identifiers. Y axes
use a logarithmic scale, labelled with the metric unit. When zeros are present,
a linear region near zero preserves those observations (`symlog`).

Add `--show-requested` to overlay each module's **maximum requested limit** on
the runtime, memory and mean-cores plots, using trace `time`, `memory` and `cpus`.
These maxima also appear in `process_summary.tsv`, alongside observed use.
They include the selected executions and retries, so they can reflect resource
increases between attempts. No requested limit is drawn on consumed CPU hours:
CPU allocation × time limit is only a theoretical budget. With large reservations,
the shared axis deliberately shows the gap between requested and observed use.

**Elapsed time needs care:** `observed_span_hours` includes gaps between resumes;
`observed_active_hours` is the union of observed task start/end intervals;
`cumulative_task_hours` sums runtime across tasks, including parallel work.
None is an exact reconstruction of a fresh, uninterrupted pipeline run. The
largest task memory peak is not the simultaneous memory of the whole pipeline.
Controller resources are not included. See the generated `README.txt` for details.

For existing traces without `run_group`, pass explicit paths or omit the group
filter in a directory containing only the experiment's traces:

```bash
python bin/summarize_execution.py execution_trace_*.tsv --output execution_stats
```

The active tool modules, preserved options and compatibility patches are listed
in [nf-core integration](docs/nfcore_modules.md).

## Pipeline output
Results are organized as `species/<genome_id>/` with `genome`, `repeats`,
`protein_coding`, `noncoding` and `submission` directories. Multi-species
statistics and figures are in `comparative/`; execution and provenance reports
are in `pipeline_info/`. A single results README describes these files.
See [annotation/export documentation](docs/ncrna_and_export.md) for submission
metadata, genetic codes and optional local NCBI validation.

## Credits

Originally written by Dorine Merlat (dorine.merlat@etu.unistra.fr).
Thanks to Odile Lecompte, Arnaud Kress and Yannis Nevers for their assistance.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For issues or support, please open an issue on the pipeline GitHub repository.

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/test for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->


An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
