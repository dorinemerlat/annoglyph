![nf-core/exogap](./docs/images/exogap_logo_light.png#gh-light-mode-only)
![nf-core/exogap](./docs/images/exogap_logo_dark.png#gh-dark-mode-only)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

EXOGAP (EXotic Organism Genome Annotation Pipeline) is a Nextflow DSL2 workflow for comprehensive annotation of genomes from non-model species. It integrates tools for repetitive element identification, protein-coding gene prediction and non-coding gene prediction. The pipeline can annotate multiple related genomes in a single run and formats outputs for deposition in public databases.

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
git clone https://github.com/dorinemerlat/exogap.git
cd exogap
```

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
- fasta (required): path to fasta assembly file (absolute or relative).
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
| RNA-seq | STAR: at most 3 eligible SRA runs per species, each at most 1 billion bases | STAR: all eligible runs discovered for each species, without these two limits |
| Results/cache | Separate `minimal/` directories | Separate `full/` directories |

The RNA-seq size limit **excludes larger runs**; it does not truncate or
subsample reads. No run satisfying the filters means no SRA RNA-seq evidence
for that species. The provided samplesheets have no local RNA-seq files;
limits do not subsample local FASTQs if you add them. Both modes use STAR so
that the minimal RNA-seq limits apply throughout the test (VARUS uses a
separate adaptive sampling strategy).

Prepare evidence and regenerate statistics from already installed assemblies:

```bash
python3 test/prepare_test_sets.py
```

For a fresh download, see the [installation instructions](test/README.md).
The download script recreates dataset directories; use the preparation command
above to retain existing results and caches.

From the project root, load your Nextflow/Singularity environment. On IFB:

```bash
source /etc/profile.d/modules.sh
module load nextflow/25.10.2 singularity
export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-$PWD/.singularity}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-$PWD/.apptainer_tmp}"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
export TMPDIR="$APPTAINER_TMPDIR"
mkdir -p "$APPTAINER_TMPDIR"
```

Run one command for the desired group and mode:

```bash
# Hanseniaspora
nextflow run main.nf -profile singularity,test_hanseniaspora_minimal -c conf/ifb.config -c test/validation/slurm_runtime.config
nextflow run main.nf -profile singularity,test_hanseniaspora_full -c conf/ifb.config -c test/validation/slurm_runtime.config

# Pieris
nextflow run main.nf -profile singularity,test_pieris_minimal -c conf/ifb.config -c test/validation/slurm_runtime.config
nextflow run main.nf -profile singularity,test_pieris_full -c conf/ifb.config -c test/validation/slurm_runtime.config

# Utricularia
nextflow run main.nf -profile singularity,test_utricularia_minimal -c conf/ifb.config -c test/validation/slurm_runtime.config
nextflow run main.nf -profile singularity,test_utricularia_full -c conf/ifb.config -c test/validation/slurm_runtime.config
```

Alternatively, submit the launcher with any of these profiles:

```bash
sbatch nextflow.sbatch test_hanseniaspora_minimal
sbatch nextflow.sbatch test_pieris_full
sbatch nextflow.sbatch test_utricularia_minimal -resume
```

Add `-resume` to resume the same mode. Outputs are written to
`test/test_<group>/<minimal|full>/results`, with persistent caches alongside
in `cache`. Existing outputs at the older locations are retained, but these
new profiles use separate directories. The legacy `test_hanseniaspora`,
`test_pieris` and `test_utricularia` names now select the minimal mode.
Use your site's execution config instead of the two IFB configs elsewhere.
Bank paths and program options are configured in `nextflow.config`.

Full mode uses all **available evidence within these fixtures**, not every
protein database or every possible public RNA-seq dataset. SRA availability
can change; the workflow records run metadata. Protein evidence comes from
the tested assemblies themselves: these are integration tests, not independent
accuracy benchmarks. U. reniformis uses protein homologs from the other
Utricularia assemblies because its downloaded package has no protein annotation.
The full plant/insect runs can be substantial; minimal mode reduces evidence,
not genome size. These commands describe the configured tests, not a guarantee
that a complete run has already passed.

## Pipeline output
Results are organized as `species/<genome_id>/` with `genome`, `repeats`,
`protein_coding`, `noncoding` and `submission` directories. Multi-species
statistics and figures are in `comparative/`; execution and provenance reports
are in `pipeline_info/`. A single results README describes these files.
See [annotation/export documentation](docs/ncrna_and_export.md) for submission
metadata, genetic codes and optional local NCBI validation.

## Credits

Originally written by Dorine Merlat (dorine.merlat@etu.unistra.fr).
Thanks to Arnaud Kress and Odile Lecompte for assistance.

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
