include { GENERATE_TAXONOMY_METADATA   } from "$projectDir/modules/local/preprocessing/generate_taxonomy_metadata.nf"
include { SEQKIT_STATS                 } from "$projectDir/modules/nf-core/seqkit/stats/main"
include { NORMALIZE_SEQKIT_STATS       } from "$projectDir/modules/local/preprocessing/normalize_seqkit_stats.nf"
include { MERGE_TABLES as MERGE_SEQKIT } from "$projectDir/modules/local/shared_modules/merge_tables.nf"
include { REFORMAT_GENOME              } from "$projectDir/modules/local/preprocessing/reformat_genome.nf"
include { SELECT_BUSCO_DATASET         } from "$projectDir/modules/local/preprocessing/select_busco_dataset.nf"
include { BUSCO_DOWNLOAD               } from "$projectDir/modules/nf-core/busco/download/main"
include { BUSCO_BUSCO                  } from "$projectDir/modules/nf-core/busco/busco/main"


workflow PREPROCESSING {

    take:
    genomes

    main:

    // Create a single taxonomy input table for all genomes.
    genomes
        .map { meta, genome ->
            [meta.id, meta.name, meta.taxid, meta.mnemonic ?: ''].join('\t')
        }
        .collectFile(
            name: 'taxonomy_input.tsv',
            seed: 'id\tname\ttaxid\tmnemonic\n',
            newLine: true,
            sort: true
        )
        .set { taxonomy_input }

    // Retrieve taxonomy metadata for all genomes.
    GENERATE_TAXONOMY_METADATA(taxonomy_input)

    GENERATE_TAXONOMY_METADATA.out.json
        .flatMap { taxonomyFile ->
            def taxonomyById = new groovy.json.JsonSlurperClassic().parseText(taxonomyFile.text)

            taxonomyById.collect { id, taxonomy ->
                [id, [lineage: taxonomy.lineage, mnemonic: taxonomy.mnemonic]]
            }
        }
        .set { taxonomy_metadata }

    // Compute assembly statistics.
    SEQKIT_STATS(genomes)
    NORMALIZE_SEQKIT_STATS(SEQKIT_STATS.out.stats)

    NORMALIZE_SEQKIT_STATS.out.stats
        .map { meta, seqkitFile ->
            [meta.id, PreprocessingUtils.parseAssemblySize(seqkitFile)]
        }
        .set { assembly_sizes }

    NORMALIZE_SEQKIT_STATS.out.stats
        .map { meta, seqkitFile -> ['size', seqkitFile] }
        .groupTuple(by: 0)
        .map { id, files -> [[id: id], files, 'all_seqkit.stats'] }
        .set { seqkit_tables }

    MERGE_SEQKIT(seqkit_tables)

    // Add taxonomy metadata and assembly size.
    genomes
        .map { meta, genome -> [meta.id, meta, genome] }
        .join(taxonomy_metadata, failOnDuplicate: true, failOnMismatch: true)
        .join(assembly_sizes, failOnDuplicate: true, failOnMismatch: true)
        .map { id, meta, genome, taxonomyMeta, assemblySize ->
            [meta + taxonomyMeta + [assembly_size: assemblySize], genome]
        }
        .set { genomes_with_metadata }

    REFORMAT_GENOME(genomes_with_metadata)

    // Select or use one BUSCO dataset for all genomes.
    if (params.busco_lineage?.toString()?.trim()) {
        log.info "Using the provided BUSCO lineage: ${params.busco_lineage}"

        busco_dataset = Channel.value(params.busco_lineage)
    } else {
        log.info 'Selecting a BUSCO dataset from the shared taxonomy'

        SELECT_BUSCO_DATASET(GENERATE_TAXONOMY_METADATA.out.table)

        SELECT_BUSCO_DATASET.out
            .map { file -> file.text.trim() }
            .set { busco_dataset }
    }

    busco_dataset.view { dataset -> "BUSCO dataset: ${dataset}" }

    // Add the shared BUSCO dataset to genome metadata.
    REFORMAT_GENOME.out.fasta
        .combine(busco_dataset)
        .map { meta, genome, buscoDataset ->
            [meta + [busco_dataset: buscoDataset], genome]
        }
        .set { preprocessed_genomes }

    // Download once and broadcast the root directory to all BUSCO tasks.
    BUSCO_DOWNLOAD(busco_dataset)
    busco_download_dir = BUSCO_DOWNLOAD.out.download_dir.first()

    // Preserve the dataset directory needed by later annotation stages.
    busco_dataset_dir = busco_download_dir.combine(busco_dataset)
        .map { root, dataset -> [[id: dataset], root.resolve("lineages/${dataset}")] }

    BUSCO_BUSCO(preprocessed_genomes, 'genome', busco_dataset.first(), busco_download_dir, [], false)

    emit:
    versions = SEQKIT_STATS.out.versions.mix(BUSCO_DOWNLOAD.out.versions).mix(BUSCO_BUSCO.out.versions)
    assembly_stats = NORMALIZE_SEQKIT_STATS.out.stats
    genome_mapping = REFORMAT_GENOME.out.mapping
    busco_summaries = BUSCO_BUSCO.out.short_summaries_json.mix(BUSCO_BUSCO.out.short_summaries_txt)
    preprocessed_genomes = preprocessed_genomes
    taxonomy_table = GENERATE_TAXONOMY_METADATA.out.table
    taxonomy_tree = GENERATE_TAXONOMY_METADATA.out.minimal_tree
    genome_stats = MERGE_SEQKIT.out
    busco_dataset = busco_dataset
    busco_dataset_dir = busco_dataset_dir
    busco_assembly_results = BUSCO_BUSCO.out.short_summaries_json
}