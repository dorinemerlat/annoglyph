include { GENERATE_TAXONOMY_METADATA   } from "$projectDir/modules/preprocessing/generate_taxonomy_metadata.nf"
include { SEQKIT_STATS                 } from "$projectDir/modules/preprocessing/seqkit_stats.nf"
include { MERGE_TABLES as MERGE_SEQKIT } from "$projectDir/modules/shared_modules/merge_tables.nf"
include { REFORMAT_GENOME              } from "$projectDir/modules/preprocessing/reformat_genome.nf"
include { SELECT_BUSCO_DATASET         } from "$projectDir/modules/preprocessing/select_busco_dataset.nf"
include { DOWNLOAD_BUSCO_DATASET       } from "$projectDir/modules/preprocessing/download_busco_dataset.nf"
include { BUSCO                        } from "$projectDir/modules/shared_modules/busco.nf"


workflow PREPROCESSING {

    take:
    genomes

    main:

    // Create a single taxonomy input table for all genomes.
    genomes
        .map { id, meta, genome ->
            [id, meta.name, meta.taxid, meta.mnemonic ?: ''].join('\t')
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
            def taxonomyById = new groovy.json.JsonSlurper().parseText(taxonomyFile.text)

            taxonomyById.collect { id, taxonomy ->
                [id, [lineage: taxonomy.lineage, mnemonic: taxonomy.mnemonic]]
            }
        }
        .set { taxonomy_metadata }

    // Compute assembly statistics.
    SEQKIT_STATS(genomes, 'dna')

    SEQKIT_STATS.out
        .map { id, meta, seqkitFile ->
            [id, PreprocessingUtils.parseAssemblySize(seqkitFile)]
        }
        .set { assembly_sizes }

    SEQKIT_STATS.out
        .map { id, meta, seqkitFile -> ['size', seqkitFile] }
        .groupTuple(by: 0)
        .map { id, files -> [id, files, 'all_seqkit.stats'] }
        .set { seqkit_tables }

    MERGE_SEQKIT(seqkit_tables)

    // Add taxonomy metadata and assembly size.
    genomes
        .join(taxonomy_metadata)
        .join(assembly_sizes)
        .map { id, meta, genome, taxonomyMeta, assemblySize ->
            [id, meta + taxonomyMeta + [assembly_size: assemblySize], genome]
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
        .map { id, meta, genome, buscoDataset ->
            [id, meta + [busco_dataset: buscoDataset], genome]
        }
        .set { preprocessed_genomes }

    // Download the shared BUSCO dataset once.
    DOWNLOAD_BUSCO_DATASET(busco_dataset)

    // Associate every genome with the downloaded BUSCO dataset.
    preprocessed_genomes
        .map { id, meta, genome ->
            [meta.busco_dataset, id, meta, genome]
        }
        .combine(DOWNLOAD_BUSCO_DATASET.out, by: 0)
        .map { buscoDataset, id, meta, genome, buscoDatasetPath ->
            [id, meta, genome, 'genome', "${id}_on_assembly", buscoDataset, buscoDatasetPath]
        }
        .set { genomes_for_busco }

    BUSCO(genomes_for_busco)
 
    emit:
    preprocessed_genomes = preprocessed_genomes
    taxonomy_table = GENERATE_TAXONOMY_METADATA.out.table
    taxonomy_tree = GENERATE_TAXONOMY_METADATA.out.minimal_tree
    genome_stats = MERGE_SEQKIT.out
    busco_dataset = busco_dataset
    busco_dataset_path = DOWNLOAD_BUSCO_DATASET.out
    busco_assembly_results = BUSCO.out.json
}