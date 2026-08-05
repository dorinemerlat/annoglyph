include { FETCH_TAXONOMY_METADATA               } from "$projectDir/modules/preprocessing/fetch_taxonomy_metadata.nf"
include { SEQKIT_STATS                          } from "$projectDir/modules/preprocessing/seqkit_stats.nf"
include { MERGE_TABLES as MERGE_SEQKIT          } from "$projectDir/modules/shared_modules/merge_tables.nf"
include { REFORMAT_GENOME                       } from "$projectDir/modules/preprocessing/reformat_genome.nf"
include { DOWNLOAD_NEWICK                       } from "$projectDir/modules/preprocessing/download_newick.nf"
include { SELECT_BUSCO_DATASET                  } from "$projectDir/modules/preprocessing/select_busco_dataset.nf"
include { DOWNLOAD_BUSCO_DATASET                } from "$projectDir/modules/preprocessing/download_busco_dataset.nf"
include { BUSCO                                 } from "$projectDir/modules/shared_modules/busco.nf"

workflow PREPROCESSING {

    take:
    genomes

    main:

    // Retrieve NCBI taxonomy and UniProt mnemonics.
    FETCH_TAXONOMY_METADATA(genomes)

    // Compute assembly statistics.
    SEQKIT_STATS(genomes, 'dna')

    FETCH_TAXONOMY_METADATA.out.lineage
        .map { id, meta, lineageFile ->
            [ id, PreprocessingUtils.parseLineage(lineageFile) ]
        }
        .set { lineage }

    FETCH_TAXONOMY_METADATA.out.mnemonic
        .map { id, meta, mnemonicFile ->
            [ id, PreprocessingUtils.readFirstLine(mnemonicFile) ]
        }
        .set { official_mnemonics }

    official_mnemonics
        .map { id, mnemonic ->
            def exogapMnemonic = mnemonic.startsWith('9')
                ? PreprocessingUtils.generateFallbackMnemonic(id)
                : mnemonic

            [id, exogapMnemonic]
        }
        .groupTuple(by: 1)
        .flatMap { genomeIds, mnemonic ->
            PreprocessingUtils.makeUniqueMnemonics(genomeIds, mnemonic)
        }
        .set { exogap_mnemonics }

    SEQKIT_STATS.out
        .map { id, meta, seqkitFile ->
            [ id, PreprocessingUtils.parseAssemblySize(seqkitFile) ]
        }
        .set { assembly_sizes }

    SEQKIT_STATS.out
        .map { id, meta, seqkitFile ->
            [ 'size', seqkitFile ]
        }
        .groupTuple(by: 0)
        .map { id, files ->
            [ id, files, 'all_seqkit.stats' ]
        }
        .set { seqkit_tables }

    MERGE_SEQKIT(seqkit_tables)

    genomes
        .join(lineage)
        .join(official_mnemonics)
        .join(exogap_mnemonics)
        .join(assembly_sizes)
        .map { id, meta, genome, lineageData, officialMnemonic, exogapMnemonic, assemblySize ->

            def mnemonics = PreprocessingUtils.selectMnemonics(
                meta,
                officialMnemonic,
                exogapMnemonic
            )

            def updatedMeta = meta + [
                lineage: lineageData,
                mnemonic: mnemonics,
                assembly_size: assemblySize
            ]

            [ id, updatedMeta, genome ]
        }
        .set { genomes_with_metadata }

    REFORMAT_GENOME(genomes_with_metadata)

    newick = Channel.empty()
    
    if (params.comparative_genomics) {
        genomes
            .map { id, meta, genome ->
                "\"${meta.name}\"/${meta.taxid}"
            }
            .collect()
            .map { taxonEntries ->
                taxonEntries.join(' ')
            }
            .set { taxa }

        DOWNLOAD_NEWICK(taxa)

        newick = DOWNLOAD_NEWICK.out
    }

    if (!params.busco_lineage?.toString()?.trim()) {
        log.info 'Selecting BUSCO datasets from the input taxonomy'

        SELECT_BUSCO_DATASET(FETCH_TAXONOMY_METADATA.out.lineage)
        SELECT_BUSCO_DATASET.out
            .map { id, meta, datasetFile ->
                [
                    id,
                    PreprocessingUtils.readFirstLine(datasetFile)
                ]
            }
            .set { busco_datasets }
        
        busco_datasets
            .join(REFORMAT_GENOME.out.fasta)
            .map { id, busco_dataset, meta, genome ->
                [
                    id,
                    meta + [busco_dataset: busco_dataset],
                    genome
                ]
            }
            .set { preprocessed_genomes }
        
        // Log the selected BUSCO datasets for each genome.
        busco_datasets
            .map { id, busco_dataset -> busco_dataset }
            .unique()
            .collect()
            .view { datasets ->
                "Selected BUSCO datasets: ${datasets.join(', ')}"
        }

    } else {
        log.info "Using the provided BUSCO lineage: ${params.busco_lineage}"

        REFORMAT_GENOME.out.fasta
            .map { id, meta, genome ->
                [
                    id,
                    meta + [busco_dataset: params.busco_lineage],
                    genome
                ]
            }
            .set { preprocessed_genomes }

        busco_datasets = preprocessed_genomes
            .map { id, meta, genome ->
                [id, meta.busco_dataset]
            }
    }

    // Download each required BUSCO dataset once.
    busco_datasets
        .map { id, busco_dataset -> busco_dataset }
        .unique()
        .set { busco_datasets_to_download }

    DOWNLOAD_BUSCO_DATASET(busco_datasets_to_download)

    // Associate each genome with its downloaded BUSCO dataset.
    preprocessed_genomes
        .map { id, meta, genome ->
            [ meta.busco_dataset, id, meta, genome ]
        }
        .combine(DOWNLOAD_BUSCO_DATASET.out.dataset, by: 0)
        .map { busco_dataset, id, meta, genome, busco_dataset_path ->
            [ id, meta, genome, 'genome', "${id}_on_assembly", busco_dataset, busco_dataset_path ]
        }
        .set { genomes_for_busco }

    BUSCO( genomes_for_busco )

    emit:
    genomes = preprocessed_genomes
    newick = newick
    genome_stats = SEQKIT_STATS.out
    busco_datasets = DOWNLOAD_BUSCO_DATASET.out
    busco_on_assembly = BUSCO.out.json
}