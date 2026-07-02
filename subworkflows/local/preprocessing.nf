include { GET_TAXONOMY } from "$projectDir/modules/local/module_preprocessing/get_taxonomy.nf"
include { SEQKIT_STATS } from "$projectDir/modules/local/module_preprocessing/seqkit_stats.nf"
include { REFORMAT_GENOME } from "$projectDir/modules/local/module_preprocessing/reformat_genome.nf"
include { DOWNLOAD_NEWICK } from "$projectDir/modules/local/module_preprocessing/download_newick.nf"
include { SELECT_BUSCO_DATASET } from "$projectDir/modules/local/module_preprocessing/select_busco_dataset.nf"

def extract_lineage(lineage_file) {
    def lines = lineage_file.text.readLines()

    def dict = [:]
    lines.each { line ->
        line = line.trim()
        if (!line || line.startsWith('#')) return

        def cols = line.split('\t', -1)
        if (cols[0] == 'rank') return   // header

        def rank  = cols[0]
        def taxid = cols[1]
        def name  = cols[2]

        dict[rank] = [ name: name, taxid: taxid ]
    }
    return dict
}


def transformList(input) {
    def ids = input[0]
    def mnemonic = input[1]

    def range = 1..ids.size()
    return range.collect { i -> [ids[i-1], mnemonic + i]}
}

workflow PREPROCESSING {
    take:
        genomes

    main:
        // fetch lineage from NCBI taxonomy and mnemonic from Uniprot
        GET_TAXONOMY(genomes)

        // compute genome stats with seqkit
        SEQKIT_STATS(genomes, 'dna')

        // add lineage, mnemonic and genome size to the meta map
        GET_TAXONOMY.out.lineage
            .map { id, meta, lineage -> [ id, extract_lineage(lineage) ] }
            .set { lineage }

        GET_TAXONOMY.out.mnemonic
            .map { id, meta, mnemonic -> [ id, mnemonic.splitText()[0].replaceAll('\n', '') ] }
            .set { official_mnemonic }

        official_mnemonic
            .map { id, mnemonic -> [ id, mnemonic.startsWith("9") ? id.split("-")[0].substring(0, 3).toUpperCase() + id.split("-")[1].substring(0, 2).toUpperCase() : mnemonic ] }
            // check if there is two identical mnemotic, add for the first a 1 at end of mnemonic, add 2 for the second, 3 for the last...
            .groupTuple(by: 1)
            .map { it[0].size() > 1 ?  transformList(it) : [it[0][0], it[1]] }
            .flatMap { item ->
                // If the element is already flat, return it as is
                if (!(item[0] instanceof List)) {
                    return [item]
                } else {
                    // If the element is a nested list, flatten it
                    return item
                }
            }
            .set { exogap_mnemonic }

        SEQKIT_STATS.out
            .map { id, meta, seqkit -> [ id, file(seqkit).splitText()[1].replaceAll('\n', '').split('\t')[4] ]}
            .set { size }
        genomes
            .join(lineage)
            .join(official_mnemonic)
            .join(exogap_mnemonic)
            .join(size)
            .map { id, meta, genome, lineage, official_mnemonic, exogap_mnemonic, size ->
                def updatedMeta = meta + ["lineage": lineage, "mnemonic": ['official': official_mnemonic, 'exogap': exogap_mnemonic], "assembly_size": size]
                [ id, updatedMeta, genome ]
            }
            .set { genomes }

        REFORMAT_GENOME(genomes)

        if (params.comparative_genomics) {
            // download newick
            genomes
                .map {id, meta, genome -> [id, meta.taxid].join('/')}
                .collect()
                .map { it.join(' ') }
                .set { taxids }

            DOWNLOAD_NEWICK(taxids)
        }

        println params.busco_lineage
        // If `busco_lineage` is not provided (missing), null, or an empty string
        if (!params.busco_lineage?.toString()?.trim()) {
            print "Selecting BUSCO dataset based on taxid"
            SELECT_BUSCO_DATASET(GET_TAXONOMY.out.lineage)
            SELECT_BUSCO_DATASET.out
                .map { id, meta, dataset -> [id, dataset.text.readLines()[0]] }
                .join (REFORMAT_GENOME.out)
                .map { id, busco_dataset, meta, genome -> [id, meta + ["busco_dataset": busco_dataset], genome] }
                .set { genomes }
        } else {
            print "Using provided BUSCO lineage: ${params.busco_lineage}"
            genomes
                .map { id, meta, genome -> [id, meta + ["busco_dataset": params.busco_lineage], genome] }
                .set { genomes }
        }
        
    emit:
        genomes = genomes
        newick = DOWNLOAD_NEWICK.out
        genome_stats = SEQKIT_STATS.out
}