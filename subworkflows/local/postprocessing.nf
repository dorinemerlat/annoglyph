/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CHECK IF FASTA FILES ARE VALIDS, REFORMATE THEM AND CALCULATE THEIR SIZE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { GET_INFO  } from "$projectDir/modules/local/module_preprocessing/get_info.nf"

workflow POSTPROCESSING {
    // take:

    main:
        println "Get infos from genomes"
        // GET_INFO // download lineage, mnemonic, genome size
        // CLEAN_GENOME_HEADERS
    // emit:
}

