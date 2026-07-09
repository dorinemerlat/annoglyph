/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nextflow workflow to annotate non-coding genes using Infernal, Barrnap and RNAmmer.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { DOWNLOAD_RFAM } from "$projectDir/modules/local/module_ncgenes/download_rfam.nf"
// include { INFERNAL } from "$projectDir/modules/local/module_ncgenes/infernal.nf"
include { BARRNAP } from "$projectDir/modules/local/module_ncgenes/barrnap.nf"
include { RNAMMER } from "$projectDir/modules/local/module_ncgenes/rnammer.nf"
include { TRNASCANSE } from "$projectDir/modules/local/module_ncgenes/trnascanse.nf"
include { INFERNAL } from "$projectDir/modules/local/module_ncgenes/infernal.nf"

workflow NCGENES_ANNOTATION {
    take:
        genomes

    main:
        println "Run non-coding genes annotation"
        // DOWNLOAD_RFAM()

        // INFERNAL(genomes.combine(DOWNLOAD_RFAM.out.models))

        genomes.map { id, meta, genome -> [id, meta, genome, 'euk']}
            .concat ( genomes.map { id, meta, genome -> [id, meta, genome, 'mito']})
            .set { barrnap_inputs }

        BARRNAP(barrnap_inputs)
        // RNAMMER(genomes)

        TRNASCANSE(genomes)

        // INFERNAL(genomes)

    // emit:
        //     infernal = INFERNAL.out
        //     barrnap_nucl = BARRNAP_NUCL.out
        //     barrnap_mito = BARRNAP_MITO.out
        //     rnammer = RNAMMER.out
}


