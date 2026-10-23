include { DOWNLOAD_RFAM; RFAM_CMPRESS } from "$projectDir/modules/local/module_ncgenes/download_rfam.nf"
include { TRNASCANSE } from "$projectDir/modules/nf-core/trnascanse/main.nf"
include { INFERNAL_CMSCAN as INFERNAL } from "$projectDir/modules/nf-core/infernal/cmscan/main.nf"
include { NORMALIZE_INFERNAL_OUTPUTS } from "$projectDir/modules/local/module_ncgenes/normalize_infernal.nf"
include { BARRNAP as BARRNAP_EUK; BARRNAP as BARRNAP_MITO } from "$projectDir/modules/nf-core/barrnap/main.nf"
include { MERGE_NCRNA } from "$projectDir/modules/local/module_ncgenes/merge_ncrna.nf"

include { NCRNA_STATS } from "$projectDir/modules/local/module_ncgenes/ncrna_stats.nf"

workflow NCGENES_ANNOTATION {
    take:
    genomes
    taxonomy_table
    taxonomy_tree
    main:
    if (params.rfam_db) {
        rfam = file(params.rfam_db, checkIfExists: true)
        ['Rfam.cm', 'Rfam.clanin', 'family.txt', 'Rfam.cm.i1f', 'Rfam.cm.i1i', 'Rfam.cm.i1m', 'Rfam.cm.i1p'].each { name ->
            if (!rfam.resolve(name).exists()) error "Missing Rfam database file: ${rfam}/${name}"
        }
        rfam_models = Channel.value(rfam)
    } else {
        DOWNLOAD_RFAM(Channel.value(params.rfam_release))
        RFAM_CMPRESS(DOWNLOAD_RFAM.out.models)
        rfam_models = RFAM_CMPRESS.out.models.first()
    }
    TRNASCANSE(genomes, true, true, true, false)
    BARRNAP_EUK(genomes.map { meta, genome -> [meta, genome, 'euk'] })
    BARRNAP_MITO(genomes.map { meta, genome -> [meta, genome, 'mito'] })
    INFERNAL(genomes, rfam_models.map { directory -> [[id: 'rfam'], directory] }, rfam_models.map { directory -> directory.resolve('Rfam.clanin') })
    NORMALIZE_INFERNAL_OUTPUTS(INFERNAL.out.tbl.join(INFERNAL.out.txt, failOnDuplicate: true, failOnMismatch: true))
    merged_inputs = TRNASCANSE.out.gff
        .join(BARRNAP_EUK.out.gff, failOnDuplicate: true, failOnMismatch: true)
        .join(BARRNAP_MITO.out.gff, failOnDuplicate: true, failOnMismatch: true)
        .join(NORMALIZE_INFERNAL_OUTPUTS.out.tblout, failOnDuplicate: true, failOnMismatch: true)
    MERGE_NCRNA(merged_inputs, rfam_models)
    statistics_manifest = MERGE_NCRNA.out.gff
        .map { meta, gff -> [meta.id, meta.assembly_size, gff.name].join('\t') }
        .collectFile(name: 'ncrna_manifest.tsv', seed: 'genome_id\tgenome_size_bp\tannotation\n', newLine: true, sort: true)
    statistics_annotations = MERGE_NCRNA.out.gff.map { meta, gff -> gff }.collect()
    NCRNA_STATS(statistics_annotations, statistics_manifest, taxonomy_table.first(), taxonomy_tree.first())
    emit:
    versions = TRNASCANSE.out.versions.mix(INFERNAL.out.versions).mix(BARRNAP_EUK.out.versions).mix(BARRNAP_MITO.out.versions)
    database = rfam_models.map { path -> ['rfam', path] }
    tables = NCRNA_STATS.out.tables
    figures = NCRNA_STATS.out.figures
    gff = MERGE_NCRNA.out.gff
    filtering = MERGE_NCRNA.out.filtering
}
