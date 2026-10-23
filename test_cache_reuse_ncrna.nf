#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { TRNASCANSE } from './modules/local/module_ncgenes/trnascanse'
include { BARRNAP as BARRNAP_EUK; BARRNAP as BARRNAP_MITO } from './modules/local/module_ncgenes/barrnap'
include { INFERNAL } from './modules/local/module_ncgenes/infernal'

workflow {
    def imported = LegacyCacheUtils.load(params.legacy_cache_manifest, params.cacheDir)
    genomes = Channel.fromList(imported.genomes.collect { id, entry ->
        [[id: id, assembly_size: entry.validation.masked.bases_verified,
          mnemonic: [annoglyph: entry.mnemonic]], file(entry.reference)]
    })
    TRNASCANSE(genomes)
    BARRNAP_EUK(genomes.map { meta, genome -> [meta, genome, 'euk'] })
    BARRNAP_MITO(genomes.map { meta, genome -> [meta, genome, 'mito'] })
    INFERNAL(genomes.map { meta, genome -> [meta, genome, file(params.validation_database)] })
    TRNASCANSE.out.gff.map { meta, gff -> 'trna\t' + meta.id }
        .concat(BARRNAP_EUK.out.gff.map { meta, gff -> 'barrnap_euk\t' + meta.id })
        .concat(BARRNAP_MITO.out.gff.map { meta, gff -> 'barrnap_mito\t' + meta.id })
        .concat(INFERNAL.out.tblout.map { meta, tbl -> 'infernal\t' + meta.id })
        .collectFile(name: 'cache_hits.tsv', storeDir: params.validation_out, newLine: true, sort: true)
}
