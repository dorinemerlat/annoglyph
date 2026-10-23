#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { STRUCTURAL_ANNOTATION } from './subworkflows/genes_annotation/structural_annotation'
include { REUSED_REPEATS_ANNOTATION } from './subworkflows/reused_repeats_annotation'

workflow {
    def imported = LegacyCacheUtils.load(params.legacy_cache_manifest, params.cacheDir)
    assert !LegacyCacheUtils.needsRnaseq(imported, 'fixture', 'both')
    assert LegacyCacheUtils.needsRnaseq([:], 'fixture', 'both')
    def partial = [genomes: [fixture: [predictions: ['braker-with-star': imported.genomes.fixture.predictions['braker-with-star']]]]]
    assert !LegacyCacheUtils.needsRnaseq(partial, 'fixture', 'star')
    assert LegacyCacheUtils.needsRnaseq(partial, 'fixture', 'varus')
    assert LegacyCacheUtils.needsRnaseq(partial, 'fixture', 'both')
    def meta = [id: 'fixture', name: 'Fixture species', assembly_size: 8,
                mnemonic: [annoglyph: 'FIXTR', official: 'FIXTR'],
                lineage: [species: [name: 'Fixture species']]]
    def entry = imported.genomes.fixture
    genomes = Channel.of([meta, file(entry.masked)])
    REUSED_REPEATS_ANNOTATION(genomes)
    STRUCTURAL_ANNOTATION(genomes, Channel.empty(), Channel.value(file(entry.family_table)), Channel.value('arthropoda_odb10'))
    STRUCTURAL_ANNOTATION.out.raw_gff
        .map { metadata, gff -> [metadata.id, metadata.annotation_method, gff.name].join('\t') }
        .collectFile(name: 'imported_predictions.tsv', storeDir: params.validation_out, newLine: true, sort: true)
    REUSED_REPEATS_ANNOTATION.out.masked_genomes
        .map { metadata, genome -> metadata.id }
        .collectFile(name: 'imported_repeats.tsv', storeDir: params.validation_out, newLine: true)
}
