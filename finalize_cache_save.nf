nextflow.enable.dsl = 2

// Explicitly finalize historical results. This entry point never calls coding
// predictors, functional annotation, BUSCO, OMArk, cmscan, tRNAscan or barrnap.
include { MERGE_NCRNA } from './modules/local/module_ncgenes/merge_ncrna.nf'
include { NCRNA_STATS } from './modules/local/module_ncgenes/ncrna_stats.nf'
include { PLOT_REPEAT_LANDSCAPE } from './modules/local/repeats_annotation/plot_landscape.nf'
include { POSTPROCESSING } from './subworkflows/postprocessing.nf'

process IMPORT_SELECTED_CODING {
    tag "${meta.id}"
    label 'annoglyph_python'
    memory 8.GB
    input:
    tuple val(meta), path(source), path(reference), path(mapping)
    output:
    tuple val(meta), path('selected_coding.gff3'), emit: gff
    tuple val(meta), path('coding_import.json'), emit: audit
    script:
    """
    python ${projectDir}/bin/normalize_cached_coding.py --source ${source} \\
        --reference ${reference} --mapping ${mapping} --output selected_coding.gff3
    """
}

workflow {
    def imported = new groovy.json.JsonSlurperClassic().parse(file(params.legacy_cache_manifest).toFile())
    def metadata = new groovy.json.JsonSlurperClassic().parse(file("${projectDir}/outputs/cache_reuse_20261023/taxonomy_validation/taxonomy.json").toFile())
    def selections = new groovy.json.JsonSlurperClassic().parse(file("${projectDir}/outputs/cache_reuse_20261023/final_coding_sources.json").toFile())
    def activeIds = selections.findAll { id, selection -> selection.source }.keySet()
    def taxonomy = Channel.value(file("${projectDir}/outputs/cache_reuse_20261023/taxonomy_validation/taxonomy.tsv", checkIfExists: true))
    def tree = Channel.value(file("${projectDir}/outputs/cache_reuse_20261023/taxonomy_validation/tree.nwk", checkIfExists: true))
    records = Channel.fromList(imported.genomes.findAll { id, entry -> id in activeIds }.collect { id, entry ->
        def meta = metadata[id] + [assembly_size: entry.validation.masked.bases_verified]
        [meta, entry]
    })
    selected = records.filter { meta, entry -> selections[meta.id]?.source }
        .map { meta, entry ->
            def selection = selections[meta.id]
            [meta + [annotation_method: selection.method], file(selection.source, checkIfExists: true),
                file(entry.reference, checkIfExists: true), file(entry.mapping.path, checkIfExists: true)]
        }
    IMPORT_SELECTED_CODING(selected)
    ncinputs = records.map { meta, entry ->
        [meta, file("${entry.trnascanse}/${meta.id}_trna.gff", checkIfExists: true),
            file(entry.barrnap.euk, checkIfExists: true), file(entry.barrnap.mito, checkIfExists: true),
            file("${params.cacheDir}/ncgenes_annotation/infernal/${meta.id}/${meta.id}_infernal.tblout", checkIfExists: true)]
    }
    MERGE_NCRNA(ncinputs, Channel.value(file("${params.cacheDir}/annotation_export_v1/ncgenes_annotation/download_rfam/15.0/rfam", checkIfExists: true)))
    manifest = MERGE_NCRNA.out.gff.map { meta, gff -> [meta.id, meta.assembly_size, gff.name].join('\t') }
        .collectFile(name: 'ncrna_manifest.tsv', seed: 'genome_id\tgenome_size_bp\tannotation\n', newLine: true, sort: true)
    NCRNA_STATS(MERGE_NCRNA.out.gff.map { meta, gff -> gff }.collect(), manifest, taxonomy, tree)
    PLOT_REPEAT_LANDSCAPE(records.map { meta, entry -> [meta, file(entry.landscape, checkIfExists: true)] })
    species_assets = records.flatMap { meta, entry -> [
        [meta.id, [file(entry.family_table, checkIfExists: true), 'repeats/stats/families.tsv']],
        [meta.id, [file(entry.landscape, checkIfExists: true), 'repeats/stats/landscape.tsv']]
    ] }.mix(PLOT_REPEAT_LANDSCAPE.out.plots.flatMap { meta, plots -> plots.collect { p -> [meta.id, [p, "repeats/plots/${p.name}"]] } })
      .mix(MERGE_NCRNA.out.filtering.map { meta, audit -> [meta.id, [audit, 'noncoding/stats/filtering.tsv']] })
      .mix(IMPORT_SELECTED_CODING.out.audit.map { meta, audit -> [meta.id, [audit, 'protein_coding/selection/cache_import.json']] })
    comparative_assets = NCRNA_STATS.out.tables.map { p -> [p, 'noncoding/stats'] }
        .mix(NCRNA_STATS.out.figures.map { p -> [p, 'noncoding/plots'] })
        .mix(NCRNA_STATS.out.readme.map { p -> [p, 'noncoding/README.txt'] })
        .mix(Channel.of([file("${projectDir}/outputs/cache_reuse_20261023/final_coding_sources.json"), 'protein_coding/stats/cache_sources.json']))
    POSTPROCESSING(records.map { meta, entry -> [meta, file(entry.masked, checkIfExists: true)] },
        IMPORT_SELECTED_CODING.out.gff,
        records.map { meta, entry -> [meta, file(entry.repeats_gff, checkIfExists: true)] },
        MERGE_NCRNA.out.gff, species_assets, comparative_assets, taxonomy, tree)
}
