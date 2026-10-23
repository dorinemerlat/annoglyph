include { PLOT_REPEAT_LANDSCAPE } from "$projectDir/modules/local/repeats_annotation/plot_landscape.nf"
include { MERGE_TABLES as MERGE_FAMILY_TABLES; MERGE_TABLES as MERGE_LANDSCAPES } from "$projectDir/modules/local/shared_modules/merge_tables.nf"

/* Resume a validated historical repeat experiment without inventing missing
 * RepeatModeler/MCHelper logs, alignments or software-version records. */
workflow REUSED_REPEATS_ANNOTATION {
    take:
    genomes

    main:
    def imported = LegacyCacheUtils.load(params.legacy_cache_manifest, params.cacheDir)
    records = genomes.map { meta, genome ->
        def entry = LegacyCacheUtils.genome(imported, meta.id)
        if (meta.mnemonic.annoglyph != entry.mnemonic) error "Imported cache mnemonic mismatch for ${meta.id}"
        [meta, entry]
    }
    masked = records.map { meta, entry -> [meta, file(entry.masked, checkIfExists: true)] }
    gff = records.map { meta, entry -> [meta, file(entry.repeats_gff, checkIfExists: true)] }
    families = records.map { meta, entry -> [meta, file(entry.family_table, checkIfExists: true)] }
    landscapes = records.map { meta, entry -> [meta, file(entry.landscape, checkIfExists: true)] }
    libraries = records.map { meta, entry -> [meta, file(imported.repeat_library, checkIfExists: true)] }
    PLOT_REPEAT_LANDSCAPE(landscapes)
    MERGE_FAMILY_TABLES(families.map { meta, table -> ['all', table] }.groupTuple()
        .map { id, tables -> [[id: id], tables, 'all_families.tsv'] })
    MERGE_LANDSCAPES(landscapes.map { meta, table -> ['all', table] }.groupTuple()
        .map { id, tables -> [[id: id], tables, 'all_repeatlandscape.tsv'] })

    emit:
    versions = Channel.empty()
    consensus_libraries = libraries
    families_table = families
    landscapes = landscapes
    landscape_plots = PLOT_REPEAT_LANDSCAPE.out.plots
    comparative_tables = MERGE_FAMILY_TABLES.out.mix(MERGE_LANDSCAPES.out)
    masked_genomes = masked
    repeat_gff = gff
}
