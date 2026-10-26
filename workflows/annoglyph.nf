include { PREPROCESSING } from "$projectDir/subworkflows/preprocessing.nf"
include { REUSED_REPEATS_ANNOTATION } from "$projectDir/subworkflows/reused_repeats_annotation.nf"
include { REPEATS_ANNOTATION } from "$projectDir/subworkflows/repeats_annotation.nf"
include { GENES_ANNOTATION } from "$projectDir/subworkflows/genes_annotation/genes_annotation.nf"
include { NCGENES_ANNOTATION } from "$projectDir/subworkflows/ncgenes_annotation.nf"
include { POSTPROCESSING } from "$projectDir/subworkflows/postprocessing.nf"

include { COLLECT_PROVENANCE } from "$projectDir/modules/local/postprocessing/provenance.nf"

workflow AnnoGlyph {
    take:
    genomes
    main:
    PREPROCESSING(genomes)
    software_versions = PREPROCESSING.out.versions
    database_artifacts = PREPROCESSING.out.busco_dataset_dir.map { meta, directory -> ["busco_${meta.id}", directory] }
    species_assets = PREPROCESSING.out.assembly_stats.flatMap { meta, files -> ResultsLayoutUtils.assets(meta, files, 'genome/stats') }
        .mix(PREPROCESSING.out.genome_mapping.flatMap { meta, files -> ResultsLayoutUtils.assets(meta, files, 'genome/stats') })
        .mix(PREPROCESSING.out.busco_summaries.flatMap { meta, files -> ResultsLayoutUtils.assets(meta, files, 'genome/stats/busco') })
    comparative_assets = PREPROCESSING.out.genome_stats.map { meta, stats -> [stats, "genome/stats/${stats.name}"] }
    if (params.module_repeats) {
        if (params.legacy_cache_manifest) {
            REUSED_REPEATS_ANNOTATION(PREPROCESSING.out.preprocessed_genomes)
            repeats_outputs = REUSED_REPEATS_ANNOTATION.out
        } else {
            REPEATS_ANNOTATION(PREPROCESSING.out.preprocessed_genomes,
                PREPROCESSING.out.taxonomy_tree, PREPROCESSING.out.genome_stats,
                PREPROCESSING.out.busco_dataset_dir)
            repeats_outputs = REPEATS_ANNOTATION.out
        }
        masked_genomes = repeats_outputs.masked_genomes
        repeat_gff = repeats_outputs.repeat_gff
        software_versions = software_versions.mix(repeats_outputs.versions)
        species_assets = species_assets
            .mix(repeats_outputs.consensus_libraries.flatMap { meta, files -> ResultsLayoutUtils.assets(meta, files, 'repeats/sequences') })
            .mix(repeats_outputs.families_table.flatMap { meta, files -> ResultsLayoutUtils.assets(meta, files, 'repeats/stats') })
            .mix(repeats_outputs.landscapes.flatMap { meta, files -> ResultsLayoutUtils.assets(meta, files, 'repeats/stats') })
            .mix(repeats_outputs.landscape_plots.flatMap { meta, files -> ResultsLayoutUtils.assets(meta, files, 'repeats/plots') })
        comparative_assets = comparative_assets.mix(repeats_outputs.comparative_tables.map { meta, table -> [table, "repeats/stats/${table.name}"] })
    } else {
        masked_genomes = PREPROCESSING.out.preprocessed_genomes
        repeat_gff = Channel.empty()
    }
    if (params.module_genes) {
        GENES_ANNOTATION(PREPROCESSING.out.preprocessed_genomes, masked_genomes, repeat_gff,
            PREPROCESSING.out.taxonomy_table, PREPROCESSING.out.taxonomy_tree,
            PREPROCESSING.out.genome_stats, PREPROCESSING.out.busco_dataset,
            PREPROCESSING.out.busco_dataset_dir)
        coding_gff = GENES_ANNOTATION.out.gff
        software_versions = software_versions.mix(GENES_ANNOTATION.out.versions)
        database_artifacts = database_artifacts.mix(GENES_ANNOTATION.out.databases)
        species_assets = species_assets
            .mix(GENES_ANNOTATION.out.candidate_tables.flatMap { meta, summary, proteins, introns -> ResultsLayoutUtils.assets(meta, [summary, proteins, introns], "protein_coding/selection/candidates/${meta.annotation_method}/${meta.filtering_status}") })
            .mix(GENES_ANNOTATION.out.selected_qc.flatMap { meta, files -> ResultsLayoutUtils.assets(meta, files, 'protein_coding/stats/quality') })
        comparative_assets = comparative_assets.mix(GENES_ANNOTATION.out.comparative_assets)
            .mix(GENES_ANNOTATION.out.tables.flatMap { status, summary, proteins, introns -> [summary, proteins, introns].collect { file -> [file, "protein_coding/stats/candidates/${status}/${file.name}"] } })
            .mix(GENES_ANNOTATION.out.best_annotations.map { status, table -> [table, 'protein_coding/stats/selection_details.tsv'] })
    } else {
        coding_gff = Channel.empty()
    }
    if (params.module_ncgenes) {
        // Use unmasked/preprocessed sequence to avoid suppressing genuine ncRNAs.
        NCGENES_ANNOTATION(PREPROCESSING.out.preprocessed_genomes,
            PREPROCESSING.out.taxonomy_table, PREPROCESSING.out.taxonomy_tree)
        ncrna_gff = NCGENES_ANNOTATION.out.gff
        software_versions = software_versions.mix(NCGENES_ANNOTATION.out.versions)
        database_artifacts = database_artifacts.mix(NCGENES_ANNOTATION.out.database)
        species_assets = species_assets.mix(NCGENES_ANNOTATION.out.filtering.flatMap { meta, files -> ResultsLayoutUtils.assets(meta, files, 'noncoding/stats') })
        comparative_assets = comparative_assets.mix(NCGENES_ANNOTATION.out.tables.map { folder -> [folder, 'noncoding/stats'] })
            .mix(NCGENES_ANNOTATION.out.figures.map { folder -> [folder, 'noncoding/plots'] })
            .mix(NCGENES_ANNOTATION.out.metadata_resolution.map { report -> [report, 'noncoding/stats/rfam_metadata_resolution.json'] })
    } else {
        ncrna_gff = Channel.empty()
    }
    POSTPROCESSING(masked_genomes, coding_gff, repeat_gff, ncrna_gff,
        species_assets, comparative_assets, PREPROCESSING.out.taxonomy_table, PREPROCESSING.out.taxonomy_tree)
    provenance_inputs = genomes.flatMap { meta, genome ->
        def reads = ((meta['RNASeq-file']?.single_end ?: []) + (meta['RNASeq-file']?.paired_end ?: [])).flatten()
        [[meta.id, genome]] + reads.collect { path -> ["${meta.id}/${path.name}", path] }
    }
        .mix(params.legacy_cache_manifest ? Channel.of(['legacy_cache_manifest', file(params.legacy_cache_manifest, checkIfExists: true)]) : Channel.empty())
        .mix(Channel.of(['samplesheet', file(params.input, checkIfExists: true)], ['modules_lock', file("${projectDir}/modules.json", checkIfExists: true)]))
        .mix(Channel.fromList([params.submission_metadata, params.submission_template, params.table2asn_binary].findAll { it }).map { path -> [path.toString(), file(path, checkIfExists: true)] })
        .collect(flat: false)
    used_databases = database_artifacts.collect(flat: false)
    COLLECT_PROVENANCE(ProvenanceUtils.snapshot(params, workflow, nextflow.version),
        provenance_inputs.map { rows -> [rows.collect { it[0] }, rows.collect { it[1].toAbsolutePath().toString() }, rows.collect { it[1] }] },
        software_versions.collect(), used_databases.map { rows -> [rows.collect { it[0] }, rows.collect { it[1] }] },
        params.database_lock_file ? Channel.value(file(params.database_lock_file, checkIfExists: true)) : [],
        Channel.value(file("${projectDir}/bin/collect_provenance.py", checkIfExists: true)),
        POSTPROCESSING.out.comparative.collect())
    emit:
    final_annotations = POSTPROCESSING.out.bundles
}
