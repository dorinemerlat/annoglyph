include { GENERATE_ANNOTATION_TABLES } from "$projectDir/modules/local/genes_annotation/annotation_comparison_and_selection/generate_annotation_tables.nf"
include { CONCAT_ANNOTATION_TABLES } from "$projectDir/modules/local/genes_annotation/annotation_comparison_and_selection/concat_annotation_tables.nf"
include { SELECT_BEST_ANNOTATION } from "$projectDir/modules/local/genes_annotation/annotation_comparison_and_selection/select_best_annotation.nf"
include { PLOT_ANNOTATION_COMPARISON } from "$projectDir/modules/local/genes_annotation/annotation_comparison_and_selection/plot_annotation_comparison.nf"
include { PLOT_ANNOTATION_ANALYSIS } from "$projectDir/modules/local/genes_annotation/annotation_comparison_and_selection/plot_annotation_analysis.nf"

workflow ANNOTATION_COMPARISON_AND_SELECTION {
    take:
    before_metrics
    after_metrics
    newick
    genome_stats
    taxonomy
    main:
    GENERATE_ANNOTATION_TABLES(before_metrics.concat(after_metrics))
    all_tables = GENERATE_ANNOTATION_TABLES.out
        .map { meta, summary, lengths, introns -> [meta.filtering_status, summary, lengths, introns] }
        .groupTuple(by: 0)
    CONCAT_ANNOTATION_TABLES(all_tables)
    SELECT_BEST_ANNOTATION(CONCAT_ANNOTATION_TABLES.out
        .filter { status, summary, lengths, introns -> status == 'after-contamination-filtering' }
        .map { status, summary, lengths, introns -> [status, summary] })

    selected_keys = SELECT_BEST_ANNOTATION.out.selection.flatMap { status, selection ->
        selection.readLines().drop(1).findAll { it.trim() }.collect { line ->
            def cols = line.split('\t', -1)
            [[cols[0], cols[1]], true]
        }
    }
    selected_files = after_metrics
        .map { meta, fasta, annotation, busco, omark -> [[meta.id, meta.annotation_method], meta, fasta, annotation, busco, omark] }
        .join(selected_keys, failOnDuplicate: true, failOnMismatch: false)
        .map { key, meta, fasta, annotation, busco, omark, selected -> [meta, fasta, annotation, busco, omark] }

    comparison_assets = Channel.empty()
    if (params.comparative_genomics) {
        PLOT_ANNOTATION_COMPARISON(CONCAT_ANNOTATION_TABLES.out.combine(newick.first()).combine(taxonomy.first()))
        analysis_inputs = CONCAT_ANNOTATION_TABLES.out
            .filter { status, summary, lengths, introns -> status == 'after-contamination-filtering' }
            .map { status, summary, lengths, introns -> [status, summary, introns] }
            .combine(newick.first())
            .combine(genome_stats.map { meta, stats -> stats }.first())
            .combine(taxonomy.first())
            .combine(SELECT_BEST_ANNOTATION.out.selection.map { status, selection -> selection }.first())
        PLOT_ANNOTATION_ANALYSIS(analysis_inputs)
        comparison_assets = PLOT_ANNOTATION_COMPARISON.out.comparison
            .flatMap { status, files -> (files instanceof List ? files : [files]).collect { file -> [file, "protein_coding/plots/comparison/${status}/${file.name}"] } }
            .mix(PLOT_ANNOTATION_ANALYSIS.out.structure.flatMap { status, files -> (files instanceof List ? files : [files]).collect { file -> [file, "protein_coding/plots/structure/${file.name}"] } })
            .mix(PLOT_ANNOTATION_ANALYSIS.out.correlations.flatMap { status, files -> (files instanceof List ? files : [files]).collect { file -> [file, "protein_coding/plots/correlations/${file.name}"] } })
    }
    emit:
    candidate_tables = GENERATE_ANNOTATION_TABLES.out
    selected_qc = selected_files.map { meta, fasta, annotation, busco, omark -> [meta, [busco, omark]] }
    comparative_assets = comparison_assets
        .mix(SELECT_BEST_ANNOTATION.out.scores.map { file -> [file, 'protein_coding/stats/selection_scores.tsv'] })
        .mix(SELECT_BEST_ANNOTATION.out.thresholds.map { file -> [file, 'protein_coding/stats/selection_thresholds.json'] })
    best_gff = selected_files.map { meta, fasta, annotation, busco, omark -> [meta, annotation] }
    best_proteins = selected_files.map { meta, fasta, annotation, busco, omark -> [meta, fasta] }
    best_annotations = SELECT_BEST_ANNOTATION.out.selection
    tables = CONCAT_ANNOTATION_TABLES.out
}
