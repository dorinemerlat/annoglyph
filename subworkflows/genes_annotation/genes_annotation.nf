include { RNASEQ_ASSEMBLY } from "$projectDir/subworkflows/genes_annotation/rnaseq_assembly.nf"
include { STRUCTURAL_ANNOTATION } from "$projectDir/subworkflows/genes_annotation/structural_annotation.nf"
include { FUNCTIONAL_ANNOTATION } from "$projectDir/subworkflows/genes_annotation/functional_annotation.nf"
include { ANNOTATION_QC_AND_CONTAMINATION_FILTERING } from "$projectDir/subworkflows/genes_annotation/annotation_qc_and_contamination_filtering.nf"
include { ANNOTATION_COMPARISON_AND_SELECTION } from "$projectDir/subworkflows/genes_annotation/annotation_comparison_and_selection.nf"

workflow GENES_ANNOTATION {
    take:
    preprocessed_genomes
    masked_genomes
    repeat_gff
    taxonomy_table
    taxonomy_tree
    genome_stats
    busco_dataset
    busco_dataset_dir
    main:
    RNASEQ_ASSEMBLY(preprocessed_genomes)
    STRUCTURAL_ANNOTATION(masked_genomes, RNASEQ_ASSEMBLY.out.bams, taxonomy_table, busco_dataset)
    FUNCTIONAL_ANNOTATION(STRUCTURAL_ANNOTATION.out.gff, STRUCTURAL_ANNOTATION.out.proteins)
    ANNOTATION_QC_AND_CONTAMINATION_FILTERING(
        FUNCTIONAL_ANNOTATION.out.gff, STRUCTURAL_ANNOTATION.out.proteins,
        busco_dataset_dir.map { meta, directory -> directory }.first()
    )
    ANNOTATION_COMPARISON_AND_SELECTION(
        ANNOTATION_QC_AND_CONTAMINATION_FILTERING.out.qc_metrics_before_filtering,
        ANNOTATION_QC_AND_CONTAMINATION_FILTERING.out.qc_metrics_after_filtering,
        taxonomy_tree, genome_stats
    )
    emit:
    databases = ANNOTATION_QC_AND_CONTAMINATION_FILTERING.out.database.map { path -> ['oma_luca', path] }.mix(STRUCTURAL_ANNOTATION.out.protein_evidence.map { path -> ['protein_evidence', path] }).mix(STRUCTURAL_ANNOTATION.out.evidence_metadata)
    versions = RNASEQ_ASSEMBLY.out.versions.mix(STRUCTURAL_ANNOTATION.out.versions).mix(FUNCTIONAL_ANNOTATION.out.versions)
    candidate_tables = ANNOTATION_COMPARISON_AND_SELECTION.out.candidate_tables
    selected_qc = ANNOTATION_COMPARISON_AND_SELECTION.out.selected_qc
    comparative_assets = ANNOTATION_COMPARISON_AND_SELECTION.out.comparative_assets
    gff = ANNOTATION_COMPARISON_AND_SELECTION.out.best_gff
    proteins = ANNOTATION_COMPARISON_AND_SELECTION.out.best_proteins
    best_annotations = ANNOTATION_COMPARISON_AND_SELECTION.out.best_annotations
    tables = ANNOTATION_COMPARISON_AND_SELECTION.out.tables
    method_gff = FUNCTIONAL_ANNOTATION.out.gff
    method_proteins = STRUCTURAL_ANNOTATION.out.proteins
}
