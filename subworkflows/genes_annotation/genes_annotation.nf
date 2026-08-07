include { RNASEQ_ASSEMBLY                           } from "$projectDir/subworkflows/genes_annotation/rnaseq_assembly.nf"
include { STRUCTURAL_ANNOTATION                     } from "$projectDir/subworkflows/genes_annotation/structural_annotation.nf"
include { FUNCTIONAL_ANNOTATION                     } from "$projectDir/subworkflows/genes_annotation/functional_annotation.nf"
// include { ANNOTATION_QC_AND_CONTAMINATION_FILTERING } from "$projectDir/subworkflows/genes_annotation/annotation_qc_and_contamination_filtering.nf"
// include { ANNOTATION_COMPARISON_AND_SELECTION       } from "$projectDir/subworkflows/genes_annotation/annotation_comparison_and_selection.nf"

workflow GENES_ANNOTATION {
    take:
    preprocessed_genomes
    masked_genomes
    repeat_gff
    taxonomy_table
    taxonomy_tree
    genome_stats
    busco_dataset

    main:
        RNASEQ_ASSEMBLY(preprocessed_genomes)

        STRUCTURAL_ANNOTATION(
            masked_genomes, 
            RNASEQ_ASSEMBLY.out.bams,
            taxonomy_table,
            busco_dataset
        )

        FUNCTIONAL_ANNOTATION(
            STRUCTURAL_ANNOTATION.out.gff,
            STRUCTURAL_ANNOTATION.out.proteins
        )

        ANNOTATION_QC_AND_CONTAMINATION_FILTERING(
            FUNCTIONAL_ANNOTATION.out.gff,
            STRUCTURAL_ANNOTATION.out.proteins
        )

        // ANNOTATION_COMPARISON_AND_SELECTION(
        //     ANNOTATION_QC_AND_CONTAMINATION_FILTERING.out.qc_metrics_before_filtering,
        //     ANNOTATION_QC_AND_CONTAMINATION_FILTERING.out.qc_metrics_after_filtering_filtered_only,
        //     newick,
        //     genome_stats
        // )

//         PUBLISH_BEST_ANNOTATIONS(
//             ANNOTATION_COMPARISON_AND_SELECTION.out.best_gff
//                 .concat(ANNOTATION_COMPARISON_AND_SELECTION.out.best_proteins)
//         )

        // COMPLEMENTARY_PROTEIN_ANALYSES(
        //     ANNOTATION_COMPARISON_AND_SELECTION.out.best_proteins
        // )

//     emit:
//         best_gff      = ANNOTATION_COMPARISON_AND_SELECTION.out.best_gff
//         best_proteins = ANNOTATION_COMPARISON_AND_SELECTION.out.best_proteins
//         tables        = ANNOTATION_COMPARISON_AND_SELECTION.out.tables
}