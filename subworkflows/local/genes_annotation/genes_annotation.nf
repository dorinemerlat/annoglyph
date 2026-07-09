include { RNASEQ_ASSEMBLY                           } from "$projectDir/subworkflows/local/genes_annotation/rnaseq_assembly.nf"
include { STRUCTURAL_ANNOTATION                     } from "$projectDir/subworkflows/local/genes_annotation/structural_annotation.nf"
include { FUNCTIONAL_ANNOTATION                     } from "$projectDir/subworkflows/local/genes_annotation/functionnal_annotation.nf"
include { ANNOTATION_QC_AND_CONTAMINATION_FILTERING } from "$projectDir/subworkflows/local/genes_annotation/annotation_qc_and_contamination_filtering.nf"
include { ANNOTATION_COMPARISON_AND_SELECTION       } from "$projectDir/subworkflows/local/genes_annotation/annotation_comparison_and_selection.nf"

workflow GENES_ANNOTATION {
    take:
        masked_genomes
        repeats_gff
        newick
        all_genome_stats

    main:
        RNASEQ_ASSEMBLY(masked_genomes)
        STRUCTURAL_ANNOTATION(masked_genomes, RNASEQ_ASSEMBLY.out.bams)

        FUNCTIONAL_ANNOTATION(
            STRUCTURAL_ANNOTATION.out.gff,
            STRUCTURAL_ANNOTATION.out.proteins
        )

        ANNOTATION_QC_AND_CONTAMINATION_FILTERING(
            FUNCTIONAL_ANNOTATION.out.gff,
            STRUCTURAL_ANNOTATION.out.proteins
        )

        ANNOTATION_COMPARISON_AND_SELECTION(
            ANNOTATION_QC_AND_CONTAMINATION_FILTERING.out.qc_metrics_before_filtering,
            ANNOTATION_QC_AND_CONTAMINATION_FILTERING.out.qc_metrics_after_filtering,
            newick,
            all_genome_stats
        )

//         PUBLISH_BEST_ANNOTATIONS(
//             ANNOTATION_COMPARISON_AND_SELECTION.out.best_gff
//                 .concat(ANNOTATION_COMPARISON_AND_SELECTION.out.best_proteins)
//         )

        COMPLEMENTARY_PROTEIN_ANALYSES(
            ANNOTATION_COMPARISON_AND_SELECTION.out.best_proteins
        )

//     emit:
//         best_gff      = ANNOTATION_COMPARISON_AND_SELECTION.out.best_gff
//         best_proteins = ANNOTATION_COMPARISON_AND_SELECTION.out.best_proteins
//         tables        = ANNOTATION_COMPARISON_AND_SELECTION.out.tables
}