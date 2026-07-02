/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNASEQ_ASSEMBLY                           } from "$projectDir/subworkflows/local/genes_annotation/rnaseq_assembly.nf"
include { STRUCTURAL_ANNOTATION                     } from "$projectDir/subworkflows/local/genes_annotation/structural_annotation.nf"
include { FUNCTIONAL_ANNOTATION                     } from "$projectDir/subworkflows/local/genes_annotation/functionnal_annotation.nf"
include { ANNOTATION_QC_AND_CONTAMINATION_FILTERING } from "$projectDir/subworkflows/local/genes_annotation/annotation_qc_and_contamination_filtering.nf"
include { ANNOTATION_COMPARISON_AND_SELECTION       } from "$projectDir/subworkflows/local/genes_annotation/annotation_comparison_and_selection.nf"
include { COMPLEMENTARY_PROTEIN_ANALYSES            } from "$projectDir/subworkflows/local/genes_annotation/complementary_protein_analyses.nf"


workflow GENES_ANNOTATION {

    take:
        /*
         * Repeat-masked genomes used for structural gene annotation.
         *
         * Expected format:
         *   [id, meta, fasta]
         */
        masked_genomes

        /*
         * Repeat annotation GFF.
         *
         * Expected format:
         *   [id, meta, repeat_gff]
         */
        repeat_gff

        /*
         * Species tree.
         */
        newick

        /*
         * Genome statistics.
         */
        all_genome_stats


    main:

        /*
         * Assemble and align RNA-seq data.
         */
        RNASEQ_ASSEMBLY(masked_genomes)


        // /*
        //  * Structural annotation.
        //  */
        // STRUCTURAL_ANNOTATION(
        //     masked_genomes,
        //     RNASEQ_ASSEMBLY.out.bams
        // )


        // /*
        //  * Functional annotation.
        //  */
        // FUNCTIONAL_ANNOTATION(
        //     STRUCTURAL_ANNOTATION.out.gff,
        //     STRUCTURAL_ANNOTATION.out.proteins
        // )


        // /*
        //  * Compute annotation quality metrics and perform OMARK-based
        //  * contamination filtering.
        //  */
        // ANNOTATION_QC_AND_CONTAMINATION_FILTERING(
        //     FUNCTIONAL_ANNOTATION.out.gff,
        //     STRUCTURAL_ANNOTATION.out.proteins
        // )


        // /*
        //  * Compare all annotation strategies and select the best annotation
        //  * for each species.
        //  */
        // ANNOTATION_COMPARISON_AND_SELECTION(
        //     ANNOTATION_QC_AND_CONTAMINATION_FILTERING.out.qc_metrics_before_filtering,
        //     ANNOTATION_QC_AND_CONTAMINATION_FILTERING.out.qc_metrics_after_filtering,
        //     newick,
        //     all_genome_stats
        // )


        // /*
        //  * Run complementary analyses on the selected proteomes.
        //  */
        // COMPLEMENTARY_PROTEIN_ANALYSES(
        //     ANNOTATION_COMPARISON_AND_SELECTION.out.best_proteins
        // )


    // emit:

    //     /*
    //      * Best annotation selected for each species.
    //      */
    //     best_annotations = ANNOTATION_COMPARISON_AND_SELECTION.out.best_annotations

    //     /*
    //      * Best GFF files.
    //      */
    //     best_gff = ANNOTATION_COMPARISON_AND_SELECTION.out.best_gff

    //     /*
    //      * Best canonical protein FASTA files.
    //      */
    //     best_proteins = ANNOTATION_COMPARISON_AND_SELECTION.out.best_proteins
}