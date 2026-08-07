include { GENERATE_ANNOTATION_TABLES                                } from "$projectDir/modules/genes_annotation/annotation_comparison_and_selection/generate_annotation_tables.nf"
include { CONCAT_ANNOTATION_TABLES                                  } from "$projectDir/modules/genes_annotation/annotation_comparison_and_selection/concat_annotation_tables.nf"
include { CONCAT_ANNOTATION_TABLES as CONCAT_BEST_ANNOTATION_TABLES } from "$projectDir/modules/genes_annotation/annotation_comparison_and_selection/concat_annotation_tables.nf"
include { SELECT_BEST_ANNOTATION                                    } from "$projectDir/modules/genes_annotation/annotation_comparison_and_selection/select_best_annotation.nf"
include { PLOT_ANNOTATION_COMPARISON                                } from "$projectDir/modules/genes_annotation/annotation_comparison_and_selection/plot_annotation_comparison.nf"
include { PLOT_ANNOTATION_ANALYSIS                                  } from "$projectDir/modules/genes_annotation/annotation_comparison_and_selection/plot_annotation_analysis.nf"
include { ANNOTATE_PROTEIN_FASTA                                    } from "$projectDir/modules/genes_annotation/annotation_comparison_and_selection/annotate_protein_fasta.nf"


def collect_best_annotation_tables(annotation_tables, best_annotations, label) {
    annotation_tables
        .combine(best_annotations, by: [0, 2])
        .map { id, annotation_method, meta, status, summary, protein_lengths, intron_lengths, selected_method ->
            [label, summary, protein_lengths, intron_lengths]
        }
        .groupTuple(by: 0)
}


def collect_best_annotation_files(qc_metrics, best_annotations) {
    qc_metrics
        .combine(best_annotations, by: [0, 2])
        .map { id, annotation_method, meta, filtering_status, fasta, gff, busco_json, omark_sum, selected_method ->
            [id, meta, annotation_method, fasta, gff]
        }
}


workflow ANNOTATION_COMPARISON_AND_SELECTION {

    take:
        qc_metrics_before_filtering
        qc_metrics_after_filtering
        newick
        all_genome_stats

    main:
        print "Running module annotation comparison and selection"

        // Generate annotation tables before and after filtering.
        qc_metrics_before_filtering
            .concat(qc_metrics_after_filtering)
            .set { qc_metrics_for_table_generation }

        GENERATE_ANNOTATION_TABLES(qc_metrics_for_table_generation)

        GENERATE_ANNOTATION_TABLES.out
            .filter { id, meta, annotation_method, filtering_status, summary, protein_lengths, intron_lengths ->
                filtering_status == 'before-contamination-filtering'
            }
            .set { annotation_tables_before_filtering }

        GENERATE_ANNOTATION_TABLES.out
            .filter { id, meta, annotation_method, filtering_status, summary, protein_lengths, intron_lengths ->
                filtering_status == 'after-contamination-filtering'
            }
            .set { annotation_tables_after_filtering }

        // Concatenate tables by filtering status.
        annotation_tables_before_filtering
            .concat(annotation_tables_after_filtering)
            .map { id, meta, annotation_method, filtering_status, summary, protein_lengths, intron_lengths ->
                [filtering_status, summary, protein_lengths, intron_lengths]
            }
            .groupTuple(by: 0)
            .set { all_tables }

        CONCAT_ANNOTATION_TABLES(all_tables)

        // Select the best annotation for each species.
        CONCAT_ANNOTATION_TABLES.out
            .filter { filtering_status, summary, protein_lengths, intron_lengths ->
                filtering_status == 'after-contamination-filtering'
            }
            .map { filtering_status, summary, protein_lengths, intron_lengths ->
                [filtering_status, summary]
            }
            .set { summary_for_selection }

        SELECT_BEST_ANNOTATION(summary_for_selection)

        // Convert the selection table to a channel.
        SELECT_BEST_ANNOTATION.out
            .flatMap { status, file ->
                file.text.readLines()
                    .drop(1)
                    .findAll { it.trim() }
                    .collect { line ->
                        def cols = line.split('\t', -1)
                        [cols[0], null, cols[1]]
                    }
            }
            .set { best_annotations }

        // Collect tables for selected annotations.
        collect_best_annotation_tables(
            annotation_tables_after_filtering,
            best_annotations,
            'best-annotations'
        ).set { best_annotation_tables }

        CONCAT_BEST_ANNOTATION_TABLES(best_annotation_tables)

        // Recover selected FASTA and GFF files.
        collect_best_annotation_files(
            qc_metrics_after_filtering,
            best_annotations
        ).set { best_annotation_files }

        best_annotation_files
            .map { id, meta, annotation_method, fasta, gff ->
                [id, meta, fasta]
            }
            .set { best_proteins }

        best_annotation_files
            .map { id, meta, annotation_method, fasta, gff ->
                [id, meta, gff]
            }
            .set { best_gff }

        // Plot annotation comparisons.
        CONCAT_BEST_ANNOTATION_TABLES.out
            .concat(CONCAT_ANNOTATION_TABLES.out)
            .combine(newick)
            .set { annotation_tables_for_comparison }

        PLOT_ANNOTATION_COMPARISON(annotation_tables_for_comparison)

        // Analyze the selected annotations.
        CONCAT_BEST_ANNOTATION_TABLES.out
            .combine(newick)
            .combine(all_genome_stats.map { id, file -> file })
            .set { annotation_tables_for_analysis }

        PLOT_ANNOTATION_ANALYSIS(annotation_tables_for_analysis)

        best_gff
            .join(best_proteins, by: [0, 1])
            .map { id, meta, gff, proteins ->
                [id, meta, gff, proteins]
            }
            .set { best_gff_and_proteins }

        ANNOTATE_PROTEIN_FASTA(best_gff_and_proteins)

    // emit:
    //     best_annotations = best_annotations
    //     best_gff         = best_gff
    //     best_proteins    = best_proteins
}