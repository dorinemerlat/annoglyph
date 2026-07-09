include { GENERATE_ANNOTATION_TABLES                                } from "$projectDir/modules/local/module_genes/annotation_comparison_and_selection/generate_annotation_tables.nf"
include { CONCAT_ANNOTATION_TABLES                                  } from "$projectDir/modules/local/module_genes/annotation_comparison_and_selection/concat_annotation_tables.nf"
include { CONCAT_ANNOTATION_TABLES as CONCAT_BEST_ANNOTATION_TABLES } from "$projectDir/modules/local/module_genes/annotation_comparison_and_selection/concat_annotation_tables.nf"
include { SELECT_BEST_ANNOTATION                                    } from "$projectDir/modules/local/module_genes/annotation_comparison_and_selection/select_best_annotation.nf"
include { PLOT_ANNOTATION_COMPARISON                                } from "$projectDir/modules/local/module_genes/annotation_comparison_and_selection/plot_annotation_comparison.nf"
include { PLOT_ANNOTATION_ANALYSIS                                  } from "$projectDir/modules/local/module_genes/annotation_comparison_and_selection/plot_annotation_analysis.nf"


def collect_best_annotation_tables(annotation_tables, best_annotations, filtering_status, label) {
    annotation_tables
        .filter { id, meta, annotation_method, status, summary, protein_lengths, intron_lengths ->
            status == filtering_status
        }
        .combine(best_annotations, by: [0, 2])
        .map { id, meta, annotation_method, status, summary, protein_lengths, intron_lengths, selected_method ->
            tuple(label, summary, protein_lengths, intron_lengths)
        }
        .groupTuple(by: 0)
}


def collect_best_annotation_files(qc_metrics, best_annotations) {
    qc_metrics
        .filter { id, meta, annotation_method, filtering_status, fasta, gff, busco_json, omark_sum ->
            filtering_status == "after-contamination-filtering"
        }
        .combine(best_annotations, by: [0, 2])
        .map { id, meta, annotation_method, filtering_status, fasta, gff, busco_json, omark_sum, selected_method ->
            tuple(id, meta, annotation_method, fasta, gff)
        }
}


workflow ANNOTATION_COMPARISON_AND_SELECTION {

    take:
        /*
         * QC metrics before contamination filtering.
         *
         * Expected format:
         *   [id, meta, annotation_method, filtering_status, fasta, gff, busco_json, omark_sum]
         */
        qc_metrics_before_filtering

        /*
         * QC metrics after contamination filtering.
         *
         * For annotations without OMARK contaminants, this channel should contain
         * the original files relabelled as "after-contamination-filtering".
         *
         * Expected format:
         *   [id, meta, annotation_method, filtering_status, fasta, gff, busco_json, omark_sum]
         */
        qc_metrics_after_filtering

        /*
         * Species tree used for plotting.
         */
        newick

        /*
         * Genome statistics used by PLOT_ANNOTATION_ANALYSIS.
         */
        all_genome_stats

    main:

        /*
         * Generate per-annotation summary, protein length and intron length tables.
         */
        qc_metrics_before_filtering
            .concat(qc_metrics_after_filtering)
            .set { all_qc_metrics }

        GENERATE_ANNOTATION_TABLES(all_qc_metrics)


        /*
         * Concatenate all generated tables by filtering status.
         *
         * Output format expected by CONCAT_ANNOTATION_TABLES:
         *   [filtering_status, [summary_files], [protein_length_files], [intron_length_files]]
         */
        GENERATE_ANNOTATION_TABLES.out
            .map { id, meta, annotation_method, filtering_status, summary, protein_lengths, intron_lengths ->
                tuple(filtering_status, summary, protein_lengths, intron_lengths)
            }
            .groupTuple(by: 0)
            .set { all_tables }

        CONCAT_ANNOTATION_TABLES(all_tables)


        /*
         * Select the best annotation method per species using the merged
         * after-contamination-filtering summary table.
         */
        CONCAT_ANNOTATION_TABLES.out
            .filter { filtering_status, summary, protein_lengths, intron_lengths ->
                filtering_status == "after-contamination-filtering"
            }
            .map { filtering_status, summary, protein_lengths, intron_lengths ->
                tuple(filtering_status, summary)
            }
            .set { summary_for_selection }

        SELECT_BEST_ANNOTATION(summary_for_selection)


        /*
         * Convert the selected best annotation table to a reusable channel.
         *
         * Output format:
         *   [id, null, annotation_method]
         *
         * The null value is only a placeholder so that this channel can be
         * combined with channels whose annotation method is at index 2.
         */
        SELECT_BEST_ANNOTATION.out
            .flatMap { status, file ->
                file.text.readLines()
                    .drop(1)
                    .findAll { it.trim() }
                    .collect { line ->
                        def cols = line.split('\t', -1)
                        tuple(cols[0], null, cols[1])
                    }
            }
            .set { best_annotations }


        /*
         * Collect summary, protein length and intron length tables corresponding
         * only to the selected best annotation of each species.
         */
        collect_best_annotation_tables(
            GENERATE_ANNOTATION_TABLES.out,
            best_annotations,
            "after-contamination-filtering",
            "best-annotations"
        ).set { best_annotation_tables }

        CONCAT_BEST_ANNOTATION_TABLES(best_annotation_tables)


        /*
         * Recover the final best protein FASTA and GFF files.
         */
        collect_best_annotation_files(
            qc_metrics_after_filtering,
            best_annotations
        ).set { best_annotation_files }

        best_annotation_files
            .map { id, meta, annotation_method, fasta, gff ->
                tuple(id, meta, fasta)
            }
            .set { best_proteins }

        best_annotation_files
            .map { id, meta, annotation_method, fasta, gff ->
                tuple(id, meta, gff)
            }
            .set { best_gff }


        /*
         * Plot all annotation sets plus the selected best annotation set.
         */
        CONCAT_BEST_ANNOTATION_TABLES.out
            .concat(CONCAT_ANNOTATION_TABLES.out)
            .combine(newick)
            .set { annotation_tables_for_comparison }

        PLOT_ANNOTATION_COMPARISON(annotation_tables_for_comparison)


        /*
         * Run global annotation analysis only on the selected best annotation
         * set across all species.
         */
        CONCAT_BEST_ANNOTATION_TABLES.out
            .combine(newick)
            .combine(all_genome_stats.map { id, file -> file })
            .set { annotation_tables_for_analysis }

        PLOT_ANNOTATION_ANALYSIS(annotation_tables_for_analysis)


    emit:
        best_annotations = best_annotations
        best_gff         = best_gff
        best_proteins    = best_proteins
}