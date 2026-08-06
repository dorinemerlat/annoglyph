include { GENERATE_ANNOTATION_TABLES                                } from "$projectDir/modules/genes_annotation/annotation_comparison_and_selection/generate_annotation_tables.nf"
include { CONCAT_ANNOTATION_TABLES                                  } from "$projectDir/modules/genes_annotation/annotation_comparison_and_selection/concat_annotation_tables.nf"
include { CONCAT_ANNOTATION_TABLES as CONCAT_BEST_ANNOTATION_TABLES } from "$projectDir/modules/genes_annotation/annotation_comparison_and_selection/concat_annotation_tables.nf"
include { SELECT_BEST_ANNOTATION                                    } from "$projectDir/modules/genes_annotation/annotation_comparison_and_selection/select_best_annotation.nf"
include { PLOT_ANNOTATION_COMPARISON                                } from "$projectDir/modules/genes_annotation/annotation_comparison_and_selection/plot_annotation_comparison.nf"
include { PLOT_ANNOTATION_ANALYSIS                                  } from "$projectDir/modules/genes_annotation/annotation_comparison_and_selection/plot_annotation_analysis.nf"
include { ANNOTATE_PROTEIN_FASTA                                    } from "$projectDir/modules/genes_annotation/annotation_comparison_and_selection/annotate_protein_fasta.nf"

def collect_best_annotation_tables(annotation_tables, best_annotations, filtering_status, label) {
    annotation_tables
        .filter { id, meta, annotation_method, status, summary, protein_lengths, intron_lengths ->
            status == filtering_status
        }
        .combine(best_annotations, by: [0, 2])
        .map { id, annotation_method, meta, status, summary, protein_lengths, intron_lengths, selected_method ->
            [ label, summary, protein_lengths, intron_lengths ]
        }
        .groupTuple(by: 0)
}


def collect_best_annotation_files(qc_metrics, best_annotations) {
    qc_metrics
        .filter { id, meta, annotation_method, filtering_status, fasta, gff, busco_json, omark_sum ->
            filtering_status == "after-contamination-filtering"
        }
        .combine(best_annotations, by: [0, 2])
        .map { id, annotation_method, meta, filtering_status, fasta, gff, busco_json, omark_sum, selected_method ->
            [ id, meta, annotation_method, fasta, gff ]
        }
}


workflow ANNOTATION_COMPARISON_AND_SELECTION {

    take:
        /*
        * Metrics for every annotation before contamination filtering.
        *
        * Format:
        * [id, meta, annotation_method, filtering_status,
        *  fasta, gff, busco_json, omark_sum]
        */
        qc_metrics_before_filtering

        /*
        * Metrics only for annotations that were genuinely filtered.
        *
        * Annotations without contaminants are absent from this channel.
        *
        * Format:
        * [id, meta, annotation_method, filtering_status,
        *  fasta, gff, busco_json, omark_sum]
        */
        qc_metrics_after_filtering_filtered_only

        /*
        * Species tree used for plotting.
        */
        newick

        /*
        * Genome statistics used by PLOT_ANNOTATION_ANALYSIS.
        */
        all_genome_stats

    main:

        print "Running module annotation comparison and selection"

        /*
        * Generate annotation tables for:
        *
        *   1. every annotation before contamination filtering;
        *   2. annotations that were genuinely filtered.
        *
        * Annotations without contaminants do not yet have an
        * "after-contamination-filtering" table at this stage.
        */
        qc_metrics_before_filtering
            .concat(qc_metrics_after_filtering_filtered_only)
            .set { qc_metrics_for_table_generation }

        GENERATE_ANNOTATION_TABLES(qc_metrics_for_table_generation)

        GENERATE_ANNOTATION_TABLES.out
            .set { generated_annotation_tables }


        /*
        * Separate tables generated before filtering from tables generated
        * from genuinely contamination-filtered annotations.
        *
        * Format:
        * [id, meta, annotation_method, filtering_status,
        *  summary, protein_lengths, intron_lengths]
        */
        generated_annotation_tables
            .filter { id, meta, annotation_method, filtering_status, summary, protein_lengths, intron_lengths ->
                filtering_status == 'before-contamination-filtering'
            }
            .set { annotation_tables_before_filtering }

        generated_annotation_tables
            .filter { id, meta, annotation_method, filtering_status, summary, protein_lengths, intron_lengths ->
                filtering_status == 'after-contamination-filtering'
            }
            .set { annotation_tables_after_filtering_filtered_only }


        /*
        * Construct a complete after-contamination-filtering table channel.
        *
        * For annotations with contaminants:
        *   use the tables generated from the filtered FASTA/GFF.
        *
        * For annotations without contaminants:
        *   reuse the before-filtering tables and relabel them as
        *   after-contamination-filtering.
        */
        annotation_tables_before_filtering
            .join(annotation_tables_after_filtering_filtered_only, by: [0, 1, 2], remainder: true)
            .map { row ->

                /*
                * Matching before and after entries:
                * 7 left-side columns + 4 non-key right-side columns = 11.
                */
                if (row.size() == 11) {
                    def (id, meta, annotation_method,
                        before_status, before_summary, before_protein_lengths, before_intron_lengths,
                        after_status, after_summary, after_protein_lengths, after_intron_lengths) = row

                    return [id, meta + [contamination_filtering: true], annotation_method, after_status, after_summary, after_protein_lengths, after_intron_lengths]
                }

                /*
                * No matching filtered annotation:
                * 7 left-side columns + one null remainder value = 8.
                */
                if (row.size() == 8) {
                    def (id, meta, annotation_method,
                        before_status, before_summary, before_protein_lengths, before_intron_lengths,
                        missing_after) = row

                    return [id, meta + [contamination_filtering: false], annotation_method, 'after-contamination-filtering', before_summary, before_protein_lengths, before_intron_lengths]
                }

                throw new IllegalStateException(
                    "Unexpected joined annotation-table row with ${row.size()} elements: ${row}"
                )
            }
            .set { annotation_tables_after_filtering }


        /*
        * Concatenate tables by filtering status.
        *
        * Input expected by CONCAT_ANNOTATION_TABLES:
        * [filtering_status, [summary], [protein_lengths], [intron_lengths]]
        */
        annotation_tables_before_filtering
            .concat(annotation_tables_after_filtering)
            .map { id, meta, annotation_method, filtering_status, summary, protein_lengths, intron_lengths ->
                [filtering_status, summary, protein_lengths, intron_lengths]
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

        // /*
        //  * Collect summary, protein length and intron length tables corresponding
        //  * only to the selected best annotation of each species.
        //  */
        collect_best_annotation_tables(
            annotation_tables_after_filtering,
            best_annotations,
            "after-contamination-filtering",
            "best-annotations"
        ).set { best_annotation_tables }

        CONCAT_BEST_ANNOTATION_TABLES(best_annotation_tables)


        /*
         * Construct a complete after-contamination-filtering QC channel.
         *
         * For annotations with contaminants:
         *   use the genuinely filtered FASTA/GFF and associated QC results.
         *
         * For annotations without contaminants:
         *   reuse the before-filtering files and QC results, but relabel the
         *   filtering status as after-contamination-filtering.
         *
         * Output format:
         * [id, meta, annotation_method, filtering_status,
         *  fasta, gff, busco_json, omark_sum]
         */
        qc_metrics_before_filtering
            .join(
                qc_metrics_after_filtering_filtered_only,
                by: [0, 1, 2],
                remainder: true
            )
            .map { row ->

                /*
                 * Matching before/after entries:
                 * 8 left-side columns + 5 non-key right-side columns = 13.
                 */
                if (row.size() == 13) {
                    def (
                        id,
                        meta,
                        annotation_method,
                        before_status,
                        before_fasta,
                        before_gff,
                        before_busco_json,
                        before_omark_sum,
                        after_status,
                        after_fasta,
                        after_gff,
                        after_busco_json,
                        after_omark_sum
                    ) = row

                    return [ id,
                        meta + [contamination_filtering: true],
                        annotation_method,
                        after_status,
                        after_fasta,
                        after_gff,
                        after_busco_json,
                        after_omark_sum
                    ]
                }

                /*
                 * No matching filtered annotation:
                 * 8 left-side columns + one null remainder value = 9.
                 */
                if (row.size() == 9) {
                    def (
                        id,
                        meta,
                        annotation_method,
                        before_status,
                        before_fasta,
                        before_gff,
                        before_busco_json,
                        before_omark_sum,
                        missing_after
                    ) = row

                    return [ id, 
                        meta + [contamination_filtering: false], 
                        annotation_method, 
                        'after-contamination-filtering',
                        before_fasta, 
                        before_gff, 
                        before_busco_json,
                        before_omark_sum ]
                }
            }
            .set { qc_metrics_after_filtering }

        /*
         * Recover the final best protein FASTA and GFF files.
         */
        collect_best_annotation_files(
            qc_metrics_after_filtering,
            best_annotations
        ).set { best_annotation_files }

        best_annotation_files
            .map { id, meta, annotation_method, fasta, gff ->
                 [ id, meta, fasta ]
            }
            .set { best_proteins }

        best_annotation_files
            .map { id, meta, annotation_method, fasta, gff ->
                [ id, meta, gff ]
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
        
        best_gff
            .join(best_proteins, by: [0, 1])
            .map { id, meta, gff, proteins -> [id, meta, gff, proteins] }
            .set { best_gff_and_proteins }

        ANNOTATE_PROTEIN_FASTA(best_gff_and_proteins)

    // emit:
    //     best_annotations = best_annotations
    //     best_gff         = best_gff
    //     best_proteins    = best_proteins
}