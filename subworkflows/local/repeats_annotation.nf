include { UNMASK_GENOME                                 } from "$projectDir/modules/repeats_annotation/unmask_genome.nf"
include { REPEATMODELER                                 } from "$projectDir/modules/repeats_annotation/repeatmodeler.nf"
include { CLEAN_REPEAT_FAMILIES                         } from "$projectDir/modules/repeats_annotation/clean_repeat_families.nf"
include { HMMSCAN_RRNA                                  } from "$projectDir/modules/repeats_annotation/hmmscan_rrna.nf"
include { HMMSCAN_BUSCO_GENES                           } from "$projectDir/modules/repeats_annotation/hmmscan_busco_genes.nf"
include { MCHELPER                                      } from "$projectDir/modules/repeats_annotation/mchelper.nf"
include { POSTPROCESS_MCHELPER                          } from "$projectDir/modules/repeats_annotation/postprocess_mchelper.nf"
include { MERGE_TABLES as MERGE_POSTPROCESS_MCHELPER    } from "$projectDir/modules/shared_modules/merge_tables.nf"
include { CD_HIT_EST                                    } from "$projectDir/modules/repeats_annotation/cd_hit_est.nf"
include { REPEATMASKER                                  } from "$projectDir/modules/repeats_annotation/repeatmasker.nf"
include { REPEATMASKER_OUT_TO_GFF                       } from "$projectDir/modules/repeats_annotation/repeatmasker_out_to_gff.nf"
// include { PLOT_REPEATS } from "$projectDir/modules/repeats_annotation/plot_repeats.nf"
include { SUMMARIZE_REPEATS_BY_LEVEL                    } from "$projectDir/modules/repeats_annotation/summarize_repeats_by_level.nf"
include { MERGE_TABLES as MERGE_SUMMARIES               } from "$projectDir/modules/shared_modules/merge_tables.nf"
include { REPEATLANDSCAPE                               } from "$projectDir/modules/repeats_annotation/repeatlandscape.nf"
include { REFORMAT_REPEATLANDSCAPE                      } from "$projectDir/modules/repeats_annotation/reformat_repeatlandscape.nf"
include { MERGE_TABLES as MERGE_REPEATLANDSCAPE         } from "$projectDir/modules/shared_modules/merge_tables.nf"

workflow REPEATS_ANNOTATION {
    take:
        genomes
        newick
        genome_stats
        busco_datasets

    main:
        // Some genomes have softmasked regions, which can correspond to low complexity or low coverage regions.
        // These softmasked regions can cause problems for repeat annotation, so we will unmask them before running RepeatModeler.
        UNMASK_GENOME(genomes)
        genomes = UNMASK_GENOME.out

        // Run RepeatModeler to identify de novo repeats
        REPEATMODELER(genomes)

        // Standardize repeat family identifiers and classifications
        repeat_classification_table = file(params.repeat_classification_table, checkIfExists: true)
        CLEAN_REPEAT_FAMILIES(REPEATMODELER.out.families, repeat_classification_table)

        // Run hmmscan to identify rRNA and BUSCO genes in the genomes
        CLEAN_REPEAT_FAMILIES.out.for_mchelper
            .map {id, meta, families -> [meta.busco_dataset, id, meta, families] }
            .combine( busco_datasets, by: 0 )
            .map { busco_dataset, id, meta, families, busco_dataset_dir -> [id, meta, families, busco_dataset, busco_dataset_dir] }
            .set { repeats_families_and_busco_datasets }

        HMMSCAN_BUSCO_GENES(repeats_families_and_busco_datasets)
        HMMSCAN_RRNA(CLEAN_REPEAT_FAMILIES.out.for_mchelper)

        // Run MCHelper
        repeats_families_and_busco_datasets
            .join(genomes, by: [0, 1])
            .map {id, meta, families, busco_dataset, busco_dataset_dir, genome -> [id, meta, genome, families, busco_dataset, busco_dataset_dir] }
            .set { mchelper_inputs }

        MCHELPER(mchelper_inputs)

        CLEAN_REPEAT_FAMILIES.out.table
            .join(CLEAN_REPEAT_FAMILIES.out.for_mchelper, by: [0, 1])
            .join(MCHELPER.out.main_results, by: [0, 1])
            .join(HMMSCAN_BUSCO_GENES.out, by: [0, 1])
            .join(HMMSCAN_RRNA.out, by: [0, 1])
            .map {id, meta, repeats_table, families, mchelper_classif, classified_module, unclassified_module, gene_tbl, rrna_tbl ->
                [ id, meta, repeats_table, families, mchelper_classif, classified_module, unclassified_module, gene_tbl, rrna_tbl, repeat_classification_table ]
            }
            .set { postprocess_mchelper_inputs }

        POSTPROCESS_MCHELPER( postprocess_mchelper_inputs )

        POSTPROCESS_MCHELPER.out.table
            .map { id, meta, table -> ['all', table] }
            .groupTuple(by: [0])
            .map { id, tables -> [id, tables, "all_families.tsv"] }
            .set { tables_for_collecting }

        MERGE_POSTPROCESS_MCHELPER( tables_for_collecting )

        // CD-HIT-EST to remove redundancy
        if (params.group_consensus_sequences) {
            all_repeat_libraries = POSTPROCESS_MCHELPER.out.families
                .map { id, meta, families -> ['repeats_families', families] }
                .groupTuple(by: 0)
                .map { id, families_list ->
                    def libraries = (families_list + [params.reference_library]).findAll { it != null }
                    def inputs = libraries.size() == 1 ? libraries.first() : libraries

                    [id, inputs, 0.8, 80, 0.8]
                }
        } else {
            all_repeat_libraries = POSTPROCESS_MCHELPER.out.families
                .map { id, meta, families ->
                    def libraries = [families, params.reference_library].findAll { it != null }
                    def inputs = libraries.size() == 1 ? libraries.first() : libraries

                    [id, inputs, 0.8, 80, 0.8]
                }
        }

        CD_HIT_EST(all_repeat_libraries)

        if (params.group_consensus_sequences) {
            def grouped_library = CD_HIT_EST.out.library.map { id, library -> library }

            genome_with_libraries = genomes.combine(grouped_library)
        } else {
            genome_with_libraries = genomes.join(CD_HIT_EST.out.library)
        }

        REPEATMASKER(genome_with_libraries)

        REPEATMASKER.out.out
            .combine( MERGE_POSTPROCESS_MCHELPER.out )
            .map { id, meta, out, collect_id, collect_table -> [id, meta, out, collect_table, repeat_classification_table] }
            .set { repeatmasker_outputs_for_gff }

        REPEATMASKER_OUT_TO_GFF( repeatmasker_outputs_for_gff )

        REPEATMASKER_OUT_TO_GFF.out
            .map { id, meta, gff -> ['all', gff] }
            .groupTuple()
            .set { gff_files }

        REPEATLANDSCAPE(REPEATMASKER.out.align)

        REPEATLANDSCAPE.out.tsv.map { id, meta, tsv -> [id, meta, tsv, file("$projectDir/data/repeats_classification.csv")] }
            .set { reformat_repeatlandscape_inputs }

        REFORMAT_REPEATLANDSCAPE(reformat_repeatlandscape_inputs)

        REFORMAT_REPEATLANDSCAPE.out
            .map { id, meta, table -> ['all', table] }
            .groupTuple(by: [0])
            .map { id, tables -> [id, tables, "all_repeatlandscape.tsv"] }
            .set { tables_for_merge_repeatlandscape }

        MERGE_REPEATLANDSCAPE( tables_for_merge_repeatlandscape )

    //     SUMMARIZE_REPEATS_BY_LEVEL( REPEATMASKER_OUT_TO_GFF.out )

    //     // gather summaries
    //     summaries_by_all            = RepeatAnnotationUtils.gatherSummaryTables(SUMMARIZE_REPEATS_BY_LEVEL.out.all, 'repeats_summary_by_all')
    //     summaries_by_all_filtered   = RepeatAnnotationUtils.gatherSummaryTables(SUMMARIZE_REPEATS_BY_LEVEL.out.all, 'repeats_summary_all_only_repetitive_elements')
    //     summaries_by_feature        = RepeatAnnotationUtils.gatherSummaryTables(SUMMARIZE_REPEATS_BY_LEVEL.out.by_feature, 'repeats_summary_by_feature')
    //     summaries_by_type           = RepeatAnnotationUtils.gatherSummaryTables(SUMMARIZE_REPEATS_BY_LEVEL.out.by_type, 'repeats_summary_by_type')
    //     summaries_by_class          = RepeatAnnotationUtils.gatherSummaryTables(SUMMARIZE_REPEATS_BY_LEVEL.out.by_class, 'repeats_summary_by_class')
    //     summaries_by_order          = RepeatAnnotationUtils.gatherSummaryTables(SUMMARIZE_REPEATS_BY_LEVEL.out.by_order, 'repeats_summary_by_order')
    //     summaries_by_superfamily    = RepeatAnnotationUtils.gatherSummaryTables(SUMMARIZE_REPEATS_BY_LEVEL.out.by_superfamily, 'repeats_summary_by_superfamily')
    //     summaries_by_family         = RepeatAnnotationUtils.gatherSummaryTables(SUMMARIZE_REPEATS_BY_LEVEL.out.by_family, 'repeats_summary_by_family')
    //     summaries_by_subfamily      = RepeatAnnotationUtils.gatherSummaryTables(SUMMARIZE_REPEATS_BY_LEVEL.out.by_subfamily, 'repeats_summary_by_subfamily')
    //     summaries_by_subfamily2     = RepeatAnnotationUtils.gatherSummaryTables(SUMMARIZE_REPEATS_BY_LEVEL.out.by_subfamily2, 'repeats_summary_by_subfamily2')


    //     summaries_by_all
    //         .concat( summaries_by_all_filtered )
    //         .concat( summaries_by_feature )
    //         .concat( summaries_by_type )
    //         .concat( summaries_by_class )
    //         .concat( summaries_by_order )
    //         .concat( summaries_by_superfamily )
    //         .concat( summaries_by_family )
    //         .concat( summaries_by_subfamily )
    //         .concat( summaries_by_subfamily2 )
    //         .set { all_summaries }

    //     // all_summaries.map{ id, files, output -> [id, files.size()]}.view()
    //     MERGE_SUMMARIES( all_summaries )

    //     // genome_stats.map { id, meta, stats -> ['all', stats] }
    //     //     .groupTuple()
    //     //     .set { stats_files }

    //     // newick
    //     //     .map {file -> ['all', file]}
    //     //     .set { newick }

    //     // gff_files
    //     //     .join( stats_files )
    //     //     .join( newick )
    //     //     .set { inputs_for_plot_repeats }

    //     // PLOT_REPEATS(inputs_for_plot_repeats)

    // emit:
    //     masked_genomes = REPEATMASKER.out.masked
    //     repeat_gff = REPEATMASKER_OUT_TO_GFF.out
}

