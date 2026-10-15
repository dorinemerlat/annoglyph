include { UNMASK_GENOME                                 } from "$projectDir/modules/local/repeats_annotation/unmask_genome.nf"
include { REPEATMODELER_BUILDDATABASE } from "$projectDir/modules/nf-core/repeatmodeler/builddatabase/main"
include { REPEATMODELER_REPEATMODELER } from "$projectDir/modules/nf-core/repeatmodeler/repeatmodeler/main"
include { CLEAN_REPEAT_FAMILIES                         } from "$projectDir/modules/local/repeats_annotation/clean_repeat_families.nf"
include { HMMSCAN_RRNA                                  } from "$projectDir/modules/local/repeats_annotation/hmmscan_rrna.nf"
include { HMMSCAN_BUSCO_GENES                           } from "$projectDir/modules/local/repeats_annotation/hmmscan_busco_genes.nf"
include { MCHELPER                                      } from "$projectDir/modules/local/repeats_annotation/mchelper.nf"
include { POSTPROCESS_MCHELPER                          } from "$projectDir/modules/local/repeats_annotation/postprocess_mchelper.nf"
include { MERGE_TABLES as MERGE_POSTPROCESS_MCHELPER    } from "$projectDir/modules/local/shared_modules/merge_tables.nf"
include { CDHIT_CDHITEST } from "$projectDir/modules/nf-core/cdhit/cdhitest/main"
include { MERGE_FILES as MERGE_REPEAT_LIBRARIES } from "$projectDir/modules/local/shared_modules/merge_files.nf"
include { REPEATMASKER_REPEATMASKER } from "$projectDir/modules/nf-core/repeatmasker/repeatmasker/main"
include { REPEATMASKER_OUT_TO_GFF                       } from "$projectDir/modules/local/repeats_annotation/repeatmasker_out_to_gff.nf"
// include { PLOT_REPEATS } from "$projectDir/modules/local/repeats_annotation/plot_repeats.nf"
include { SUMMARIZE_REPEATS_BY_LEVEL                    } from "$projectDir/modules/local/repeats_annotation/summarize_repeats_by_level.nf"
include { MERGE_TABLES as MERGE_SUMMARIES               } from "$projectDir/modules/local/shared_modules/merge_tables.nf"
include { REPEATLANDSCAPE                               } from "$projectDir/modules/local/repeats_annotation/repeatlandscape.nf"
include { REFORMAT_REPEATLANDSCAPE                      } from "$projectDir/modules/local/repeats_annotation/reformat_repeatlandscape.nf"
include { MERGE_TABLES as MERGE_REPEATLANDSCAPE         } from "$projectDir/modules/local/shared_modules/merge_tables.nf"

include { PLOT_REPEAT_LANDSCAPE } from "$projectDir/modules/local/repeats_annotation/plot_landscape.nf"

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
        REPEATMODELER_BUILDDATABASE(genomes)
        REPEATMODELER_REPEATMODELER(REPEATMODELER_BUILDDATABASE.out.db)

        // Standardize repeat family identifiers and classifications
        repeat_classification_table = file(params.repeat_classification_table, checkIfExists: true)
        CLEAN_REPEAT_FAMILIES(REPEATMODELER_REPEATMODELER.out.fasta, repeat_classification_table)

        // Run hmmscan to identify rRNA and BUSCO genes in the genomes
        // Associate each repeat library with its BUSCO dataset.
        CLEAN_REPEAT_FAMILIES.out.for_mchelper
            .map { meta, families -> [meta.busco_dataset, meta, families] }
            .combine(busco_datasets.map { meta, directory -> [meta.id, directory] }, by: 0)
            .map { busco_dataset, meta, families, busco_dataset_dir -> [meta, families, busco_dataset, busco_dataset_dir] }
            .set { repeats_families_and_busco_datasets }

        HMMSCAN_BUSCO_GENES(repeats_families_and_busco_datasets)
        HMMSCAN_RRNA(CLEAN_REPEAT_FAMILIES.out.for_mchelper)

        // Combine genomes, repeat libraries and BUSCO datasets to create inputs for MCHelper.
        repeats_families_and_busco_datasets
            .join(genomes, by: 0)
            .map { meta, families, busco_dataset, busco_dataset_dir, genome -> [meta, genome, families, busco_dataset, busco_dataset_dir] }
            .set { mchelper_inputs }

        MCHELPER(mchelper_inputs)

        // Merge RepeatModeler, MCHelper and HMMER results into a single repeat annotation table.
        CLEAN_REPEAT_FAMILIES.out.table
            .join(CLEAN_REPEAT_FAMILIES.out.for_mchelper, by: 0)
            .join(MCHELPER.out.main_results, by: 0)
            .join(HMMSCAN_BUSCO_GENES.out, by: 0)
            .join(HMMSCAN_RRNA.out, by: 0)
            .map { meta, repeats_table, families, mchelper_classif, classified_module, unclassified_module, gene_tbl, rrna_tbl ->
                [meta, repeats_table, families, mchelper_classif, classified_module, unclassified_module, gene_tbl, rrna_tbl, repeat_classification_table ]
            }
            .set { postprocess_mchelper_inputs }

        POSTPROCESS_MCHELPER( postprocess_mchelper_inputs )
        
        // Collect all processed repeat family tables into a single project-wide summary.
        POSTPROCESS_MCHELPER.out.table
            .map { meta, table -> ['all', table] }
            .groupTuple(by: [0])
            .map { id, tables -> [[id: id], tables, "all_families.tsv"] }
            .set { tables_for_collecting }

        MERGE_POSTPROCESS_MCHELPER( tables_for_collecting )

        // Remove redundant consensus sequences before RepeatMasker.
        // Consensus sequences can either be clustered independently for each species
        // or pooled across all species, depending on group_consensus_sequences.
        if (params.group_consensus_sequences) {
            // Build one shared repeat library for all species.
            all_repeat_libraries = POSTPROCESS_MCHELPER.out.families
                .map { meta, families -> ['repeats_families', families] }
                .groupTuple(by: 0)
                .map { id, families_list ->
                    def libraries = (families_list + [params.reference_library]).findAll { it != null }
                    def inputs = libraries.size() == 1 ? libraries.first() : libraries

                    [[id: id], inputs, "${id}_merged_input.fa"]
                }
        } else {
            // Build one repeat library per species.
            all_repeat_libraries = POSTPROCESS_MCHELPER.out.families
                .map { meta, families ->
                    def libraries = [families, params.reference_library].findAll { it != null }
                    def inputs = libraries.size() == 1 ? libraries.first() : libraries

                    [meta, inputs, "${meta.id}_merged_input.fa"]
                }
        }

        MERGE_REPEAT_LIBRARIES(all_repeat_libraries)
        CDHIT_CDHITEST(MERGE_REPEAT_LIBRARIES.out)

        // Associate each genome with the appropriate repeat library.
        if (params.group_consensus_sequences) {
            def grouped_library = CDHIT_CDHITEST.out.fasta.map { meta, library -> library }

            genome_with_libraries = genomes.combine(grouped_library)
        } else {
            genome_with_libraries = genomes.join(CDHIT_CDHITEST.out.fasta)
        }

        REPEATMASKER_REPEATMASKER(genome_with_libraries)

        // Convert RepeatMasker annotations into a standardized GFF3 using the processed repeat family classifications.
        REPEATMASKER_REPEATMASKER.out.out
            .combine( MERGE_POSTPROCESS_MCHELPER.out )
            .map { meta, out, collect_meta, collect_table -> [meta, out, collect_table, repeat_classification_table] }
            .set { repeatmasker_outputs_for_gff }

        REPEATMASKER_OUT_TO_GFF( repeatmasker_outputs_for_gff )

        REPEATMASKER_OUT_TO_GFF.out
            .map { meta, gff -> ['all', gff] }
            .groupTuple()
            .set { gff_files }

        // Compute repeat landscapes from RepeatMasker alignments.
        REPEATLANDSCAPE(REPEATMASKER_REPEATMASKER.out.align)

        REPEATLANDSCAPE.out.tsv.map { meta, tsv -> [meta, tsv, file("$projectDir/data/repeats_classification.csv")] }
            .set { reformat_repeatlandscape_inputs }

        // Replace RepeatMasker classifications with the standardized EXOGAP classification.
        REFORMAT_REPEATLANDSCAPE(reformat_repeatlandscape_inputs)

        PLOT_REPEAT_LANDSCAPE(REFORMAT_REPEATLANDSCAPE.out)

        // Merge repeat landscape tables from all genomes.
        REFORMAT_REPEATLANDSCAPE.out
            .map { meta, table -> ['all', table] }
            .groupTuple(by: [0])
            .map { id, tables -> [[id: id], tables, "all_repeatlandscape.tsv"] }
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

    //     // genome_stats.map { meta, stats -> ['all', stats] }
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

    emit:
    versions = REPEATMODELER_BUILDDATABASE.out.versions.mix(REPEATMODELER_REPEATMODELER.out.versions).mix(CDHIT_CDHITEST.out.versions).mix(REPEATMASKER_REPEATMASKER.out.versions)
        consensus_libraries = genome_with_libraries.map { meta, genome, library -> [meta, library] }
        families_table = POSTPROCESS_MCHELPER.out.table
        landscapes = REFORMAT_REPEATLANDSCAPE.out
        landscape_plots = PLOT_REPEAT_LANDSCAPE.out.plots
        comparative_tables = MERGE_POSTPROCESS_MCHELPER.out.mix(MERGE_REPEATLANDSCAPE.out)
        masked_genomes = REPEATMASKER_REPEATMASKER.out.masked
        repeat_gff = REPEATMASKER_OUT_TO_GFF.out
}

