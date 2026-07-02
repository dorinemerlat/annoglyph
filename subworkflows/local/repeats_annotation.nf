include { UNMASK_GENOME } from "$projectDir/modules/local/module_repeats/unmask_genome.nf"
include { REPEATMODELER } from "$projectDir/modules/local/module_repeats/repeatmodeler.nf"
include { CLEAN_REPEAT_FAMILIES } from "$projectDir/modules/local/module_repeats/clean_repeat_families.nf"
include { DOWNLOAD_RRNA_EUKARYOTA } from "$projectDir/modules/local/module_repeats/download_rrna_eukaryota.nf"
include { DOWNLOAD_BUSCO_HMM } from "$projectDir/modules/local/module_repeats/download_busco_hmm.nf"
include { MCHELPER } from "$projectDir/modules/local/module_repeats/mchelper.nf"
include { HMMSCAN } from "$projectDir/modules/local/module_repeats/hmmscan.nf"
include { POSTPROCESS_MCHELPER } from "$projectDir/modules/local/module_repeats/postprocess_mchelper.nf"
include { MERGE_TABLES } from "$projectDir/modules/local/merge_tables.nf"
include { CD_HIT } from "$projectDir/modules/local/module_repeats/cd_hit.nf"
include { REPEATMASKER } from "$projectDir/modules/local/module_repeats/repeatmasker.nf"
include { REPEATMASKER_OUT_TO_GFF } from "$projectDir/modules/local/module_repeats/repeatmasker_out_to_gff.nf"
// include { PLOT_REPEATS } from "$projectDir/modules/local/module_repeats/plot_repeats.nf"
include { SUMMARIZE_REPEATS_BY_LEVEL } from "$projectDir/modules/local/module_repeats/summarize_repeats_by_level.nf"
include { MERGE_TABLES as MERGE_SUMMARIES} from "$projectDir/modules/local/merge_tables.nf"
include { REPEATLANDSCAPE } from "$projectDir/modules/local/module_repeats/repeatlandscape.nf"
include { REFORMAT_REPEATLANDSCAPE } from "$projectDir/modules/local/module_repeats/reformat_repeatlandscape.nf"
include { MERGE_TABLES as MERGE_REPEATLANDSCAPE} from "$projectDir/modules/local/merge_tables.nf"

def gather_summary_ch(ch, prefix) {
    ch
      .map { id, meta, file -> [prefix, file] }
      .groupTuple()
      .map { key, files -> [key, files, prefix + '.tsv'] }
}

workflow REPEATS_ANNOTATION {
    take:
        genomes
        newick
        genome_stats

    main:
        println "Run repeats annotation"

        // Some genomes have softmasked regions, which can correspond to low complexity or low coverage regions.
        // These softmasked regions can cause problems for repeat annotation, so we will unmask them before running RepeatModeler.
        UNMASK_GENOME(genomes)
        UNMASK_GENOME.out.set { genomes }

        // Run RepeatModeler to identify de novo repeats
        REPEATMODELER( genomes )

        // Clean up repeat libraries name
        CLEAN_REPEAT_FAMILIES( REPEATMODELER.out, file("$projectDir/data/repeats_classification.csv") )

        // Download rRNA HMM database and prepare it for hmmscan
        DOWNLOAD_RRNA_EUKARYOTA()
        DOWNLOAD_BUSCO_HMM(genomes.map {id, meta, genome -> [meta.busco_dataset, id] }.groupTuple(by: 0) )

        CLEAN_REPEAT_FAMILIES.out.for_mchelper
            .combine( DOWNLOAD_RRNA_EUKARYOTA.out.map { it -> [it]} )
            .map { id, meta, families, rRNA_hmm -> [id, 'rRNA', families, rRNA_hmm] }
            .set { hmm_rRNA_inputs }
            
        DOWNLOAD_BUSCO_HMM.out
            .flatMap { dataset_name, species_list, dataset_hmms -> species_list.collect { id -> [id, dataset_hmms] }}
            .set { busco_hmms }

        CLEAN_REPEAT_FAMILIES.out.for_mchelper
            .join(busco_hmms, by: 0)
            .map { id, meta, families, busco_hmm -> [id, 'genes', families, busco_hmm] }
            .set { hmm_genes_inputs }
 
        hmm_rRNA_inputs.concat( hmm_genes_inputs )
            .set { all_hmmscan_inputs }

        HMMSCAN( all_hmmscan_inputs )

        // Run MCHelper
        CLEAN_REPEAT_FAMILIES.out.for_mchelper
            .join(genomes, by: [0, 1])
            .map {id, meta, families, genome -> [id, meta, families, genome] }
            .set { all_families }

        MCHELPER( all_families )

        // Collect all results for process MChelper
        HMMSCAN.out.branch {id, molecule, tblout ->
            gene: molecule == 'genes'
            rrna: molecule == 'rRNA'
            }
            .set { hmmscan_results }

        CLEAN_REPEAT_FAMILIES.out.table
            .join(CLEAN_REPEAT_FAMILIES.out.for_mchelper, by: [0,1])
            .join(MCHELPER.out.main_results, by: [0,1])
            .join(hmmscan_results.gene, by: [0])
            .join(hmmscan_results.rrna, by: [0])
            .map {id, meta, table, rm2_families, mchelper_classif, classified_module, unclassified_module, gene_name, gene_tbl, rrna_name, rrna_tbl ->
                [id, meta, table, rm2_families, mchelper_classif, classified_module, unclassified_module, gene_tbl, rrna_tbl, file("$projectDir/data/repeats_classification.csv")]
            }
            .set { postprocess_mchelper_inputs }

        POSTPROCESS_MCHELPER( postprocess_mchelper_inputs )

        POSTPROCESS_MCHELPER.out.table
            .map { id, meta, table -> ['all', table] }
            .groupTuple(by: [0])
            .map { id, tables -> [id, tables, "all_families.tsv"] }
            .set { tables_for_collecting }

        MERGE_TABLES( tables_for_collecting )

        // CD-HIT to remove redundancy
        POSTPROCESS_MCHELPER.out.families
            .map { id, meta, families -> ['repeats_families', families] }
            .groupTuple(by: [0])
            .map{id, families_list -> [id, 'meta', families_list + params.reference_library, "nucleotide", "0.8", "80", "0.8"]}
            .set { all_repeat_libraries }

        CD_HIT(all_repeat_libraries)

    //     // // ANALYZE_FAMILIES_CLUSTER(CD_HIT.out.clusters.map { id, meta, cluster -> [id, meta, cluster, 'repeat_cluster' ] } )

        genomes
            .combine( CD_HIT.out.library.map { id, meta, library -> library } )
            .set { genome_with_libraries }

        REPEATMASKER(genome_with_libraries)

        REPEATMASKER.out.out
            .combine( MERGE_TABLES.out )
            .map { id, meta, out, collect_id, collect_table -> [id, meta, out, collect_table, file("$projectDir/data/repeats_classification.csv")] }
            .set { repeatmasker_outputs_for_gff }

        REPEATMASKER_OUT_TO_GFF( repeatmasker_outputs_for_gff )

        REPEATMASKER_OUT_TO_GFF.out
            .map { id, meta, gff -> ['all', gff] }
            .groupTuple()
            .set { gff_files }

        SUMMARIZE_REPEATS_BY_LEVEL( REPEATMASKER_OUT_TO_GFF.out )

        // gather summaries
        summaries_by_all            = gather_summary_ch(SUMMARIZE_REPEATS_BY_LEVEL.out.all, 'repeats_summary_by_all')
        summaries_by_all_filtered   = gather_summary_ch(SUMMARIZE_REPEATS_BY_LEVEL.out.all, 'repeats_summary_all_only_repetitive_elements')
        summaries_by_feature        = gather_summary_ch(SUMMARIZE_REPEATS_BY_LEVEL.out.by_feature, 'repeats_summary_by_feature')
        summaries_by_type           = gather_summary_ch(SUMMARIZE_REPEATS_BY_LEVEL.out.by_type, 'repeats_summary_by_type')
        summaries_by_class          = gather_summary_ch(SUMMARIZE_REPEATS_BY_LEVEL.out.by_class, 'repeats_summary_by_class')
        summaries_by_order          = gather_summary_ch(SUMMARIZE_REPEATS_BY_LEVEL.out.by_order, 'repeats_summary_by_order')
        summaries_by_superfamily    = gather_summary_ch(SUMMARIZE_REPEATS_BY_LEVEL.out.by_superfamily, 'repeats_summary_by_superfamily')
        summaries_by_family         = gather_summary_ch(SUMMARIZE_REPEATS_BY_LEVEL.out.by_family, 'repeats_summary_by_family')
        summaries_by_subfamily      = gather_summary_ch(SUMMARIZE_REPEATS_BY_LEVEL.out.by_subfamily, 'repeats_summary_by_subfamily')
        summaries_by_subfamily2     = gather_summary_ch(SUMMARIZE_REPEATS_BY_LEVEL.out.by_subfamily2, 'repeats_summary_by_subfamily2')


        summaries_by_all
            .concat( summaries_by_all_filtered )
            .concat( summaries_by_feature )
            .concat( summaries_by_type )
            .concat( summaries_by_class )
            .concat( summaries_by_order )
            .concat( summaries_by_superfamily )
            .concat( summaries_by_family )
            .concat( summaries_by_subfamily )
            .concat( summaries_by_subfamily2 )
            .set { all_summaries }

        // all_summaries.map{ id, files, output -> [id, files.size()]}.view()
        MERGE_SUMMARIES( all_summaries )

        // genome_stats.map { id, meta, stats -> ['all', stats] }
        //     .groupTuple()
        //     .set { stats_files }

        // newick
        //     .map {file -> ['all', file]}
        //     .set { newick }

        // gff_files
        //     .join( stats_files )
        //     .join( newick )
        //     .set { inputs_for_plot_repeats }

        // PLOT_REPEATS(inputs_for_plot_repeats)

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

    emit:
        unmasked_genomes = UNMASK_GENOME.out
        masked_genomes = REPEATMASKER.out.masked
        repeat_gff = REPEATMASKER_OUT_TO_GFF.out
}