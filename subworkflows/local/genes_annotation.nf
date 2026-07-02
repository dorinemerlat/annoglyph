// include modules
include { BRAKER2                                   } from "$projectDir/modules/local/module_genes/braker2.nf"
include { BRAKER3                                   } from "$projectDir/modules/local/module_genes/braker3.nf"
include { GALBA                                     } from "$projectDir/modules/local/module_genes/galba.nf"
include { REFORMAT_GFF3                             } from "$projectDir/modules/local/module_genes/reformat_gff3.nf"
include { AGAT_STATISTICS                           } from "$projectDir/modules/local/module_genes/agat_statistics.nf"
include { AGAT_FILTER_GENE_BY_LENGTH                } from "$projectDir/modules/local/module_genes/agat_filter_gene_by_length.nf"
include { EXTRACT_CANONICAL_PROTEINS                } from "$projectDir/modules/local/module_genes/extract_canonical_proteins.nf"
include { DOWNLOAD_OMAMER_DB                        } from "$projectDir/modules/local/module_genes/download_omamer_db.nf"
include { OMARK                                     } from "$projectDir/modules/local/module_genes/omark.nf"
include { BUSCO                                     } from "$projectDir/modules/local/module_genes/busco.nf"
include { FILTER_OMARK_CONTAMINATIONS               } from "$projectDir/modules/local/module_genes/filter_omark_contaminations.nf"
include { FILTER_OMARK_CONTAMINATIONS_IN_GFF        } from "$projectDir/modules/local/module_genes/filter_omark_contaminations_in_gff.nf"
include { OMARK as OMARK_ON_CONTAMINATION_FILTERED  } from "$projectDir/modules/local/module_genes/omark.nf"
include { BUSCO as BUSCO_ON_CONTAMINATION_FILTERED  } from "$projectDir/modules/local/module_genes/busco.nf"
include { CLASSIFY_GENES_VS_REPEATS                 } from "$projectDir/modules/local/module_genes/classify_genes_vs_repeats.nf"
include { BLASTP                                    } from "$projectDir/modules/local/module_genes/blastp.nf"
include { BUILD_BLASTP_TAXONOMY                     } from "$projectDir/modules/local/module_genes/build_blastp_taxonomy.nf"
include { ADD_TAXONOMY_TO_BLASTP                    } from "$projectDir/modules/local/module_genes/add_taxonomy_to_blastp.nf"
include { INTERPROSCAN                              } from "$projectDir/modules/local/module_genes/interproscan.nf"
include { EGGNOG_MAPPER                             } from "$projectDir/modules/local/module_genes/eggnog_mapper.nf"
include { TRANSFER_FUNCTIONAL_ANNOTATIONS           } from "$projectDir/modules/local/module_genes/transfer_functional_annotations.nf"
include { AGAT_MANAGE_FUNCTIONAL_ANNOTATION         } from "$projectDir/modules/local/module_genes/agat_manage_functional_annotation.nf"
include { DOWNLOAD_GO_DB                            } from "$projectDir/modules/local/module_genes/download_go_db.nf"
include { GENERATE_ANNOTATION_TABLES                } from "$projectDir/modules/local/module_genes/generate_annotation_tables.nf"
include { MERGE_TABLES                              } from "$projectDir/modules/local/merge_tables.nf"
include { MERGE_TABLES as MERGE_BEST_TABLES         } from "$projectDir/modules/local/merge_tables.nf"
include { SELECT_BEST_ANNOTATION                    } from "$projectDir/modules/local/module_genes/select_best_annotation.nf"
include { PLOT_ANNOTATION_COMPARISON                } from "$projectDir/modules/local/module_genes/plot_annotation_comparison.nf"
include { PLOT_ANNOTATION_ANALYSIS                  } from "$projectDir/modules/local/module_genes/plot_annotation_analysis.nf"
include { PUBLISH_BEST_ANNOTATIONS                  } from "$projectDir/modules/local/module_genes/publish_best_annotations.nf"
def branch_by_filtering_status(channel) {
    channel.branch {
        before_contamination_filtering: it[3] == "before-contamination-filtering"
        after_contamination_filtering:  it[3] == "after-contamination-filtering"
    }
}

def prepare_for_comparison(channel, label) {
    channel.map { id, meta, annotation_method, filtering_status, file -> [filtering_status + '-' + label, file]}
        .groupTuple(by: 0)
}

def collect_and_merge_without_contamination_tables(before, after, label) {
    before
        .concat (after)
        .groupTuple(by: [0, 1, 2])
        .map { id, meta, annotation_method, filtering_status, file -> [ 'after-contamination-filtering-' + label, file[-1]] }
        .groupTuple(by: 0)
}

def collect_without_contamination_tables(before, after, label) {
    before
        .concat (after)
        .groupTuple(by: [0, 1, 2])
        .map { id, meta, annotation_method, filtering_status, file -> [ id, annotation_method, 'after-contamination-filtering-' + label, file[-1]] }
}

def collect_best_tables(best, before, after, label) {
    best
        .join(collect_without_contamination_tables(before, after, label), by: [0, 1])
        .map { id, annotation_method, status, file -> [label, file] }
        .groupTuple(by: 0)
} 


def collect_best_annotation_tables(best_annotations, summaries, protein_lengths, intron_lengths, label) {
    def best_summaries = collect_best_tables(
        best_annotations,
        summaries.before_contamination_filtering,
        summaries.after_contamination_filtering,
        label + '-summary'
    )

    def best_protein_lengths = collect_best_tables(
        best_annotations,
        protein_lengths.before_contamination_filtering,
        protein_lengths.after_contamination_filtering,
        label + '-proteins'
    )

    def best_intron_lengths = collect_best_tables(
        best_annotations,
        intron_lengths.before_contamination_filtering,
        intron_lengths.after_contamination_filtering,
        label + '-introns'
    )

    return best_summaries
        .concat(best_protein_lengths, best_intron_lengths)
        .map { status, file ->
            tuple(status, file, "${status}.tsv")
        }
}

def collect_best_file(best_annotations, before_files, after_files) {
    def before_clean = before_files
        .map { id, meta, annotation_method, file ->
            tuple(id, annotation_method, meta, "before", file)
        }

    def after_01 = after_files
        .filter { id, meta, annotation_method, threshold, file -> threshold == "0.1" }
        .map { id, meta, annotation_method, threshold, file ->
            tuple(id, annotation_method, meta, "after", file)
        }

    return before_clean
        .concat(after_01)
        .groupTuple(by: [0, 1])
        .map { id, annotation_method, meta_list, status_list, file_list ->
            def idx = status_list.findIndexOf { it == "after" }
            def selected_idx = idx >= 0 ? idx : 0

            tuple(
                id,
                annotation_method,
                meta_list[selected_idx],
                file_list[selected_idx]
            )
        }
        .join(best_annotations, by: [0, 1])
        .map { id, annotation_method, meta, file ->
            tuple(id, meta, file)
        }
}

workflow GENES_ANNOTATION {
    take:
        masked_genomes
        bams
        repeats_gff
        newick
        all_genome_stats

    main:
        masked_genomes
            .combine(bams, by:[0, 1])
            .map { id, meta, genome, aligner, bam -> [id, meta, aligner, genome, bam, file(params.protein_set)] }
            .set { braker3_input }
        
        BRAKER3(braker3_input)
        
        masked_genomes
            .map { id, meta, genome -> [id, meta, genome, file(params.protein_set)] }
            .set { braker2_and_galba_input }

        BRAKER2(braker2_and_galba_input)
        GALBA(braker2_and_galba_input)

        BRAKER3.out.gff.map { id, meta, rnaseq_aligner, gff -> [id, meta, "braker-with-" + rnaseq_aligner, gff] }
            .concat(BRAKER2.out.gff.map { id, meta, gff -> [id, meta, "braker2", gff] })
            .concat(GALBA.out.gff.map { id, meta, gff -> [id,  meta, "galba", gff] })
            // filter empty gff files
            .filter { id, meta, annotation_method, gff -> gff.size() > 0 }
            .set { all_gff }

        REFORMAT_GFF3(all_gff)

        // filtering annotation
        AGAT_FILTER_GENE_BY_LENGTH(REFORMAT_GFF3.out.map { id, meta, annotation_method, gff -> [id, meta, annotation_method, gff, 90] })
        
        AGAT_FILTER_GENE_BY_LENGTH.out
            .combine(masked_genomes, by: [0, 1])
            .map { id, meta, annotation_method, gff, genome -> [id, meta, annotation_method, genome, gff] }
            .set { gff_and_genome }

        EXTRACT_CANONICAL_PROTEINS(gff_and_genome)

        // Functional annotation
        EXTRACT_CANONICAL_PROTEINS.out.fasta
            .map { id, meta, annotation_method, fasta -> [id, meta, annotation_method, fasta, "uniprot_sprot", "1e-6", "3", "20"] }
            .set { blastp_input }

        BLASTP(blastp_input)

        // BLASTP.out.tsv
        //     .map{id, meta, annotation_method, blastp -> blastp }
        //     .collect()
        //     .set { blastp_tsv_outputs }

        // BUILD_BLASTP_TAXONOMY(blastp_tsv_outputs)
        // ADD_TAXONOMY_TO_BLASTP(BLASTP.out.tsv, BUILD_BLASTP_TAXONOMY.out.taxonomy)
        INTERPROSCAN(EXTRACT_CANONICAL_PROTEINS.out.fasta)
        EGGNOG_MAPPER(EXTRACT_CANONICAL_PROTEINS.out.fasta)

        // Add functional annotation to GFF
        EXTRACT_CANONICAL_PROTEINS.out.gff
            .join(BLASTP.out.tsv, by: [0, 1, 2])
            .join(INTERPROSCAN.out.tsv, by: [0, 1, 2])
            .join(EGGNOG_MAPPER.out.annotations, by: [0, 1, 2])
            .set { transfer_functional_annotations_input }
        TRANSFER_FUNCTIONAL_ANNOTATIONS(transfer_functional_annotations_input)


        // statistics on annotation before filtering incomplete genes
        AGAT_STATISTICS(AGAT_FILTER_GENE_BY_LENGTH.out)

        // Run OMARK and BUSCO
        DOWNLOAD_OMAMER_DB()
        EXTRACT_CANONICAL_PROTEINS.out.fasta
            .combine(DOWNLOAD_OMAMER_DB.out)
            .set { omark_input }

        OMARK(omark_input)
        BUSCO(EXTRACT_CANONICAL_PROTEINS.out.fasta)

        OMARK.out.contaminants_count
            .filter { id, meta, annotation_method, contaminants_count -> contaminants_count.readLines()[0].toInteger() > 0 }
            .join ( EXTRACT_CANONICAL_PROTEINS.out.gff, by: [0, 1, 2] )
            .join ( EXTRACT_CANONICAL_PROTEINS.out.fasta, by: [0, 1, 2] )
            .join ( OMARK.out.ump, by: [0, 1, 2] )
            .map { id, meta, annotation_method, contaminants_count, gff, fasta, ump -> [id, meta, annotation_method, gff, fasta, ump, "0.33" ] }
            .set { contamination_filtering_input }

        contamination_filtering_input
            .concat (contamination_filtering_input.map{id, meta, annotation_method, gff, fasta, ump, threshold -> [id, meta, annotation_method, gff, fasta, ump, "0.1"]})
            .concat (contamination_filtering_input.map{id, meta, annotation_method, gff, fasta, ump, threshold -> [id, meta, annotation_method, gff, fasta, ump, "0"]})
            .concat (contamination_filtering_input.map{id, meta, annotation_method, gff, fasta, ump, threshold -> [id, meta, annotation_method, gff, fasta, ump, "1"]})
            .map {id, meta, annotation_method, gff, fasta, ump, threshold -> [ id, meta, annotation_method, gff, fasta, ump, threshold ] }
            .set{contamination_filtering_input}

       // Check if OMARK filters enough contaminations
        FILTER_OMARK_CONTAMINATIONS(contamination_filtering_input)
        FILTER_OMARK_CONTAMINATIONS.out.contaminations
            .combine ( TRANSFER_FUNCTIONAL_ANNOTATIONS.out, by: [0, 1, 2] )
            .set { filter_omark_contaminations_in_gff_input }

        FILTER_OMARK_CONTAMINATIONS_IN_GFF(filter_omark_contaminations_in_gff_input)

        FILTER_OMARK_CONTAMINATIONS.out.fasta
            .map { id, meta, annotation_method, threshold, fasta -> [id, meta, annotation_method + "_c-" + threshold, fasta] }
            .set { contamination_filtered_input }

        BUSCO_ON_CONTAMINATION_FILTERED(contamination_filtered_input)
        OMARK_ON_CONTAMINATION_FILTERED(contamination_filtered_input.combine(DOWNLOAD_OMAMER_DB.out))
        
        // compare predictions
        EXTRACT_CANONICAL_PROTEINS.out.fasta
            .join(TRANSFER_FUNCTIONAL_ANNOTATIONS.out, by: [0, 1, 2])
            .join(BUSCO.out.json, by: [0, 1, 2])
            .join(OMARK.out.sum, by: [0, 1, 2])
            .map { id, meta, annotation_method, fasta, gff, busco, omark -> [id, meta, annotation_method, 'before-contamination-filtering', fasta, gff, busco, omark] }
            .set { compare_annotations_before_contamination_filtering_input }
        
        FILTER_OMARK_CONTAMINATIONS_IN_GFF.out
            .filter { id, meta, annotation_method, threshold, gff -> threshold == "0.1" }
            .map { id, meta, annotation_method, threshold, gff -> [id, meta, annotation_method, gff] }
            .set { filtering_gff }

        BUSCO_ON_CONTAMINATION_FILTERED.out.json
            .map { id, meta, annotation_method, busco_json ->
                def parts = annotation_method.split("_c-", 2)
                def method = parts[0]
                def threshold = parts[1]
                [id, meta, method, threshold, busco_json]
            }
            .filter { id, meta, annotation_method, threshold, busco_json -> threshold == "0.1" }
            .map { id, meta, annotation_method, threshold, busco_json -> [id, meta, annotation_method, busco_json] }
            .set { busco_on_contamination_filtered_json }
        
        OMARK_ON_CONTAMINATION_FILTERED.out.sum
            .map { id, meta, annotation_method, sum ->
                def parts = annotation_method.split("_c-", 2)
                def method = parts[0]
                def threshold = parts[1]
                [id, meta, method, threshold, sum]
            }
            .filter { id, meta, annotation_method, threshold, sum -> threshold == "0.1" }
            .map { id, meta, annotation_method, threshold, sum -> [id, meta, annotation_method, sum] }
            .set { omark_on_contamination_filtered_sum }

        FILTER_OMARK_CONTAMINATIONS.out.fasta.filter { id, meta, annotation_method, threshold, fasta -> threshold == "0.1" }
            .map { id, meta, annotation_method, threshold, fasta -> [id, meta, annotation_method, fasta] }
            .join(filtering_gff, by: [0, 1, 2])
            .join(busco_on_contamination_filtered_json, by: [0, 1, 2])
            .join(omark_on_contamination_filtered_sum, by: [0, 1, 2])
            .map { id, meta, annotation_method, fasta, gff, busco_json, omark_sum -> [id, meta, annotation_method, 'after-contamination-filtering', fasta, gff, busco_json, omark_sum] }
            .set { compare_annotations_after_contamination_filtering_input }

        compare_annotations_before_contamination_filtering_input
            .concat(compare_annotations_after_contamination_filtering_input)
            .set { generate_annotation_tables_input }

        GENERATE_ANNOTATION_TABLES(generate_annotation_tables_input)
        
        branch_by_filtering_status(GENERATE_ANNOTATION_TABLES.out.summary)
            .set { summaries }

        branch_by_filtering_status(GENERATE_ANNOTATION_TABLES.out.protein_lengths)
            .set { protein_lengths }

        branch_by_filtering_status(GENERATE_ANNOTATION_TABLES.out.intron_lengths)
            .set { intron_lengths }

        // Group tables by filtering status and prepare them for comparison
        prepare_for_comparison(summaries.before_contamination_filtering, 'summary').set { summaries_before }
        prepare_for_comparison(protein_lengths.before_contamination_filtering, 'proteins').set { protein_lengths_before }
        prepare_for_comparison(intron_lengths.before_contamination_filtering, 'introns').set { intron_lengths_before }
        collect_and_merge_without_contamination_tables(summaries.before_contamination_filtering, summaries.after_contamination_filtering, 'summary').set { summaries_after }        
        collect_and_merge_without_contamination_tables(protein_lengths.before_contamination_filtering, protein_lengths.after_contamination_filtering, 'proteins').set { protein_lengths_after }
        collect_and_merge_without_contamination_tables(intron_lengths.before_contamination_filtering, intron_lengths.after_contamination_filtering, 'introns').set { intron_lengths_after }

        summaries_before
            .concat(protein_lengths_before, intron_lengths_before, summaries_after, protein_lengths_after, intron_lengths_after)
            .map { status, file -> [status, file, status + '.tsv'] }
            .set { tables_for_merging }

        MERGE_TABLES(tables_for_merging)
        MERGE_TABLES.out
            .filter {it[0] == 'after-contamination-filtering-summary'}
            .set { merged_after_summary }

        SELECT_BEST_ANNOTATION ( merged_after_summary )

        SELECT_BEST_ANNOTATION.out
            .flatMap { status, file ->
                file.text.readLines()
                    .drop(1)
                    .findAll { it.trim() }
                    .collect { line ->
                        def cols = line.split('\t', -1)

                        tuple(
                            cols[0],  // genome_id
                            cols[1],  // method_annotation
                            cols[2],  // sufficient_quality_for_orthoinspector
                            cols[3]   // sufficient_quality_for_orthofinder
                        )
                    }
            }
            .set { best_annotations_for_each_genomes }

        best_annotations_for_each_genomes
            .map { id, annotation_method, sufficient_quality_for_orthoinspector, sufficient_quality_for_orthofinder -> [id, annotation_method] }
            .set { all_best_annotations }

        best_annotations_for_each_genomes
            .filter { id, annotation_method, sufficient_quality_for_orthoinspector, sufficient_quality_for_orthofinder -> sufficient_quality_for_orthoinspector == "yes" }
            .map { id, annotation_method, sufficient_quality_for_orthoinspector, sufficient_quality_for_orthofinder -> [id, annotation_method] }
            .set { best_annotations_for_orthoinspector }
        
        // all annotations
        collect_best_annotation_tables(all_best_annotations, summaries, protein_lengths, intron_lengths, 'all-best-annotations')
            .set { bests_tables_all }

        // best annotations for orthoinspector
        collect_best_annotation_tables(best_annotations_for_orthoinspector, summaries, protein_lengths, intron_lengths, 'best-annotations-for-orthoinspector')
            .set { bests_tables_for_orthoinspector }

        bests_tables_all
            .concat(bests_tables_for_orthoinspector) 
            .set { bests_tables }

        MERGE_BEST_TABLES(bests_tables)
        
        MERGE_TABLES.out
            .concat(MERGE_BEST_TABLES.out)
            .map { id, file -> 
                def items = id.tokenize('-')
                def filtering_status = items[0..-2].join('-')
                def type = items[-1]

                tuple(filtering_status, type, file)
            }
            .groupTuple(by: 0)
            .map { filtering_status, types, files ->
                def file_by_type = [:]

                types.eachWithIndex { type, i ->
                    file_by_type[type] = files[i]
                }

                tuple(
                    filtering_status,
                    file_by_type['summary'],
                    file_by_type['proteins'],
                    file_by_type['introns']
                )
            }
            .combine(newick)
            .set{ annotation_tables_for_comparison }

        PLOT_ANNOTATION_COMPARISON(annotation_tables_for_comparison)

        annotation_tables_for_comparison
            .filter {it[0] == 'all-best-annotations' }
            .combine(all_genome_stats.map {id, file -> file} )
            .set {annotation_tables_for_analysis}

        PLOT_ANNOTATION_ANALYSIS(annotation_tables_for_analysis)


        // best_annotations_for_each_genomes.view { "BEST_ALL = $it" }

        // best_annotations_for_orthoinspector.view { "BEST_ORTHOINSPECTOR = $it" }

        // TRANSFER_FUNCTIONAL_ANNOTATIONS.out.view { id, meta, method, gff ->
        //     "GFF_BEFORE = ${id} ${method} ${gff}"
        // }

        // FILTER_OMARK_CONTAMINATIONS_IN_GFF.out.view { id, meta, method, threshold, gff ->
        //     "GFF_AFTER = ${id} ${method} c=${threshold} ${gff}"
        // }

        collect_best_file(
            all_best_annotations,
            TRANSFER_FUNCTIONAL_ANNOTATIONS.out,
            FILTER_OMARK_CONTAMINATIONS_IN_GFF.out
        ).set { best_gff_all }

        collect_best_file(
            best_annotations_for_orthoinspector,
            TRANSFER_FUNCTIONAL_ANNOTATIONS.out,
            FILTER_OMARK_CONTAMINATIONS_IN_GFF.out
        ).set { best_gff_for_orthoinspector }

        collect_best_file(
            best_annotations_for_orthoinspector,
            TRANSFER_FUNCTIONAL_ANNOTATIONS.out,
            FILTER_OMARK_CONTAMINATIONS_IN_GFF.out
        ).set { best_gff_for_orthoinspector }
    

        collect_best_file(
            all_best_annotations,
            EXTRACT_CANONICAL_PROTEINS.out.fasta,
            FILTER_OMARK_CONTAMINATIONS.out.fasta
        ).set { best_proteins_all }

        collect_best_file(
            best_annotations_for_orthoinspector,
            EXTRACT_CANONICAL_PROTEINS.out.fasta,
            FILTER_OMARK_CONTAMINATIONS.out.fasta
        ).set { best_proteins_for_orthoinspector }

        collect_best_file(
            best_annotations_for_orthoinspector,
            EXTRACT_CANONICAL_PROTEINS.out.fasta,
            FILTER_OMARK_CONTAMINATIONS.out.fasta
        ).set { best_proteins_for_orthoinspector }

        PUBLISH_BEST_ANNOTATIONS(best_gff_for_orthoinspector.concat(best_proteins_for_orthoinspector))
        
/////////////////
        // FILTER_OMARK_CONTAMINATIONS.out.fasta
        //     .map { id, meta, annotation_method, fasta, threshold -> [id, meta, annotation_method + "_c-" + threshold, fasta] }
        //     .join(FILTER_OMARK_CONTAMINATIONS_IN_GFF.out, by: [0, 1, 2])
        //     .join(BUSCO_ON_CONTAMINATION_FILTERED.out.json, by: [0, 1, 2])
        //     .join(OMARK_ON_CONTAMINATION_FILTERED.out.ump, by: [0, 1, 2])
        //     .set { compare_annotations_after_contamination_filtering_input }
        
        // COMPAPARE_ANNOTATIONS(compare_annotations_after_contamination_filtering_input)
        // // If repeat annotation is available, classify genes as overlapping or not with repeats
        // if (params.module_repeats == true) {
        //     // Analysis overlapping between genes annotations and repeats
        //     EXTRACT_CANONICAL_PROTEINS.out.gff
        //         .join(repeats_gff, by: [0, 1]) // -> id, meta, annotation_method, proteins_gff, repeat_gff
        //         .set { genes_vs_repeats_input }

        //     CLASSIFY_GENES_VS_REPEATS(genes_vs_repeats_input)
        // }

        // SIGNALP(EXTRACT_CANONICAL_PROTEINS.out.fasta)
        // EGGNOG_MAPPER(EXTRACT_CANONICAL_PROTEINS.out.fasta)
        // KOFAM_KOALA(EXTRACT_CANONICAL_PROTEINS.out.fasta)
    // emit:
}


