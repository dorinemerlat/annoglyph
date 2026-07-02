// Include local modules
include { AGAT_STATISTICS                          } from "$projectDir/modules/local/module_genes/annotation_qc_and_contamination_filtering/agat_statistics.nf"
include { DOWNLOAD_OMAMER_DB                       } from "$projectDir/modules/local/module_genes/annotation_qc_and_contamination_filtering/download_omamer_db.nf"
include { OMARK                                    } from "$projectDir/modules/local/module_genes/annotation_qc_and_contamination_filtering/omark.nf"
include { BUSCO                                    } from "$projectDir/modules/local/module_genes/annotation_qc_and_contamination_filtering/busco.nf"
include { FILTER_OMARK_CONTAMINATIONS              } from "$projectDir/modules/local/module_genes/annotation_qc_and_contamination_filtering/filter_omark_contaminations.nf"
include { FILTER_OMARK_CONTAMINATIONS_IN_GFF       } from "$projectDir/modules/local/module_genes/annotation_qc_and_contamination_filtering/filter_omark_contaminations_in_gff.nf"
include { OMARK as OMARK_ON_CONTAMINATION_FILTERED } from "$projectDir/modules/local/module_genes/annotation_qc_and_contamination_filtering/omark.nf"
include { BUSCO as BUSCO_ON_CONTAMINATION_FILTERED } from "$projectDir/modules/local/module_genes/annotation_qc_and_contamination_filtering/busco.nf"


workflow ANNOTATION_QC_AND_CONTAMINATION_FILTERING {

    take:
        /*
         * GFF files enriched with functional annotations.
         *
         * Expected format:
         *   [id, meta, annotation_method, gff]
         */
        gff

        /*
         * Canonical protein FASTA files before contamination filtering.
         *
         * Expected format:
         *   [id, meta, annotation_method, fasta]
         */
        proteins

    main:

        /*
         * Contamination filtering threshold.
         *
         * Define in config:
         *   params.omark_contamination_threshold ->  default: 0.1
         */
        contamination_threshold = params.omark_contamination_threshold.toString()


        /*
         * Compute basic annotation statistics before contamination filtering.
         */
        AGAT_STATISTICS(gff)


        /*
         * Run OMARK and BUSCO before contamination filtering.
         */
        DOWNLOAD_OMAMER_DB()

        proteins
            .combine(DOWNLOAD_OMAMER_DB.out)
            .set { omark_input }

        OMARK(omark_input)
        BUSCO(proteins)


        /*
         * Prepare contamination filtering only for annotations with OMARK contaminants.
         *
         * FILTER_OMARK_CONTAMINATIONS expects:
         *   [id, meta, annotation_method, gff, fasta, ump, threshold]
         */
        OMARK.out.contaminants_count
            .filter { id, meta, annotation_method, contaminants_count ->
                contaminants_count.readLines()[0].toInteger() > 0
            }
            .join(gff, by: [0, 1, 2])
            .join(proteins, by: [0, 1, 2])
            .join(OMARK.out.ump, by: [0, 1, 2])
            .map { id, meta, annotation_method, contaminants_count, gff_file, fasta, ump ->
                [ id, meta, annotation_method, gff_file, fasta, ump, contamination_threshold ]
            }
            .set { contamination_filtering_input }

        FILTER_OMARK_CONTAMINATIONS(contamination_filtering_input)


        /*
         * Remove OMARK contaminants from functional GFF annotations.
         *
         * Expected input:
         *   [id, meta, annotation_method, contaminations, gff]
         */
        FILTER_OMARK_CONTAMINATIONS.out.contaminations
            .combine(gff, by: [0, 1, 2])
            .set { filter_omark_contaminations_in_gff_input }

        FILTER_OMARK_CONTAMINATIONS_IN_GFF(filter_omark_contaminations_in_gff_input)


        /*
         * Prepare filtered protein FASTA files for BUSCO and OMARK.
         *
         * Expected input:
         *   [id, meta, annotation_method, fasta]
         */

        BUSCO_ON_CONTAMINATION_FILTERED(FILTER_OMARK_CONTAMINATIONS.out.fasta)

        FILTER_OMARK_CONTAMINATIONS.out.fasta
            .combine(DOWNLOAD_OMAMER_DB.out)
            .set { omark_filtered_input }

        OMARK_ON_CONTAMINATION_FILTERED(omark_filtered_input)


        /*
         * Build comparison input before contamination filtering.
         *
         * Output format:
         *   [id, meta, annotation_method, filtering_status, fasta, gff, busco_json, omark_sum]
         */
        proteins
            .join(gff, by: [0, 1, 2])
            .join(BUSCO.out.json, by: [0, 1, 2])
            .join(OMARK.out.sum, by: [0, 1, 2])
            .map { id, meta, annotation_method, fasta, gff_file, busco_json, omark_sum ->
                [ id, meta, annotation_method, 'before-contamination-filtering', fasta, gff_file, busco_json, omark_sum ]
            }
            .set { qc_metrics_before_filtering }


        /*
        * Keys present after contamination filtering.
        */
        FILTER_OMARK_CONTAMINATIONS.out.fasta
            .map { id, meta, annotation_method, fasta ->
                tuple(id, annotation_method)
            }
            .set { filtered_keys }


        /*
        * Build comparison input after contamination filtering.
        *
        * For annotations without OMARK contamination filtering, reuse the before-filtering files.
        */
        FILTER_OMARK_CONTAMINATIONS.out.fasta
            .join(FILTER_OMARK_CONTAMINATIONS_IN_GFF.out, by: [0, 1, 2])
            .join(BUSCO_ON_CONTAMINATION_FILTERED.out.json, by: [0, 1, 2])
            .join(OMARK_ON_CONTAMINATION_FILTERED.out.sum, by: [0, 1, 2])
            .map { id, meta, annotation_method, fasta, gff_file, busco_json, omark_sum ->
                tuple(id, meta, annotation_method, 'after-contamination-filtering', fasta, gff_file, busco_json, omark_sum)
            }
            .set { qc_metrics_after_filtering_filtered_only }


        /*
        * Recover annotations that were not filtered.
        *
        * These are present before filtering but absent from the filtered outputs.
        * For them, before and after are identical.
        */
        qc_metrics_before_filtering
            .map { id, meta, annotation_method, filtering_status, fasta, gff_file, busco_json, omark_sum ->
                tuple(id, annotation_method, meta, fasta, gff_file, busco_json, omark_sum)
            }
            .join(filtered_keys, by: [0, 1], remainder: true)
            .filter { id, annotation_method, meta, fasta, gff_file, busco_json, omark_sum, filtered_key ->
                filtered_key == null
            }
            .map { id, annotation_method, meta, fasta, gff_file, busco_json, omark_sum, filtered_key ->
                tuple(id, meta, annotation_method, 'after-contamination-filtering', fasta, gff_file, busco_json, omark_sum)
            }
            .set { qc_metrics_after_filtering_unfiltered }


        qc_metrics_after_filtering_filtered_only
            .concat(qc_metrics_after_filtering_unfiltered)
            .set { qc_metrics_after_filtering }


    emit:
        /*
         * AGAT statistics before contamination filtering.
         */
        agat_statistics = AGAT_STATISTICS.out.stats

        /*
         * Contamination-filtered outputs.
         */
        filtered_gff      = FILTER_OMARK_CONTAMINATIONS_IN_GFF.out
        filtered_proteins = FILTER_OMARK_CONTAMINATIONS.out.fasta

        /*
         * Standardized inputs for annotation comparison.
         */
        qc_metrics_before_filtering = qc_metrics_before_filtering
        qc_metrics_after_filtering  = qc_metrics_after_filtering
}