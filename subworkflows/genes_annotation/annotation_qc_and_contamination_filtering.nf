// Include local modules
include { AGAT_STATISTICS                          } from "$projectDir/modules/genes_annotation/annotation_qc_and_contamination_filtering/agat_statistics.nf"
include { DOWNLOAD_OMAMER_DB                       } from "$projectDir/modules/genes_annotation/annotation_qc_and_contamination_filtering/download_omamer_db.nf"
include { OMARK                                    } from "$projectDir/modules/genes_annotation/annotation_qc_and_contamination_filtering/omark.nf"
include { BUSCO                                    } from "$projectDir/modules/genes_annotation/annotation_qc_and_contamination_filtering/busco.nf"
include { FILTER_OMARK_CONTAMINATIONS              } from "$projectDir/modules/genes_annotation/annotation_qc_and_contamination_filtering/filter_omark_contaminations.nf"
include { FILTER_OMARK_CONTAMINATIONS_IN_GFF       } from "$projectDir/modules/genes_annotation/annotation_qc_and_contamination_filtering/filter_omark_contaminations_in_gff.nf"
include { OMARK as OMARK_ON_CONTAMINATION_FILTERED } from "$projectDir/modules/genes_annotation/annotation_qc_and_contamination_filtering/omark.nf"
include { BUSCO as BUSCO_ON_CONTAMINATION_FILTERED } from "$projectDir/modules/genes_annotation/annotation_qc_and_contamination_filtering/busco.nf"


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
            .map { id, meta, annotation_method, contaminants_count, gff, fasta, ump ->
                [ id, meta, annotation_method, gff, fasta, ump, contamination_threshold ]
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
            .map { id, meta, annotation_method, fasta, gff, busco_json, omark_sum ->
                [ id, meta, annotation_method, 'before-contamination-filtering', fasta, gff, busco_json, omark_sum ]
            }
            .set { qc_metrics_before_filtering }


        /*
        * Build comparison input after contamination filtering.
        *
        * For annotations without OMARK contamination filtering, reuse the before-filtering files.
        */
        FILTER_OMARK_CONTAMINATIONS.out.fasta
            .join(FILTER_OMARK_CONTAMINATIONS_IN_GFF.out, by: [0, 1, 2])
            .join(BUSCO_ON_CONTAMINATION_FILTERED.out.json, by: [0, 1, 2])
            .join(OMARK_ON_CONTAMINATION_FILTERED.out.sum, by: [0, 1, 2])
            .map { id, meta, annotation_method, fasta, gff, busco_json, omark_sum ->
                [ id, meta, annotation_method, 'after-contamination-filtering', fasta, gff, busco_json, omark_sum ]
            }
            .set { qc_metrics_after_filtering_filtered_only }


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
        * Metrics for every annotation before filtering.
        */
        qc_metrics_before_filtering = qc_metrics_before_filtering

        /*
        * Metrics only for annotations that were actually filtered.
        *
        * Annotations without detected contaminants are absent from this channel.
        */
        qc_metrics_after_filtering_filtered_only = qc_metrics_after_filtering_filtered_only
}