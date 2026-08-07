// Include local modules
include { AGAT_STATISTICS                            } from "$projectDir/modules/genes_annotation/annotation_qc_and_contamination_filtering/agat_statistics.nf"
include { DOWNLOAD_OMAMER_DB                         } from "$projectDir/modules/genes_annotation/annotation_qc_and_contamination_filtering/download_omamer_db.nf"
include { OMARK                                      } from "$projectDir/modules/genes_annotation/annotation_qc_and_contamination_filtering/omark.nf"
include { BUSCO                                      } from "$projectDir/modules/genes_annotation/annotation_qc_and_contamination_filtering/busco.nf"
include { FILTER_OMARK_CONTAMINATIONS                } from "$projectDir/modules/genes_annotation/annotation_qc_and_contamination_filtering/filter_omark_contaminations.nf"
include { FILTER_OMARK_CONTAMINATIONS_IN_GFF         } from "$projectDir/modules/genes_annotation/annotation_qc_and_contamination_filtering/filter_omark_contaminations_in_gff.nf"
include { OMARK as OMARK_ON_CONTAMINATION_FILTERED   } from "$projectDir/modules/genes_annotation/annotation_qc_and_contamination_filtering/omark.nf"
include { BUSCO as BUSCO_ON_CONTAMINATION_FILTERED   } from "$projectDir/modules/genes_annotation/annotation_qc_and_contamination_filtering/busco.nf"


workflow ANNOTATION_QC_AND_CONTAMINATION_FILTERING {

    take:
        gff
        proteins

    main:
        contamination_threshold = params.omark_contamination_threshold

        // Compute annotation statistics.
        AGAT_STATISTICS(gff)

        // Run OMArk and BUSCO before filtering.
        DOWNLOAD_OMAMER_DB()

        proteins
            .combine(DOWNLOAD_OMAMER_DB.out)
            .set { omark_input }

        OMARK(omark_input)
        BUSCO(proteins)

        // Separate contaminated and clean annotations.
        OMARK.out.contaminants_count
            .branch {
                contaminated: it[3].text.trim().toInteger() > 0
                clean: true
            }
            .set { contamination_status }

        // Prepare contamination filtering.
        contamination_status.contaminated
            .join(gff, by: [0, 1, 2])
            .join(proteins, by: [0, 1, 2])
            .join(OMARK.out.ump, by: [0, 1, 2])
            .map { id, meta, annotation_method, contaminants_count, gff, fasta, ump ->
                [id, meta, annotation_method, gff, fasta, ump, contamination_threshold]
            }
            .set { contamination_filtering_input }

        FILTER_OMARK_CONTAMINATIONS(contamination_filtering_input)

        // Remove contaminants from GFF files.
        FILTER_OMARK_CONTAMINATIONS.out.contaminations
            .join(gff, by: [0, 1, 2])
            .set { gff_filtering_input }

        FILTER_OMARK_CONTAMINATIONS_IN_GFF(gff_filtering_input)

        // Reassess filtered annotations.
        BUSCO_ON_CONTAMINATION_FILTERED(FILTER_OMARK_CONTAMINATIONS.out.fasta)

        FILTER_OMARK_CONTAMINATIONS.out.fasta
            .combine(DOWNLOAD_OMAMER_DB.out)
            .set { omark_filtered_input }

        OMARK_ON_CONTAMINATION_FILTERED(omark_filtered_input)

        // Build QC metrics before filtering.
        proteins
            .join(gff, by: [0, 1, 2])
            .join(BUSCO.out.json, by: [0, 1, 2])
            .join(OMARK.out.sum, by: [0, 1, 2])
            .map { id, meta, annotation_method, fasta, gff, busco_json, omark_sum ->
                [id, meta, annotation_method, 'before-contamination-filtering', fasta, gff, busco_json, omark_sum]
            }
            .set { qc_metrics_before_filtering }

        // Build QC metrics for filtered annotations.
        FILTER_OMARK_CONTAMINATIONS.out.fasta
            .join(FILTER_OMARK_CONTAMINATIONS_IN_GFF.out, by: [0, 1, 2])
            .join(BUSCO_ON_CONTAMINATION_FILTERED.out.json, by: [0, 1, 2])
            .join(OMARK_ON_CONTAMINATION_FILTERED.out.sum, by: [0, 1, 2])
            .map { id, meta, annotation_method, fasta, gff, busco_json, omark_sum ->
                [id, meta + [contamination_filtering: true], annotation_method, 'after-contamination-filtering', fasta, gff, busco_json, omark_sum]
            }
            .set { qc_metrics_after_filtering_filtered }

        // Reuse original metrics for clean annotations.
        contamination_status.clean
            .join(qc_metrics_before_filtering, by: [0, 1, 2])
            .map { id, meta, annotation_method, contaminants_count, filtering_status, fasta, gff, busco_json, omark_sum ->
                [id, meta + [contamination_filtering: false], annotation_method, 'after-contamination-filtering', fasta, gff, busco_json, omark_sum]
            }
            .set { qc_metrics_after_filtering_clean }

        // Complete after-filtering channel.
        qc_metrics_after_filtering_filtered
            .concat(qc_metrics_after_filtering_clean)
            .set { qc_metrics_after_filtering }

        // Build complete final GFF channel.
        FILTER_OMARK_CONTAMINATIONS_IN_GFF.out
            .map { id, meta, annotation_method, gff ->
                [id, meta + [contamination_filtering: true], annotation_method, gff]
            }
            .set { filtered_gff }

        contamination_status.clean
            .join(gff, by: [0, 1, 2])
            .map { id, meta, annotation_method, contaminants_count, gff ->
                [id, meta + [contamination_filtering: false], annotation_method, gff]
            }
            .set { clean_gff }

        filtered_gff
            .concat(clean_gff)
            .set { final_gff }

        // Build complete final protein channel.
        FILTER_OMARK_CONTAMINATIONS.out.fasta
            .map { id, meta, annotation_method, fasta ->
                [id, meta + [contamination_filtering: true], annotation_method, fasta]
            }
            .set { filtered_proteins }

        contamination_status.clean
            .join(proteins, by: [0, 1, 2])
            .map { id, meta, annotation_method, contaminants_count, fasta ->
                [id, meta + [contamination_filtering: false], annotation_method, fasta]
            }
            .set { clean_proteins }

        filtered_proteins
            .concat(clean_proteins)
            .set { final_proteins }

    emit:
        agat_statistics = AGAT_STATISTICS.out.stats

        filtered_gff      = final_gff
        filtered_proteins = final_proteins

        qc_metrics_before_filtering = qc_metrics_before_filtering
        qc_metrics_after_filtering  = qc_metrics_after_filtering
}