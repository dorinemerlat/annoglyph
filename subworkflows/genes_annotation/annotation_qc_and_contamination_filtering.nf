include { AGAT_STATISTICS } from "$projectDir/modules/local/genes_annotation/annotation_qc_and_contamination_filtering/agat_statistics.nf"
include { DOWNLOAD_OMAMER_DB } from "$projectDir/modules/local/genes_annotation/annotation_qc_and_contamination_filtering/download_omamer_db.nf"
include { OMARK; OMARK as OMARK_AFTER } from "$projectDir/modules/local/genes_annotation/annotation_qc_and_contamination_filtering/omark.nf"
include { BUSCO_PROTEINS; BUSCO_PROTEINS as BUSCO_PROTEINS_AFTER } from "$projectDir/modules/local/genes_annotation/annotation_qc_and_contamination_filtering/busco_proteins.nf"
include { FILTER_OMARK_CONTAMINATIONS } from "$projectDir/modules/local/genes_annotation/annotation_qc_and_contamination_filtering/filter_omark_contaminations.nf"
include { FILTER_OMARK_CONTAMINATIONS_IN_GFF } from "$projectDir/modules/local/genes_annotation/annotation_qc_and_contamination_filtering/filter_omark_contaminations_in_gff.nf"

workflow ANNOTATION_QC_AND_CONTAMINATION_FILTERING {
    take:
    gff
    proteins
    busco_lineage
    main:
    if (params.omamer_db) {
        omamer_database = Channel.value(file(params.omamer_db, checkIfExists: true))
    } else {
        DOWNLOAD_OMAMER_DB()
        omamer_database = DOWNLOAD_OMAMER_DB.out.first()
    }
    AGAT_STATISTICS(gff)
    OMARK(proteins.combine(omamer_database))
    BUSCO_PROTEINS(proteins, busco_lineage)

    contamination = OMARK.out.contaminants_count.branch { meta, count ->
        contaminated: count.text.trim().toInteger() > 0
        clean: true
    }
    filtering_inputs = contamination.contaminated
        .join(gff, failOnDuplicate: true, failOnMismatch: false)
        .join(proteins, failOnDuplicate: true, failOnMismatch: false)
        .join(OMARK.out.ump, failOnDuplicate: true, failOnMismatch: false)
        .map { meta, count, gff_file, fasta, ump -> [meta, gff_file, fasta, ump, params.omark_contamination_threshold] }
    FILTER_OMARK_CONTAMINATIONS(filtering_inputs)
    FILTER_OMARK_CONTAMINATIONS_IN_GFF(FILTER_OMARK_CONTAMINATIONS.out.contaminations.join(gff))

    // Do not run protein QC on an empty FASTA after removing contaminants.
    viable_filtered = FILTER_OMARK_CONTAMINATIONS.out.fasta.filter { meta, fasta ->
        def valid = fasta.size() > 0
        if (!valid) log.warn "No proteins remain after OMArk filtering: ${meta.id}/${meta.annotation_method}"
        valid
    }
    BUSCO_PROTEINS_AFTER(viable_filtered, busco_lineage)
    OMARK_AFTER(viable_filtered.combine(omamer_database))

    before_metrics = proteins
        .join(gff, failOnDuplicate: true, failOnMismatch: true)
        .join(BUSCO_PROTEINS.out.json, failOnDuplicate: true, failOnMismatch: true)
        .join(OMARK.out.sum, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, fasta, annotation, busco, omark ->
            [meta + [filtering_status: 'before-contamination-filtering'], fasta, annotation, busco, omark]
        }
    after_filtered = viable_filtered
        .join(FILTER_OMARK_CONTAMINATIONS_IN_GFF.out)
        .join(BUSCO_PROTEINS_AFTER.out.json)
        .join(OMARK_AFTER.out.sum)
        .map { meta, fasta, annotation, busco, omark ->
            [meta + [filtering_status: 'after-contamination-filtering', contamination_filtering: true], fasta, annotation, busco, omark]
        }
    after_clean = contamination.clean
        .join(proteins).join(gff).join(BUSCO_PROTEINS.out.json).join(OMARK.out.sum)
        .map { meta, count, fasta, annotation, busco, omark ->
            [meta + [filtering_status: 'after-contamination-filtering', contamination_filtering: false], fasta, annotation, busco, omark]
        }
    after_metrics = after_filtered.concat(after_clean)
    emit:
    database = omamer_database
    agat_statistics = AGAT_STATISTICS.out.stats
    qc_metrics_before_filtering = before_metrics
    qc_metrics_after_filtering = after_metrics
    filtered_gff = after_metrics.map { meta, fasta, annotation, busco, omark -> [meta, annotation] }
    filtered_proteins = after_metrics.map { meta, fasta, annotation, busco, omark -> [meta, fasta] }
}
