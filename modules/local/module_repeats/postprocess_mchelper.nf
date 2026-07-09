process POSTPROCESS_MCHELPER {
    tag "${id}"
    label "exogap_python"

    input:
    tuple val(id), val(meta), path(table), path(rm2_families), path(mchelper_classif), path(classified_module), path(unclassified_module), path(gene_tbl), path(rrna_tbl), path(classification)

    output:
    tuple val(id), val(meta), path("${id}_families_processed.tsv"), emit: table
    tuple val(id), val(meta), path("${id}_families_processed.fa"), emit: families
    tuple val(id), val(meta), path("${id}_mchelper_postprocess.log"), emit: log

    script:
    """
    postprocess_mchelper.py \
        --table ${table} \
        --fasta ${rm2_families} \
        --mchelper-classified-fasta ${classified_module} \
        --mchelper-unclassified-fasta ${unclassified_module} \
        --mchelper-classif ${mchelper_classif} \
        --hmmscan-genes ${gene_tbl} \
        --hmmscan-rrna ${rrna_tbl} \
        --classification ${classification} \
        --out ${id}_families_processed > ${id}_mchelper_postprocess.log
    """

    stub:
    """
    touch ${id}_families_processed.tsv ${id}_families_processed.fa ${id}_mchelper_postprocess.log
    """
}
