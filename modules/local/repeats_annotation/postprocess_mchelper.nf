process POSTPROCESS_MCHELPER {
    tag "${meta.id}"
    label "annoglyph_python"

    input:
    tuple val(meta), path(repeats_table), path(repeats_families), path(mchelper_classif), path(classified_module), path(unclassified_module), path(gene_tbl), path(rrna_tbl), path(classification)

    output:
    tuple val(meta), path("${meta.id}_families_processed.tsv"), emit: table
    tuple val(meta), path("${meta.id}_families_processed.fa"), emit: families
    tuple val(meta), path("${meta.id}_mchelper_postprocess.log"), emit: log

    script:
    """
    postprocess_mchelper.py \
        --table ${repeats_table} \
        --fasta ${repeats_families} \
        --mchelper-classified-fasta ${classified_module} \
        --mchelper-unclassified-fasta ${unclassified_module} \
        --mchelper-classif ${mchelper_classif} \
        --hmmscan-genes ${gene_tbl} \
        --hmmscan-rrna ${rrna_tbl} \
        --classification ${classification} \
        --out ${meta.id}_families_processed > ${meta.id}_mchelper_postprocess.log
    """

    stub:
    """
    touch ${meta.id}_families_processed.tsv
    touch ${meta.id}_families_processed.fa
    touch ${meta.id}_mchelper_postprocess.log
    """
}