process BUILD_BLASTP_TAXONOMY {
    tag "all"
    label 'exogap_python'
    cpus 1
    time '6h'

    input:
    path(blastp_files, stageAs: "input/*")

    output:
    path("all_lineages.tax"), emit: taxonomy
    path("all_lineages.tax.failed_taxids.tsv"), emit: failed_taxids, optional: true

    script:
    """
    blast_taxonomy.py \
        --mode 1 \
        -i input/* \
        -o all_lineages.tax \
        --batch-size 500
    """

    stub:
    """
    touch all_lineages.tax
    """
}