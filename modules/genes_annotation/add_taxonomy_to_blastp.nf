process ADD_TAXONOMY_TO_BLASTP {
    tag "${annotation_method}/${id}"
    label 'exogap_python'
    cpus 1
    time '6h'

    input:
    tuple val(id), val(meta), val(annotation_method), path(blastp)
    each path(taxonomy)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_with_tax.tsv")

    script:
    """
    blast_taxonomy.py \
        --mode 2 \
        -i ${blastp} \
        -o ${id}_${annotation_method}_with_tax.tsv \
        --taxonomy ${taxonomy}
    """

    stub:
    """
    touch ${id}_${annotation_method}_with_tax.tsv
    """
}