process TRANSFER_FUNCTIONAL_ANNOTATIONS {

    tag "${id}/${annotation_method}"
    label 'exogap_python'
    scratch false

    input:
    tuple val(id), val(meta), val(annotation_method), path(gff), path(blastp), path(interproscan), path(eggnog)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_with_functional_annotations.gff"), emit: gff

    script:
    """
    transfer_functional_annotations.py \\
        --gff ${gff} \\
        --blast ${blastp} \\
        --interpro ${interproscan} \\
        --eggnog ${eggnog} \\
        --output ${id}_${annotation_method}_with_functional_annotations.gff
    """

    stub:
    """
    touch ${id}_${annotation_method}_with_functional_annotations.gff
    """
}