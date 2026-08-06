process TRANSFER_FUNCTIONAL_ANNOTATIONS {
    tag "${annotation_method}/${id}"
    label 'exogap_python'

    input:
    tuple val(id), val(meta), val(annotation_method), path(gff), path(blastp), path(interproscan), path(eggnog)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_with_functional_annotations.gff")

    script:
    def output = "${id}_${annotation_method}_with_functional_annotations.gff"

    """
    # Transfer functional annotations to the GFF3.
    transfer_functional_annotations.py \
        -g ${gff} \
        -b ${blastp} \
        -i ${interproscan} \
        -e ${eggnog} \
        -o ${output}
    """

    stub:
    """
    command -v transfer_functional_annotations.py >/dev/null

    touch ${output}
    """
}
