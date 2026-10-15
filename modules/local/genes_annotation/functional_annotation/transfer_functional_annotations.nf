process TRANSFER_FUNCTIONAL_ANNOTATIONS {
    tag "${meta.annotation_method}/${meta.id}"
    label 'exogap_python'

    input:
    tuple val(meta), path(gff), path(blastp), path(interproscan), path(eggnog)

    output:
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_with_functional_annotations.gff")

    script:
    def output = "${meta.id}_${meta.annotation_method}_with_functional_annotations.gff"

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
    def output = "${meta.id}_${meta.annotation_method}_with_functional_annotations.gff"
    """
    command -v transfer_functional_annotations.py >/dev/null

    touch ${output}
    """
}
