process TRANSFER_FUNCTIONAL_ANNOTATIONS {
    tag "${annotation_method}/${id}"
    label 'exogap_python'
    memory '10GB'

    input:
    tuple val(id), val(meta), val(annotation_method), path(gff), path(blastp), path(interproscan), path(eggnog)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_with_functional_annotations.gff")

    script:
    """
    transfer_functional_annotations.py -g ${gff} -b ${blastp} -i ${interproscan} -e ${eggnog} -o ${id}_${annotation_method}_with_functional_annotations.gff
    """

    stub:
    """
    touch ${id}_${annotation_method}_with_functional_annotations.gff 
    """
}