process BLASTP {
    tag "${annotation_method}/${id}"
    label 'blast'
    cpus 10

    input:
    tuple val(id), val(meta), val(annotation_method), path(archive), val(suffix)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${suffix}")

    script:
    """
    blast_formatter -archive $archive -outfmt $format -out ${id}_${suffix}
    """

    stub:
    """
    touch ${id}_${suffix}
    """
}