process AGAT_FILTER_GENE_BY_LENGTH {
    tag "${annotation_method}/${size}/${id}"
    label 'agat'

    input:
    tuple val(id), val(meta), val(annotation_method), path(gff), val(size)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_size${size}.gff")

    script:
    """
    agat_sp_filter_gene_by_length.pl \
        -f ${gff} \
        -o ${id}_${annotation_method}_size${size}.gff \
        --test ">=" \
        -s ${size}
    """

    stub:
    """
    command -v agat_sp_filter_gene_by_length.pl >/dev/null

    touch ${id}_${annotation_method}_size${size}.gff
    """
}
