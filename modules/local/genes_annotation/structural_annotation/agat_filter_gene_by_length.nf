process AGAT_FILTER_GENE_BY_LENGTH {
    tag "${meta.annotation_method}/${size}/${meta.id}"
    label 'agat'

    input:
    tuple val(meta), path(gff), val(size)

    output:
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_size${size}.gff")

    script:
    """
    agat_sp_filter_gene_by_length.pl \
        -f ${gff} \
        -o ${meta.id}_${meta.annotation_method}_size${size}.gff \
        --test ">=" \
        -s ${size}
    """

    stub:
    """
    command -v agat_sp_filter_gene_by_length.pl >/dev/null

    touch ${meta.id}_${meta.annotation_method}_size${size}.gff
    """
}
