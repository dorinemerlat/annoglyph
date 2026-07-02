process AGAT_FILTER_GENE_BY_LENGTH {

    tag "${id}/${annotation_method}"
    label 'exogap_agat'
    scratch false

    input:
    tuple val(id), val(meta), val(annotation_method), path(gff)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_filtered.gff")

    script:
    """
    agat_sp_filter_by_gene_length.pl \
        --gff ${gff} \
        --size ${params.min_gene_length} \
        --output ${id}_${annotation_method}_filtered.gff

    test -s ${id}_${annotation_method}_filtered.gff
    """

    stub:
    """
    touch ${id}_${annotation_method}_filtered.gff
    """
}