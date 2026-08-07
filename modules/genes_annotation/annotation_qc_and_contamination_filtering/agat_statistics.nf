process AGAT_STATISTICS {
    tag "${annotation_method}/${id}"
    label 'agat'

    input:
    tuple val(id), val(meta), val(annotation_method), path(gff)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_stats.out"), emit: stats
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_stats.out_distribution_plots"), emit: plots

    script:
    """
    agat_sp_statistics.pl \
        --gff ${gff} \
        --output ${id}_${annotation_method}_stats.out \
        -d
    """

    stub:
    """
    touch ${id}_${annotation_method}_stats.out
    mkdir -p ${id}_${annotation_method}_stats.out_distribution_plots
    """
}