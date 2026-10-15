process AGAT_STATISTICS {
    tag "${meta.annotation_method}/${meta.id}"
    label 'agat'

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_stats.out"), emit: stats
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_stats.out_distribution_plots"), emit: plots

    script:
    """
    agat_sp_statistics.pl \
        --gff ${gff} \
        --output ${meta.id}_${meta.annotation_method}_stats.out \
        -d
    """

    stub:
    """
    touch ${meta.id}_${meta.annotation_method}_stats.out
    mkdir -p ${meta.id}_${meta.annotation_method}_stats.out_distribution_plots
    """
}