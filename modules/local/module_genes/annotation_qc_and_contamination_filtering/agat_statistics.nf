process AGAT_STATISTICS {
    tag "${annotation_method}/${id}"
    // label 'exogap_python'
    stageInMode 'copy'
    scratch false

    input:
    tuple val(id), val(meta), val(annotation_method), path(gff)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_stats.out"), emit: stats
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_stats.out_distribution_plots/*/*pdf"), emit: plots

 
    script:
    """
    module load agat
    agat_sp_statistics.pl --gff $gff -o ${id}_${annotation_method}_stats.out -d
    """

    stub:
    """
    touch ${id}_${annotation_method}_stats.out
    """
}