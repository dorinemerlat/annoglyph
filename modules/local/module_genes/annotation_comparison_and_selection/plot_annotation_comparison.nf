process PLOT_ANNOTATION_COMPARISON {
    tag "${status}"
    label 'exogap_python'
    scratch false

    input:
    tuple val(status), path(summary), path(protein_lengths), path(intron_lengths), path(newick)

    output:
    tuple val(status), path("${status}_annotation_comparison_summary.{png,pdf,svg}"), emit: comparison

    script:
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib
    mkdir -p "\$MPLCONFIGDIR"

    plot_annotation_comparison.py \\
        --summary ${summary} \\
        --protein-lengths ${protein_lengths} \\
        --out "${status}_annotation_comparison_summary" \\
        --tree ${newick}
    """

    stub:
    """
    mkdir -p figures/summary_plots figures/correlation_plots
    touch ${status}_annotation_comparison_summary.{png,pdf,svg}
    """
}