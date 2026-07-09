process PLOT_ANNOTATION_ANALYSIS {
    tag "${status}"
    label 'exogap_python'
    scratch false

    input:
    tuple val(status), path(summary), path(protein_lengths), path(intron_lengths), path(newick), path(seqkit_stats)

    output:
    tuple val(status), path("summary_plots/gene_structure_summary.{png,pdf,svg}"), emit: structure
    tuple val(status), path("correlation_plots/*.{png,pdf,svg}"), emit: correlations

    script:
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib
    mkdir -p "\$MPLCONFIGDIR"

    plot_annotation_analysis.py \
        --summary $summary \
        --intron-lengths $intron_lengths \
        --genome-stats ${seqkit_stats} \
        --tree ${newick} \
        --out-dir .
    """

    stub:
    """
    mkdir -p summary_plots correlation_plots
    touch summary_plots/gene_structure_summary.{png,pdf,svg}
    touch correlation_plots/test.{png,pdf,svg}
    """
}