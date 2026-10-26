process PLOT_ANNOTATION_ANALYSIS {
    tag "${status}"
    label 'annoglyph_python'

    input:
    tuple val(status), path(summary), path(intron_lengths), path(newick), path(seqkit_stats), path(taxonomy), path(selection)

    output:
    tuple val(status), path("summary_plots/gene_structure_summary.{png,pdf,svg}"), emit: structure
    tuple val(status), path("correlation_plots/*.{png,pdf,svg}"), emit: correlations

    script:
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib
    mkdir -p "\$MPLCONFIGDIR"

    plot_annotation_analysis.py \
        --summary ${summary} \
        --intron-lengths ${intron_lengths} \
        --genome-stats ${seqkit_stats} \
        --taxonomy ${taxonomy} --tree ${newick} --selection ${selection} \
        --out-dir .
    """

    stub:
    """
    command -v plot_annotation_analysis.py >/dev/null

    mkdir -p summary_plots correlation_plots

    touch \
        summary_plots/gene_structure_summary.png \
        summary_plots/gene_structure_summary.pdf \
        summary_plots/gene_structure_summary.svg \
        correlation_plots/test.png \
        correlation_plots/test.pdf \
        correlation_plots/test.svg
    """
}
