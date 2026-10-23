process PLOT_ANNOTATION_COMPARISON {
    tag "${status}"
    label 'annoglyph_python'

    input:
    tuple val(status), path(summary), path(protein_lengths), path(intron_lengths), path(newick), path(taxonomy)

    output:
    tuple val(status), path("${status}_annotation_comparison_summary.{png,pdf,svg}"), emit: comparison

    script:
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib
    mkdir -p "\$MPLCONFIGDIR"

    plot_annotation_comparison.py \
        --summary ${summary} \
        --protein-lengths ${protein_lengths} \
        --taxonomy ${taxonomy} --tree ${newick} \
        --out ${status}_annotation_comparison_summary
    """

    stub:
    """
    command -v plot_annotation_comparison.py >/dev/null

    touch \
        ${status}_annotation_comparison_summary.png \
        ${status}_annotation_comparison_summary.pdf \
        ${status}_annotation_comparison_summary.svg
    """
}