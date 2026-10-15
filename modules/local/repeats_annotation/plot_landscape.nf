process PLOT_REPEAT_LANDSCAPE {
    tag "${meta.id}"
    label 'exogap_python'
    input:
    tuple val(meta), path(table)
    output:
    tuple val(meta), path('repeat_landscape.{png,pdf,svg}'), emit: plots
    script:
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib
    plot_repeat_landscape.py --input ${table} --output repeat_landscape
    """
    stub:
    """
    touch repeat_landscape.png repeat_landscape.pdf repeat_landscape.svg
    """
}
