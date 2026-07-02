process SIGNALP {
    tag "${id}"
    scratch false
    memory { (10 * task.attempt).GB }

    input:
    tuple val(id), val(meta), path(fasta)

    output:
    tuple val(id), val(meta), path("signalp_output"), emit: outdir
    tuple val(id), val(meta), path("${id}_signalp_summary.tsv"), emit: summary

    script:
    """
    export APPTAINER_HOME="\$PWD"
    export SINGULARITY_HOME="\$PWD"
    export APPTAINER_CACHEDIR="\$PWD/.apptainer"
    export SINGULARITY_CACHEDIR="\$PWD/.singularity"
    mkdir -p "\$APPTAINER_CACHEDIR" "\$SINGULARITY_CACHEDIR"

    export HOME=\$PWD
    export TMPDIR=\$PWD/tmp
    export TMP=\$PWD/tmp
    export TEMP=\$PWD/tmp
    mkdir -p "\$TMPDIR"

    module load signalp

    signalp6 \\
        --fastafile ${fasta} \\
        --organism eukarya \\
        --output_dir signalp_output \\
        --format none \\
        --mode fast

    cp signalp_output/prediction_results.txt ${id}_signalp_summary.tsv
    test -s ${id}_signalp_summary.tsv
    """
}