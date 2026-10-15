process SIGNALP {
    tag "${meta.annotation_method}/${meta.id}"
    scratch false

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("signalp_output"), emit: outdir
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_signalp_summary.tsv"), emit: summary

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

    cp signalp_output/prediction_results.txt ${meta.id}_${meta.annotation_method}_signalp_summary.tsv
    test -s ${meta.id}_${meta.annotation_method}_signalp_summary.tsv
    """
    stub:
    """
    mkdir -p signalp_output
    echo 'ID\tPrediction' > signalp_output/prediction_results.txt
    cp signalp_output/prediction_results.txt ${meta.id}_${meta.annotation_method}_signalp_summary.tsv
    """
}
