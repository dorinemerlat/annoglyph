process DEEPMHMM {
    tag "${id}"
    label 'exogap_deepmhmm'

    input:
    tuple val(id), val(meta), path(fasta)

    output:
    tuple val(id), val(meta), path("deepmhmm_output"), emit: outdir
    tuple val(id), val(meta), path("${id}_deepmhmm.tsv"), emit: tsv

    script:
    """
    biolib run DTU/DeepTMHMM \\
        --fasta ${fasta}

    mkdir -p deepmhmm_output
    cp biolib_results/* deepmhmm_output/ || true

    if [ -s biolib_results/TMRs.gff3 ]; then
        cp biolib_results/TMRs.gff3 ${id}_deepmhmm.tsv
    elif [ -s biolib_results/predicted_topologies.3line ]; then
        cp biolib_results/predicted_topologies.3line ${id}_deepmhmm.tsv
    else
        touch ${id}_deepmhmm.tsv
    fi
    """

    stub:
    """
    mkdir -p deepmhmm_output
    touch ${id}_deepmhmm.tsv
    """
}