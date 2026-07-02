process VARUS {

    tag "${id}"
    label 'exogap_varus'
    cpus 20
    memory '40 GB'
    scratch false

    input:
    tuple val(id), val(meta), path(genome), path(runlist)

    output:
    tuple val(id), val(meta), path("${id}_varus.bam"), emit: bam

    script:
    """
    mkdir -p varus

    VARUS \\
        -g ${genome} \\
        -r ${runlist} \\
        -o varus \\
        -t ${task.cpus}

    cp varus/*.bam ${id}_varus.bam

    test -s ${id}_varus.bam
    """

    stub:
    """
    touch ${id}_varus.bam
    """
}