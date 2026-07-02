process SAMTOOLS_READ_COUNT {

    tag "${id}/${aligner}"
    label 'exogap_samtools'
    scratch false

    input:
    tuple val(id), val(meta), val(aligner), path(bam)

    output:
    tuple val(id), val(meta), val(aligner), path("${id}_${aligner}_mapped_reads.txt"), emit: out

    script:
    """
    samtools view \\
        -@ ${task.cpus} \\
        -F 4 \\
        -c ${bam} > ${id}_${aligner}_mapped_reads.txt

    test -s ${id}_${aligner}_mapped_reads.txt
    """

    stub:
    """
    echo 10000000 > ${id}_${aligner}_mapped_reads.txt
    """

}