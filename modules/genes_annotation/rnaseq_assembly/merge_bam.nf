process MERGE_BAM {
    tag "${id}"
    label 'varus'

    input:
    tuple val(id), val(meta), path(bams)

    output:
    tuple val(id), val(meta), path("${id}.bam")

    script:
    """
    samtools merge \
        -@ ${task.cpus} \
        -o ${id}.bam \
        ${bams}
    """

    stub:
    """
    command -v samtools >/dev/null

    touch ${id}.bam
    """
}
