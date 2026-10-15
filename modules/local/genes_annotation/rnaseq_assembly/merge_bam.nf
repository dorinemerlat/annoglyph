process MERGE_BAM {
    tag "${meta.id}"
    label 'varus'

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("${meta.id}.bam")

    script:
    """
    samtools merge \
        -@ ${task.cpus} \
        -o ${meta.id}.bam \
        ${bams}
    """

    stub:
    """
    command -v samtools >/dev/null

    touch ${meta.id}.bam
    """
}
