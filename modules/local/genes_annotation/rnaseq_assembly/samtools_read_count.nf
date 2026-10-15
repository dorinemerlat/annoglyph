process SAMTOOLS_READ_COUNT {
    tag "${meta.rnaseq_aligner}/${meta.id}"
    label 'varus'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}_${meta.rnaseq_aligner}_mapped_read_count.txt")

    script:
    """
    # Count mapped reads.
    samtools view \
        -c \
        -F 4 \
        ${bam} \
        > ${meta.id}_${meta.rnaseq_aligner}_mapped_read_count.txt
    """

    stub:
    """
    command -v samtools >/dev/null

    echo ${task.ext.stub_read_count ?: 0} > ${meta.id}_${meta.rnaseq_aligner}_mapped_read_count.txt
    """
}
