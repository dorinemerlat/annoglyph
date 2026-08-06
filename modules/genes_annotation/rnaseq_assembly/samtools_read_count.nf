process SAMTOOLS_READ_COUNT {
    tag "${id}"
    label 'varus'

    input:
    tuple val(id), val(meta), val(aligner), path(bam)

    output:
    tuple val(id), val(meta), val(aligner), path("${id}_${aligner}_mapped_read_count.txt")

    script:
    """
    # Count mapped reads.
    samtools view \
        -c \
        -F 4 \
        ${bam} \
        > ${id}_${aligner}_mapped_read_count.txt
    """

    stub:
    """
    command -v samtools >/dev/null

    echo 0 > ${id}_${aligner}_mapped_read_count.txt
    """
}
