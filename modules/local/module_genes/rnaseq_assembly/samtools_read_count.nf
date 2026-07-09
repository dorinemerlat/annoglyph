process SAMTOOLS_READ_COUNT { 
    tag "${aligner}/${id}"
    label 'varus'
    memory { (2 * task.attempt).GB }
    maxRetries 5

    input:
    tuple val(id), val(meta), val(aligner), path(bam)

    output:
    tuple val(id), val(meta), val(aligner), path("${id}_mapped_read_count.txt")

    script:
    """
    # calculate the number of mapped reads 
    samtools view -c -F 4  $bam > ${id}_mapped_read_count.txt
    """

    stub:
    """
    touch ${id}_mapped_read_count.txt
    """
}