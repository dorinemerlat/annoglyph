process MERGE_BAM {
    tag "${id}"
    time '2d'
    label 'varus'
    memory { (10 * task.attempt).GB }
    maxRetries 5
    cpus 4

    input:
    tuple val(id), val(meta), path(bams, stageAs: 'bams/*')

    output:
    tuple val(id), val(meta), path("${id}.bam")

    script:
    """
    samtools merge --threads ${task.cpus} -o ${id}.bam bams/*.bam

    # calculate the number of mapped reads 
    samtools view -c -F 4  ${id}.bam > ${id}_mapped_read_count.txt
    """

    stub:
    """
    touch ${id}.bam
    """
}