process SORT_RNASEQ_BAM {
    tag "${meta.genome_id}/${meta.run_id}/${meta.readset_signature}"
    cpus 4
    maxForks 4
    scratch false
    stageInMode 'symlink'
    container 'docker://dorinemerlat/varus:v1.1'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.sorted.bam'), emit: bam
    path 'versions.yml', emit: versions

    script:
    // samtools -m is per sorting thread. Keep at least 20% for overhead.
    def threads = Math.max(1, task.cpus as int)
    def availableMb = Math.floor(task.memory.toMega() * 0.8 / threads) as long
    def perThreadMb = Math.min((params.rnaseq_sort_memory_per_thread_mb as long) * task.attempt, availableMb)
    if (perThreadMb < 1) error 'Insufficient memory for samtools sort'
    """
    mkdir -p sort_tmp
    samtools sort -@ ${threads - 1} -m ${perThreadMb}M \
        -T sort_tmp/${meta.id} -o '${meta.id}.sorted.bam' '${bam}'
    samtools quickcheck '${meta.id}.sorted.bam'
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        samtools: \$(samtools --version | sed -n '1s/samtools //p')
    END_VERSIONS
    """

    stub:
    """
    touch '${meta.id}.sorted.bam'
    echo 'samtools: stub' > versions.yml
    """
}
