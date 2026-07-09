process DOWNLOAD_SRA {
    tag "${id}"
    scratch false
    label 'retry_with_backoff'
    memory {"${10 + (2 * (task.attempt - 1))} GB"}
    maxRetries 3
    cpus 10

    input:
    tuple val(id), val(sra_id)

    output:
    tuple val(id), path("${sra_id}*.fastq.gz")

    script:
    """
    module load sra-tools
    sleep \$((RANDOM % 121))
    prefetch ${sra_id}
    fasterq-dump --split-files -e 10 ${sra_id}
    pigz *.fastq
    """

    stub:
    """
    touch ${sra_id}_1.fastq.gz ${sra_id}_2.fastq.gz
    """
}
