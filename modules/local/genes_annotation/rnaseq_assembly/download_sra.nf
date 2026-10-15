process DOWNLOAD_SRA {
    tag "${sra_id}"
    label 'rna'

    input:
    tuple val(meta), val(sra_id)

    output:
    tuple val(meta), path("${sra_id}*.fastq.gz")

    script:
    """
    # Randomize requests to reduce simultaneous access to the SRA servers.
    sleep \$((RANDOM % 121))

    prefetch \
        --max-size 100G \
        ${sra_id}

    fasterq-dump \
        --split-files \
        --threads ${task.cpus} \
        ${sra_id}

    pigz \
        --processes ${task.cpus} \
        *.fastq
    """

    stub:
    """
    command -v prefetch >/dev/null
    command -v fasterq-dump >/dev/null
    command -v pigz >/dev/null

    touch ${sra_id}.fastq.gz
    """
}