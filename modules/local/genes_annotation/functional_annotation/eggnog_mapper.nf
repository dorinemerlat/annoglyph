process EGGNOG_MAPPER {
    tag "${meta.annotation_method}/${meta.id}"
    label 'eggnog'

    input:
    tuple val(meta), path(proteins)

    output:
    tuple val(meta), path("${prefix}/${prefix}.emapper.annotations"), emit: annotations
    tuple val(meta), path("${prefix}/${prefix}.emapper.seed_orthologs"), emit: seed_orthologs
    tuple val(meta), path("${prefix}/${prefix}.emapper.hits"), emit: hits
    tuple val(meta), path("${prefix}/${prefix}.emapper.orthologs"), emit: orthologs

    script:
    def prefix = "${meta.id}_${meta.annotation_method}_eggnog"

    """
    mkdir -p ${prefix}

    emapper.py \
        -i ${proteins} \
        --cpu ${task.cpus} \
        --report_orthologs \
        --output ${prefix} \
        --output_dir ${prefix} \
        --scratch_dir /tmp \
        --temp_dir /tmp \
        --data_dir ${params.eggnog_db}
    """

    stub:
    def prefix = "${meta.id}_${meta.annotation_method}_eggnog"
    """
    command -v emapper.py >/dev/null

    mkdir -p ${prefix}

    touch \
        ${prefix}/${prefix}.emapper.annotations \
        ${prefix}/${prefix}.emapper.seed_orthologs \
        ${prefix}/${prefix}.emapper.hits \
        ${prefix}/${prefix}.emapper.orthologs
    """
}