process EGGNOG_MAPPER {
    tag "${annotation_method}/${id}"
    label 'eggnog'

    input:
    tuple val(id), val(meta), val(annotation_method), path(proteins)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${prefix}/${prefix}.emapper.annotations"), emit: annotations
    tuple val(id), val(meta), val(annotation_method), path("${prefix}/${prefix}.emapper.seed_orthologs"), emit: seed_orthologs
    tuple val(id), val(meta), val(annotation_method), path("${prefix}/${prefix}.emapper.hits"), emit: hits
    tuple val(id), val(meta), val(annotation_method), path("${prefix}/${prefix}.emapper.orthologs"), emit: orthologs

    script:
    def prefix = "${id}_${annotation_method}_eggnog"

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
        --data_dir /shared/bank/emapperdb/5.0.2/
    """

    stub:
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