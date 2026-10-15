process INTERPROSCAN {
    tag "${meta.annotation_method}/${meta.id}"
    label 'interproscan'

    input:
    tuple val(meta), path(proteins)

    output:
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_ips.tsv"), emit: tsv

    script:
    """
    interproscan.sh \
        -i ${proteins} \
        -b ${meta.id}_${meta.annotation_method}_ips \
        -f TSV \
        -dp \
        -goterms \
        -iprlookup \
        -pa \
        -t p \
        --cpu ${task.cpus}
    """

    stub:
    """
    command -v interproscan.sh >/dev/null

    touch ${meta.id}_${meta.annotation_method}_ips.tsv
    """
}
