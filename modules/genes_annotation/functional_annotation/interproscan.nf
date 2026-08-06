process INTERPROSCAN {
    tag "${annotation_method}/${id}"
    label 'interproscan'

    input:
    tuple val(id), val(meta), val(annotation_method), path(proteins)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_ips.tsv"), emit: tsv

    script:
    """
    interproscan.sh \
        -i ${proteins} \
        -b ${id}_${annotation_method}_ips \
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

    touch ${id}_${annotation_method}_ips.tsv
    """
}
