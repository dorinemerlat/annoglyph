process INTERPROSCAN {
    tag "${annotation_method}/${id}"
    label 'interproscan'
    cpus 10

    input:
    tuple val(id), val(meta), val(annotation_method), path(proteins)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_ips.xml"),    emit: xml
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_ips.tsv"),    emit: tsv
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_ips.gff3"),   emit: gff

    script:
    """
    interproscan.sh -i ${proteins} -b ${id}_${annotation_method}_ips -f TSV,XML,GFF3 -dp -goterms -iprlookup -pa -t p --cpu ${task.cpus}
    """

    stub:
    """
    touch ${id}_${annotation_method}_ips.xml ${id}_${annotation_method}_ips.tsv ${id}_${annotation_method}_ips.gff3
    """
}