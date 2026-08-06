process REFORMAT_GFF3 {
    tag "${annotation_method}/${id}"
    label 'exogap_tools'

    input:
    tuple val(id), val(meta), val(annotation_method), path(gff)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_reformatted.gff3")

    script:
    def genomeCode = meta.mnemonic.exogap
    def outputFile = "${id}_${annotation_method}_reformatted.gff3"

    """
    # Prefix gene and transcript identifiers with the genome mnemonic.
    sed -E 's/(^|[;,])(ID|Parent)=g/\\1\\2=${genomeCode}.g/g' \
        ${gff} \
        > ${outputFile}
    """

    stub:
    """
    command -v sed >/dev/null

    touch ${id}_${annotation_method}_reformatted.gff3
    """
}
