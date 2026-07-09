process REFORMAT_GFF3 {
    tag "${annotation_method}/${id}"
    stageInMode 'copy'

    input:
    tuple val(id), val(meta), val(annotation_method), path(gff)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_reformatted.gff3")

    script:
    def genome_code = meta.mnemonic.exogap
    """
    sed 's/=g/=${genome_code}.g/g' "${gff}" > "${id}_${annotation_method}_reformatted.gff3"
    """

    stub:
    """
    touch "${id}_${annotation_method}_reformatted.gff3"
    """
}