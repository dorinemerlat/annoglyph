process REFORMAT_GFF3 {
    tag "${meta.annotation_method}/${meta.id}"
    label 'annoglyph_tools'

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_reformatted.gff3")

    script:
    def genomeCode = meta.mnemonic.annoglyph
    def outputFile = "${meta.id}_${meta.annotation_method}_reformatted.gff3"

    """
    # Prefix gene and transcript identifiers with the genome mnemonic.
    sed -E 's/(^|[;,])(ID|Parent)=g/\\1\\2=${genomeCode}.g/g' \
        ${gff} \
        > ${outputFile}
    """

    stub:
    """
    command -v sed >/dev/null

    touch ${meta.id}_${meta.annotation_method}_reformatted.gff3
    """
}
