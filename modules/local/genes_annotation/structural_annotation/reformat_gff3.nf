process REFORMAT_GFF3 {
    tag "${meta.annotation_method}/${meta.id}"
    label 'annoglyph_python'

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_reformatted.gff3")

    script:
    def genomeCode = meta.mnemonic.annoglyph
    def outputFile = "${meta.id}_${meta.annotation_method}_reformatted.gff3"

    """
    # Prefix gene and transcript identifiers with the genome mnemonic.
    python3 ${projectDir}/bin/prefix_gff_ids.py --input ${gff} --prefix ${genomeCode} --output ${outputFile}
    """

    stub:
    """
    command -v python3 >/dev/null

    touch ${meta.id}_${meta.annotation_method}_reformatted.gff3
    """
}
