process REPEATMASKER_OUT_TO_GFF {
    tag "${meta.id}"
    label 'annoglyph_python'
    
    input:
    tuple val(meta), path(out), path(table), path(classification_csv)

    output:
    tuple val(meta), path("${meta.id}_repeats.gff")

    script:
    """
    repeatmasker_out_to_gff.py -r ${out} -m ${table} -o ${meta.id}_repeats.gff -c ${classification_csv}
    """

    stub:
    """
    touch ${meta.id}_repeats.gff
    """
}
