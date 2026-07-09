process REPEATMASKER_OUT_TO_GFF {
    tag "${id}"
    label 'exogap_python'
    
    input:
    tuple val(id), val(meta), path(out), path(table), path(classification_csv)

    output:
    tuple val(id), val(meta), path("${id}_repeats.gff")

    script:
    """
    repeatmasker_out_to_gff.py -r ${out} -m ${table} -o ${id}_repeats.gff -c ${classification_csv}
    """

    stub:
    """
    touch ${id}_repeats.gff
    """
}
