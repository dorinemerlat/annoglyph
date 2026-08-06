process AGAT_ADD_INTRONS {
    tag "${annotation_method}/${id}"
    // label 'exogap_python'

    input:
    tuple val(id), val(meta), val(annotation_method), path(gff)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_with-introns.gff")

    script:
    """
    module load agat
    agat_sp_add_introns.pl --gff $gff --out ${id}_${annotation_method}_with-introns.gff
    """

    stub:
    """
    touch ${id}_${annotation_method}_with-introns.gff
    """
}