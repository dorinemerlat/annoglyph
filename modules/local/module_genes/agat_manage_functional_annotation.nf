process AGAT_MANAGE_FUNCTIONAL_ANNOTATION {
    tag "${annotation_method}/${id}"
    // label 'exogap_python'
    cpus 1
    time '6h'
    stageInMode 'copy'
    scratch false

    input:
    tuple val(id), val(meta), val(annotation_method), path(gff), path(fasta), path(blastp), path(inteproscan)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_with-functional-annotation.gff")

    script:
    """
    module load agat
    # reformat blastp output for agat
    awk -F '\t' 'BEGIN{OFS="\t"} NF>=20 {print \$1,\$2,\$7,\$8,\$10,\$13,\$15,\$16,\$17,\$18,\$19,\$20}' \
    ${blastp} \
    > ${id}_${annotation_method}_blastp_for_agat.tsv

    # add functional annotation to gff
    agat_sp_manage_functional_annotation.pl \
        --gff ${gff} \
        --db ${fasta} \
        --blast ${id}_${annotation_method}_blastp_for_agat.tsv \
        --interpro ${inteproscan} \
        --output ${id}_${annotation_method}_with-functional-annotation.gff
    """

    stub:
    """
    touch ${id}_${annotation_method}_with-functional-annotation.gff
    """
}