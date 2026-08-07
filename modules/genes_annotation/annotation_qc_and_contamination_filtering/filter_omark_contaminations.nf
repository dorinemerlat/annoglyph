process FILTER_OMARK_CONTAMINATIONS {
    tag "${annotation_method}/${id}"
    label 'exogap_python'

    input:
    tuple val(id), val(meta), val(annotation_method), path(gff), path(fasta), path(ump), val(threshold)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_contaminations.txt"), emit: contaminations
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_contaminations.fa"), emit: fasta
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_contaminations.log"), emit: log

    script:
    def gffForOmark = "${id}_${annotation_method}_for_omark.gff"
    def outputPrefix = "${id}_${annotation_method}_contaminations"

    """
    reformat_gff_for_omark.py \
        -i ${gff} \
        -o ${gffForOmark}

    # The OMArk .ump file is staged in the current directory and detected with -i .
    contamination_chromosome_filtering.py \
        -i . \
        -g ${gffForOmark} \
        -f ${fasta} \
        -o ${outputPrefix} \
        -t ${threshold} \
        > ${outputPrefix}.log 2>&1
    """

    stub:
    """
    touch \
        ${id}_${annotation_method}_contaminations.txt \
        ${id}_${annotation_method}_contaminations.fa \
        ${id}_${annotation_method}_contaminations.log
    """
}