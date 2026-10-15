process FILTER_OMARK_CONTAMINATIONS {
    tag "${meta.annotation_method}/${meta.id}"
    label 'annoglyph_python'
    label 'omark_filter'

    input:
    tuple val(meta), path(gff), path(fasta), path(ump), val(threshold)

    output:
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_contaminations.txt"), emit: contaminations
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_contaminations.fa"), emit: fasta
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_contaminations.log"), emit: log

    script:
    def gffForOmark = "${meta.id}_${meta.annotation_method}_for_omark.gff"
    def outputPrefix = "${meta.id}_${meta.annotation_method}_contaminations"

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
        ${meta.id}_${meta.annotation_method}_contaminations.txt \
        ${meta.id}_${meta.annotation_method}_contaminations.fa \
        ${meta.id}_${meta.annotation_method}_contaminations.log
    """
}