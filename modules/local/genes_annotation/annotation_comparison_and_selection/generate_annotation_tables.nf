process GENERATE_ANNOTATION_TABLES {
    tag "${meta.filtering_status}/${meta.annotation_method}/${meta.id}"
    label 'annoglyph_python'

    input:
    tuple val(meta), path(proteins), path(gff), path(busco), path(omark)

    output:
    tuple val(meta),
        path("${meta.id}_${meta.annotation_method}_${meta.filtering_status}_annotation_summary.tsv"),
        path("${meta.id}_${meta.annotation_method}_${meta.filtering_status}_protein_lengths.tsv"),
        path("${meta.id}_${meta.annotation_method}_${meta.filtering_status}_intron_lengths.tsv")

    script:
    def prefix = "${meta.id}_${meta.annotation_method}_${meta.filtering_status}"

    """
    generate_annotation_tables.py \
        --genome-id "${meta.id}" --annotation-method "${meta.annotation_method}" \
        --proteins ${proteins} \
        --gff ${gff} \
        --busco ${busco} \
        --omark ${omark} \
        --out ${prefix}
    """

    stub:
    """
    command -v generate_annotation_tables.py >/dev/null

    touch \
        ${meta.id}_${meta.annotation_method}_${meta.filtering_status}_annotation_summary.tsv \
        ${meta.id}_${meta.annotation_method}_${meta.filtering_status}_protein_lengths.tsv \
        ${meta.id}_${meta.annotation_method}_${meta.filtering_status}_intron_lengths.tsv
    """
}