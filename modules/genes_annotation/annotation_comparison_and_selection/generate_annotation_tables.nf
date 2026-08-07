process GENERATE_ANNOTATION_TABLES {
    tag "${filtering_status}/${annotation_method}/${id}"
    label 'exogap_python'

    input:
    tuple val(id), val(meta), val(annotation_method), val(filtering_status), path(proteins), path(gff), path(busco), path(omark)

    output:
    tuple val(id), val(meta), val(annotation_method), val(filtering_status),
        path("${id}_${annotation_method}_${filtering_status}_annotation_summary.tsv"),
        path("${id}_${annotation_method}_${filtering_status}_protein_lengths.tsv"),
        path("${id}_${annotation_method}_${filtering_status}_intron_lengths.tsv")

    script:
    def prefix = "${id}_${annotation_method}_${filtering_status}"

    """
    generate_annotation_tables.py \
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
        ${id}_${annotation_method}_${filtering_status}_annotation_summary.tsv \
        ${id}_${annotation_method}_${filtering_status}_protein_lengths.tsv \
        ${id}_${annotation_method}_${filtering_status}_intron_lengths.tsv
    """
}