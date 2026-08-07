process OMARK {
    tag "${annotation_method}/${id}"
    label 'omark'

    input:
    tuple val(id), val(meta), val(annotation_method), path(proteins), path(luca_db)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_omark_detailed_summary.txt"), emit: summary
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_omark.sum"), emit: sum
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_omark.ump"), emit: ump
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_contaminants.count"), emit: contaminants_count

    script:
    def prefix = "${id}_${annotation_method}_omark"

    """
    omamer search \
        --db ${luca_db} \
        --query ${proteins} \
        --out ${prefix}.omamer

    mkdir omark

    omark \
        --file ${prefix}.omamer \
        --database ${luca_db} \
        --outputFolder omark \
        --taxid ${meta.taxid}

    mv omark/* .

    grep '^A:' ${prefix}.sum \
        | sed -E 's/.*C:([0-9]+).*/\\1/' \
        > ${id}_${annotation_method}_contaminants.count
    """

    stub:
    """
    touch \
        ${id}_${annotation_method}_omark_detailed_summary.txt \
        ${id}_${annotation_method}_omark.sum \
        ${id}_${annotation_method}_omark.ump \
        ${id}_${annotation_method}_contaminants.count
    """
}