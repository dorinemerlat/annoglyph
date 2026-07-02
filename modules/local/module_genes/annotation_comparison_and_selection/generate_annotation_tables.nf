process GENERATE_ANNOTATION_TABLES {

    tag "${id}/${annotation_method}/${filtering_status}"
    label 'exogap_python'
    scratch false

    input:
    tuple val(id),
          val(meta),
          val(annotation_method),
          val(filtering_status),
          path(fasta),
          path(gff),
          path(busco_json),
          path(omark_sum)

    output:
    tuple val(id),
          val(meta),
          val(annotation_method),
          val(filtering_status),
          path("${id}_${annotation_method}_${filtering_status}_summary.tsv"),
          path("${id}_${annotation_method}_${filtering_status}_protein_lengths.tsv"),
          path("${id}_${annotation_method}_${filtering_status}_intron_lengths.tsv")

    script:
    """
    generate_annotation_tables.py \\
        --gff ${gff} \\
        --proteins ${fasta} \\
        --busco ${busco_json} \\
        --omark ${omark_sum} \\
        --prefix ${id}_${annotation_method}_${filtering_status}

    test -s ${id}_${annotation_method}_${filtering_status}_summary.tsv
    test -s ${id}_${annotation_method}_${filtering_status}_protein_lengths.tsv
    test -s ${id}_${annotation_method}_${filtering_status}_intron_lengths.tsv
    """

    stub:
    """
    touch ${id}_${annotation_method}_${filtering_status}_summary.tsv
    touch ${id}_${annotation_method}_${filtering_status}_protein_lengths.tsv
    touch ${id}_${annotation_method}_${filtering_status}_intron_lengths.tsv
    """
}