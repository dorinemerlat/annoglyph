process AGAT_STATISTICS {

    tag "${id}/${annotation_method}"
    label 'exogap_agat'
    scratch false

    input:
    tuple val(id),
          val(meta),
          val(annotation_method),
          path(gff)

    output:
    tuple val(id),
          val(meta),
          val(annotation_method),
          path("${id}_${annotation_method}_summary.tsv"),
          path("${id}_${annotation_method}_protein_lengths.tsv"),
          path("${id}_${annotation_method}_intron_lengths.tsv"),
          emit: stats

    script:
    """
    agat_sp_statistics.pl \\
        --gff ${gff} \\
        --output ${id}_${annotation_method}

    mv ${id}_${annotation_method}_summary.txt ${id}_${annotation_method}_summary.tsv
    mv ${id}_${annotation_method}_protein_length_distribution.txt ${id}_${annotation_method}_protein_lengths.tsv
    mv ${id}_${annotation_method}_intron_length_distribution.txt ${id}_${annotation_method}_intron_lengths.tsv

    test -s ${id}_${annotation_method}_summary.tsv
    test -s ${id}_${annotation_method}_protein_lengths.tsv
    test -s ${id}_${annotation_method}_intron_lengths.tsv
    """

    stub:
    """
    touch ${id}_${annotation_method}_summary.tsv
    touch ${id}_${annotation_method}_protein_lengths.tsv
    touch ${id}_${annotation_method}_intron_lengths.tsv
    """
}