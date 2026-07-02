process CONCAT_ANNOTATION_TABLES {

    tag "${filtering_status}"

    input:
    tuple val(filtering_status),
          path(summary),
          path(protein_lengths),
          path(intron_lengths)

    output:
    tuple val(filtering_status),
          path("${filtering_status}.summary.tsv"),
          path("${filtering_status}.protein_lengths.tsv"),
          path("${filtering_status}.intron_lengths.tsv")

    script:
    """
    concat_table.sh \
        ${filtering_status}.summary.tsv \
        ${summary}

    concat_table.sh \
        ${filtering_status}.protein_lengths.tsv \
        ${protein_lengths}

    concat_table.sh \
        ${filtering_status}.intron_lengths.tsv \
        ${intron_lengths}
    """
}