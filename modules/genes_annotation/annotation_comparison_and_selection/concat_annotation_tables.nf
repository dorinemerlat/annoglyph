process CONCAT_ANNOTATION_TABLES {
    tag "${filtering_status}"

    input:
    tuple val(filtering_status), path(summary), path(protein_lengths), path(intron_lengths)

    output:
    tuple val(filtering_status),
        path("${filtering_status}_summary.tsv"),
        path("${filtering_status}_protein_lengths.tsv"),
        path("${filtering_status}_intron_lengths.tsv")

    script:
    """
    concat_tables.sh \
        ${filtering_status}_summary.tsv \
        ${summary}

    concat_tables.sh \
        ${filtering_status}_protein_lengths.tsv \
        ${protein_lengths}

    concat_tables.sh \
        ${filtering_status}_intron_lengths.tsv \
        ${intron_lengths}
    """

    stub:
    """
    command -v concat_tables.sh >/dev/null

    touch \
        ${filtering_status}_summary.tsv \
        ${filtering_status}_protein_lengths.tsv \
        ${filtering_status}_intron_lengths.tsv
    """
}