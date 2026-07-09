process SUMMARIZE_REPEATS_BY_LEVEL {
    tag "${id}"
    memory { 10.GB }
    maxRetries 5
    label "exogap_python"

    input:
    tuple val(id), val(meta), path(gff)

    output:
    tuple val(id), val(meta), path("${id}_repeats_summary_by_all.tsv"), emit: all
    tuple val(id), val(meta), path("${id}_repeats_summary_by_all_only_repetitive_elements.tsv"), emit: all_filtered
    tuple val(id), val(meta), path("${id}_repeats_summary_by_feature.tsv"), emit: by_feature
    tuple val(id), val(meta), path("${id}_repeats_summary_*_type.tsv"), emit: by_type
    tuple val(id), val(meta), path("${id}_repeats_summary_*_class.tsv"), emit: by_class
    tuple val(id), val(meta), path("${id}_repeats_summary_*_order.tsv"), emit: by_order
    tuple val(id), val(meta), path("${id}_repeats_summary_*_superfamily.tsv"), emit: by_superfamily
    tuple val(id), val(meta), path("${id}_repeats_summary_*_family.tsv"), emit: by_family
    tuple val(id), val(meta), path("${id}_repeats_summary_*_subfamily.tsv"), emit: by_subfamily
    tuple val(id), val(meta), path("${id}_repeats_summary_*_subfamily2.tsv"), emit: by_subfamily2

    script:
    """
    summarize_repeats_by_level.py -i $gff -s ${meta.assembly_size}
    """

    stub:
    """
    touch ${id}_RE_summary.tsv
    """
}
