process SUMMARIZE_REPEATS_BY_LEVEL {
    tag "${meta.id}"
    label "annoglyph_python"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("${meta.id}_repeats_summary_by_all.tsv"), emit: all
    tuple val(meta), path("${meta.id}_repeats_summary_by_all_only_repetitive_elements.tsv"), emit: all_filtered
    tuple val(meta), path("${meta.id}_repeats_summary_by_feature.tsv"), emit: by_feature
    tuple val(meta), path("${meta.id}_repeats_summary_*_type.tsv"), emit: by_type
    tuple val(meta), path("${meta.id}_repeats_summary_*_class.tsv"), emit: by_class
    tuple val(meta), path("${meta.id}_repeats_summary_*_order.tsv"), emit: by_order
    tuple val(meta), path("${meta.id}_repeats_summary_*_superfamily.tsv"), emit: by_superfamily
    tuple val(meta), path("${meta.id}_repeats_summary_*_family.tsv"), emit: by_family
    tuple val(meta), path("${meta.id}_repeats_summary_*_subfamily.tsv"), emit: by_subfamily
    tuple val(meta), path("${meta.id}_repeats_summary_*_subfamily2.tsv"), emit: by_subfamily2

    script:
    """
    summarize_repeats_by_level.py -i ${gff} -s ${meta.assembly_size}
    """

    stub:
    """
    touch ${meta.id}_repeats_summary_by_all.tsv
    touch ${meta.id}_repeats_summary_by_all_only_repetitive_elements.tsv
    touch ${meta.id}_repeats_summary_by_feature.tsv
    touch ${meta.id}_repeats_summary_stub_type.tsv
    touch ${meta.id}_repeats_summary_stub_class.tsv
    touch ${meta.id}_repeats_summary_stub_order.tsv
    touch ${meta.id}_repeats_summary_stub_superfamily.tsv
    touch ${meta.id}_repeats_summary_stub_family.tsv
    touch ${meta.id}_repeats_summary_stub_subfamily.tsv
    touch ${meta.id}_repeats_summary_stub_subfamily2.tsv
    """
}