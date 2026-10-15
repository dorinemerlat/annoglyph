process CLASSIFY_GENES_VS_REPEATS {
    tag "${annotation_method}/${id}"
    label 'exogap_python'
    scratch false
    stageInMode 'copy'

    input:
    tuple val(id), val(meta), val(annotation_method), path(proteins), path(repeats)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_genes_vs_merged_repeats.tsv"), emit: merge
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_genes_vs_no_merged_repeats.tsv"), emit: no_merge

    script:
    """
    classify_genes_vs_repeats.py \
        --proteins ${proteins} \
        --repeats ${repeats} \
        --output ${id}_${annotation_method}_genes_vs_merged_repeats.tsv \
        --merge-repeats \
        --tmp tmp

    classify_genes_vs_repeats.py \
        --proteins ${proteins} \
        --repeats ${repeats} \
        --output ${id}_${annotation_method}_genes_vs_no_merged_repeats.tsv \
        --tmp tmp
    """

    stub:
    """
    touch ${id}_${annotation_method}_genes_vs_merged_repeats.tsv ${id}_${annotation_method}_genes_vs_no_merged_repeats.tsv 
    """
} 