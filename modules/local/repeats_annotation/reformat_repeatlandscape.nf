process REFORMAT_REPEATLANDSCAPE {
    tag "${meta.id}"
    label "exogap_python"

    input:
    tuple val(meta), path(tsv), path(classification)

    output:
    tuple val(meta), path("${meta.id}_TE_repeat_landscape_reformated.tsv")

    script:
    """
    reformat_repeatlandscape.py \
        -i ${tsv} \
        -c ${classification} \
        --id ${meta.id} \
        -o ${meta.id}_TE_repeat_landscape_reformated.tsv
    """

    stub:
    """
    touch ${meta.id}_TE_repeat_landscape_reformated.tsv
    """
}