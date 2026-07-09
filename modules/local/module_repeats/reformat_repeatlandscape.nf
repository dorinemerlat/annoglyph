process REFORMAT_REPEATLANDSCAPE {
    tag "${id}"
    label "exogap_python"

    input:
    tuple val(id), val(meta), path(tsv), path(classification)

    output:
    tuple val(id), val(meta), path("${id}_TE_repeat_landscape_reformated.tsv")

    script:
    """
    reformat_repeatlandscape.py -i $tsv -c $classification --id $id -o ${id}_TE_repeat_landscape_reformated.tsv
    """

    stub:
    """
    touch ${id}_TE_repeat_landscape_reformated.tsv
    """
}
