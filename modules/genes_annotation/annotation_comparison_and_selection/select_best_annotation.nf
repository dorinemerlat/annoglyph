process SELECT_BEST_ANNOTATION {
    label 'exogap_python'

    input:
    tuple val(status), path(summary)

    output:
    tuple val(status), path("best_annotations.tsv")

    script:
    """
    select_best_annotations.py \
        -i ${summary} \
        -o best_annotations.tsv
    """

    stub:
    """
    command -v select_best_annotations.py >/dev/null

    touch best_annotations.tsv 
    """
}