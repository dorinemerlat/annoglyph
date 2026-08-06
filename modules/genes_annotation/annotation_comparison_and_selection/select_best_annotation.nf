process SELECT_BEST_ANNOTATION {
    tag ""
    label 'exogap_python'
    scratch false

    input:
    tuple val(status), path(summary)

    output:
    tuple val(status), path("best_annotations.tsv")

    script:
    """
    select_best_annotations.py -i ${summary} -o best_annotations.tsv
    """

    stub:
    """
    touch best_annotations.tsv 
    """
}