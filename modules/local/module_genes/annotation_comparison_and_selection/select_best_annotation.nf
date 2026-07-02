process SELECT_BEST_ANNOTATION {

    tag "${status}"
    label 'exogap_python'
    scratch false

    input:
    tuple val(status), path(summary)

    output:
    tuple val(status), path("best_annotations.tsv")

    script:
    """
    select_best_annotations.py \\
        --summary ${summary} \\
        --output best_annotations.tsv

    test -s best_annotations.tsv
    """

    stub:
    """
    cat <<EOF > best_annotations.tsv
genome_id	annotation_method	selected
species_1	braker3	yes
EOF
    """
}