process MERGE_FILES {
    tag "${id}"

    input:
    tuple val(id), path(inputs), val(output_name)

    output:
    tuple val(id), path(output_name)

    script:
    """
    cat $inputs > $output_name
    """

    stub:
    """
    touch $output_name
    """
}
