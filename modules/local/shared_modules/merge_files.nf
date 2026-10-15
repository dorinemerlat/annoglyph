process MERGE_FILES {
    tag "${meta instanceof Map ? meta.id : meta}"

    input:
    tuple val(meta), path(inputs), val(output_name)

    output:
    tuple val(meta), path(output_name)

    script:
    """
    cat $inputs > $output_name
    """

    stub:
    """
    touch $output_name
    """
}
