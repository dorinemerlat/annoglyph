process MERGE_TABLES {
    tag "${meta instanceof Map ? meta.id : meta}"

    input:
    tuple val(meta), path(inputs), val(output_name)

    output:
    tuple val(meta), path(output_name)

    script:
    """
    first_file=\$(echo "$inputs" | cut -f 1 -d ' ')
    head -n 1 \$first_file > $output_name
    tail -q -n +2 $inputs >> $output_name
    """

    stub:
    """
    touch $output_name
    """
}
