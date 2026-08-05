process MERGE_TABLES {
    tag "${id}"

    input:
    tuple val(id), path(inputs), val(output_name)

    output:
    tuple val(id), path(output_name)

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
