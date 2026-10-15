process UNMASK_GENOME {
    tag "${meta.id}"

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${meta.id}_unmasked.fa")

    script:
    """
    awk '{if (\$0 ~ /^>/) print \$0; else print toupper(\$0)}' $genome > ${meta.id}_unmasked.fa
    """

    stub:
    """
    command -v awk >/dev/null

    cat > ${meta.id}_unmasked.fa << EOF
>${meta.id}
ACGTACGTACGT
EOF
    """
}