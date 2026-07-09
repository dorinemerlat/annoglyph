process UNMASK_GENOME {
    tag "${id}"

    input:
    tuple  val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path("${id}_unmasked.fa")

    script:
    """
    awk '{if (\$0 ~ /^>/) print \$0; else print toupper(\$0)}' $genome > ${id}_unmasked.fa
    """

    stub:
    """
    touch ${id}_unmasked.fa
    """
}
