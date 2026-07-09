process SEQKIT_STATS {
    tag "${id}"
    label "seqkit"

    input:
    tuple val(id), val(meta), val(genome)
    each seqtype

    output:
    tuple val(id), val(meta), path("${id}.stats")

    script:
    """
    # compute genome stats
    seqkit stats --tabular --seq-type $seqtype --all $genome \
        | awk -v id="${id}" 'BEGIN { FS=OFS="\\t" }
            NR == 1 { \$1 = "genome_id" }
            NR == 2 { \$1 = id }
            { print }
        ' > ${id}.stats
    """

    stub:
    """
    touch ${id}.stats
    """
}
