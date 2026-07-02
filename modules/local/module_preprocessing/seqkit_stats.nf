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
    seqkit stats --tabular --seq-type $seqtype --all $genome > ${id}.stats
    """

    stub:
    """
    touch ${id}.stats
    """
}