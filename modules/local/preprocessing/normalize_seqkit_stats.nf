process NORMALIZE_SEQKIT_STATS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("${meta.id}.stats"), emit: stats

    script:
    """
    awk -v id="${meta.id}" 'BEGIN { FS=OFS="\\t" }
        NR == 1 { \$1 = "genome_id" }
        NR > 1 { \$1 = id }
        { print }
    ' ${stats} > ${meta.id}.stats
    """

    stub:
    """
    printf 'genome_id\\tformat\\ttype\\tnum_seqs\\tsum_len\\n${meta.id}\\tFASTA\\tDNA\\t1\\t4\\n' > ${meta.id}.stats
    """
}
