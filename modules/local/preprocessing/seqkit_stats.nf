process SEQKIT_STATS {
    tag "${meta.id}"
    label "seqkit"

    input:
    tuple val(meta), path(genome)
    each seqtype

    output:
    tuple val(meta), path("${meta.id}.stats")

    script:
    """
    # Compute sequence statistics
    seqkit stats --tabular --seq-type $seqtype --all $genome \
        | awk -v id="${meta.id}" 'BEGIN { FS=OFS="\\t" }
            NR == 1 { \$1 = "genome_id" }
            NR == 2 { \$1 = id }
            { print }
        ' > ${meta.id}.stats
    """

    stub:
    """
    # Check that seqkit is available
    command -v seqkit >/dev/null

    cat > ${meta.id}.stats << EOF
genome_id	format	type	num_seqs	sum_len	min_len	avg_len	max_len	Q1	Q2	Q3	sum_gap	N50	Q20(%)	Q30(%)	GC(%)
${meta.id}	FASTA	DNA	0	0	0	0	0	0	0	0	0	0	0	0	0
EOF
    """
}
