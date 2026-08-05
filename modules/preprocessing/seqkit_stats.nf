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
    # Compute sequence statistics
    seqkit stats --tabular --seq-type $seqtype --all $genome \
        | awk -v id="${id}" 'BEGIN { FS=OFS="\\t" }
            NR == 1 { \$1 = "genome_id" }
            NR == 2 { \$1 = id }
            { print }
        ' > ${id}.stats
    """

    stub:
    """
    # Check that seqkit is available
    command -v seqkit >/dev/null

    cat > ${id}.stats << EOF
genome_id	format	type	num_seqs	sum_len	min_len	avg_len	max_len	Q1	Q2	Q3	sum_gap	N50	Q20(%)	Q30(%)	GC(%)
${id}	FASTA	DNA	0	0	0	0	0	0	0	0	0	0	0	0	0
EOF
    """
}
