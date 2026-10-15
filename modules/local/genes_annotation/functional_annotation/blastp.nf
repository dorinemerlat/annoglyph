process BLASTP {
    tag "${meta.annotation_method}/${meta.id}"
    label 'blast'

    label 'program_blastp'
    input:
    tuple val(meta), path(proteins), val(db), val(evalue), val(word_size), val(max_hsps)

    output:
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_blastp.asn"), emit: archive
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_blastp.tsv"), emit: tsv

    script:
    def prefix = "${meta.id}_${meta.annotation_method}_blastp"

    """
    blastp \
        -query ${proteins} \
        -db ${db} \
        -evalue ${evalue} \
        -word_size ${word_size} \
        -max_hsps ${max_hsps} \
        -outfmt 11 \
        -out ${prefix}.asn \
        -num_threads ${task.cpus} ${params.program_args.blastp ?: ''}

    blast_formatter \
        -archive ${prefix}.asn \
        -outfmt "6 qseqid sseqid sacc saccver qlen slen pident length nident mismatch positive ppos gapopen gaps qstart qend sstart send evalue bitscore score qcovhsp qcovs stitle" \
        -out ${prefix}.raw.tsv

    # The annotation-transfer parser consumes a TSV with named columns.
    cat > ${prefix}.tsv <<'HEADER'
qseqid	sseqid	sacc	saccver	qlen	slen	pident	length	nident	mismatch	positive	ppos	gapopen	gaps	qstart	qend	sstart	send	evalue	bitscore	score	qcovhsp	qcovs	stitle
HEADER
    cat ${prefix}.raw.tsv >> ${prefix}.tsv
    rm ${prefix}.raw.tsv
    """

    stub:
    def prefix = "${meta.id}_${meta.annotation_method}_blastp"
    """
    command -v blastp >/dev/null
    command -v blast_formatter >/dev/null

    touch ${prefix}.asn ${prefix}.tsv
    """
}
