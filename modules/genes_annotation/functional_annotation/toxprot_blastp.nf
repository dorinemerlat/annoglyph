process TOXPROT_BLASTP {
    tag "${id}"
    label 'blast'

    input:
    tuple val(id), val(meta), val(annotation_method), path(fasta), path(toxprot_db)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_toxprot_blastp.tsv")

    script:
    """
    blastp \
        -query ${fasta} \
        -db ${toxprot_db}/toxprot \
        -out blast.tmp.tsv \
        -evalue 1e-6 \
        -max_target_seqs 10 \
        -num_threads ${task.cpus} \
        -outfmt '6 qseqid sseqid sacc saccver qlen slen pident length nident mismatch positive ppos gapopen gaps qstart qend sstart send evalue bitscore score qcovhsp qcovs stitle'

    cat <<EOF > ${id}_${annotation_method}_toxprot_blastp.tsv
qseqid	sseqid	sacc	saccver	qlen	slen	pident	length	nident	mismatch	positive	ppos	gapopen	gaps	qstart	qend	sstart	send	evalue	bitscore	score	qcovhsp	qcovs	stitle
EOF

    cat blast.tmp.tsv >> ${id}_${annotation_method}_toxprot_blastp.tsv
    """

    stub:
    """
    command -v blastp >/dev/null
    touch ${id}_${annotation_method}_toxprot_blastp.tsv
    """
}