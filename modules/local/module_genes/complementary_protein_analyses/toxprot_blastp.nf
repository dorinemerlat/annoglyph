process TOXPROT_BLASTP {
    tag "${id}"
    scratch false
    cpus 10
    memory '5 GB'

    input:
    tuple val(id), val(meta), path(fasta), path(toxprot_db)

    output:
    tuple val(id), val(meta), path("${id}_toxprot_blastp.tsv")

    script:
    """
    blastp \\
        -query ${fasta} \\
        -db blastdb/toxprot \\
        -out blast.tmp.tsv \\
        -evalue 1e-6 \\
        -max_target_seqs 10 \\
        -num_threads ${task.cpus} \\
        -outfmt '6 qseqid sseqid sacc saccver qlen slen pident length nident mismatch positive ppos gapopen gaps qstart qend sstart send evalue bitscore score qcovhsp qcovs stitle'

    cat <<EOF > ${id}_toxprot_blastp.tsv
    qseqid	sseqid	sacc	saccver	qlen	slen	pident	length	nident	mismatch	positive	ppos	gapopen	gaps	qstart	qend	sstart	send	evalue	bitscore	score	qcovhsp	qcovs	stitle
    EOF

    cat blast.tmp.tsv >> ${id}_toxprot_blastp.tsv

    test -s ${id}_toxprot_blastp.tsv
    """
}