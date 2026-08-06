process BLASTP {
    tag "${annotation_method}/${id}"
    label 'blast'

    input:
    tuple val(id), val(meta), val(annotation_method), path(proteins), val(db), val(evalue), val(word_size), val(max_hsps)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_blastp.asn"), emit: archive
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_blastp.tsv"), emit: tsv

    script:
    def prefix = "${id}_${annotation_method}_blastp"

    """
    blastp \
        -query ${proteins} \
        -db ${db} \
        -evalue ${evalue} \
        -word_size ${word_size} \
        -max_hsps ${max_hsps} \
        -outfmt 11 \
        -out ${prefix}.asn \
        -num_threads ${task.cpus}

    blast_formatter \
        -archive ${prefix}.asn \
        -outfmt "6 qseqid sseqid sacc saccver qlen slen pident length nident mismatch positive ppos gapopen gaps qstart qend sstart send evalue bitscore score qcovhsp qcovs stitle" \
        -out ${prefix}.tsv
    """

    stub:
    """
    command -v blastp >/dev/null
    command -v blast_formatter >/dev/null

    touch ${prefix}.asn ${prefix}.tsv
    """
}
