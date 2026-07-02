process BLASTP {
    tag "${annotation_method}/${id}"
    label 'blast'
    cpus 10

    input:
    tuple val(id), val(meta), val(annotation_method), path(proteins), val(db), val(evalue), val(word_size), val(max_hsps)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_blastp.asn"), emit: archive
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_blastp.xml"), emit: xml
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_blastp.tsv"), emit: tsv

    script:
    def archive = "${id}_${annotation_method}_blastp.asn"
    """
    blastp -query $proteins -db $db -evalue $evalue -word_size $word_size -max_hsps $max_hsps -outfmt 11 -out ${archive} -num_threads ${task.cpus}
    blast_formatter -archive $archive -outfmt 5 -out ${id}_${annotation_method}_blastp.xml
    blast_formatter -archive $archive -outfmt 6 -out ${id}_${annotation_method}_blastp.tsv
    """


    stub:
    """
    touch ${id}_${annotation_method}_blastp.asn ${id}_${annotation_method}_blastp.xml ${id}_${annotation_method}_blastp.tsv
    """
}