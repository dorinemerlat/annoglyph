process FILTER_OMARK_CONTAMINATIONS {
    tag "${annotation_method}/${id}"
    label 'exogap_python'
    memory '20 GB'
    scratch false
    
    input:
    tuple val(id), val(meta), val(annotation_method), path(gff), path(fasta), path(ump), val(threshold)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_contaminations.txt"), emit: contaminations
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_contaminations.fa"), emit: fasta
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_contaminations.log"), emit: log

    script:
    def gff_for_omark = "${id}_${annotation_method}_for_omark.gff"
    def contaminations = "${id}_${annotation_method}_contaminations"
    """
    reformat_gff_for_omark.py -i $gff -o $gff_for_omark
    contamination_chromosome_filtering.py -i "." -g $gff_for_omark -f $fasta -o $contaminations -t $threshold > ${contaminations}.log
    """

    stub:
    """
    touch ${id}_${annotation_method}_contaminations.txt ${id}_${annotation_method}_contaminations.fa ${id}_${annotation_method}_contaminations.log
    """
}
