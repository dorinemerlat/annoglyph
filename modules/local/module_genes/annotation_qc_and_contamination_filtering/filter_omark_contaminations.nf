process FILTER_OMARK_CONTAMINATIONS {
    tag "${annotation_method}/c-${threshold}/${id}"
    label 'exogap_python'
    memory '20 GB'
    scratch false
    
    input:
    tuple val(id), val(meta), val(annotation_method), path(gff), path(fasta), path(ump), val(threshold)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_contaminations-c-${threshold}.txt"), emit: contaminations
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_contaminations-filtered-c-${threshold}.fa"), emit: fasta

    script:
    def gff_for_omark = "${id}_${annotation_method}_for_omark.gff"
    def contaminations = "${id}_${annotation_method}_contaminations-c-${threshold}"
    def contaminations_filtered = "${id}_${annotation_method}_contaminations-filtered-c-${threshold}"
    """
    echo \$PATH
    ls /shared/ifbstor1/projects/metainvert/exogap/bin
    reformat_gff_for_omark.py -i $gff -o $gff_for_omark
    contamination_chromosome_filtering.py -i "." -g $gff_for_omark -f $fasta -o $contaminations -t $threshold
    mv ${contaminations}.fa ${contaminations_filtered}.fa
    """

    stub:
    """
    touch ${id}_${annotation_method}_contaminations-c-${threshold}.txt ${id}_${annotation_method}_contaminations-filtered-c-${threshold}.fa
    """
}