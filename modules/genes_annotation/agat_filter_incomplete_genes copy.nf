process AGAT_START_AND_STOP_CODON {
    tag "${annotation_method}/${id}"
    // label 'exogap_python'
    
    input:
    tuple val(id), val(meta), val(annotation_method), path(genome), path(gff)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_start_stop.gff")

    script:
    """
    module load agat
    agat_sp_add_start_and_stop.pl --gff $gff --fasta $genome -o ${id}_${annotation_method}_start_stop.gff
    """

    stub:
    """
    touch ${id}_${annotation_method}_start_stop.gff
    """
}