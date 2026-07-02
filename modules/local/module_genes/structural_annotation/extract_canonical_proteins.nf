process EXTRACT_CANONICAL_PROTEINS {
    tag "${annotation_method}/${id}"
    // label 'agat'
    scratch false
    stageInMode 'copy'
    label 'retry_with_backoff'
    maxRetries 5
    
    input:
    tuple val(id), val(meta), val(annotation_method), path(genome), path(gff)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_canonical.gff"), emit: gff
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_proteins.fasta"), emit: fasta
    script:
    """
    # extract the longest isoform for each gene using AGAT
    module load agat
    agat_sp_keep_longest_isoform.pl --gff $gff -o ${id}_${annotation_method}_canonical.gff

    # extract proteins sequences
    agat_sp_extract_sequences.pl -g ${id}_${annotation_method}_canonical.gff -f $genome -o ${id}_${annotation_method}_proteins.fasta --protein
    """

    stub:
    """
    touch ${id}_${annotation_method}_canonical.gff ${id}_${annotation_method}_proteins.fasta
    """
}