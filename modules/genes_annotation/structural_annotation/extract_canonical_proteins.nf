process EXTRACT_CANONICAL_PROTEINS {
    tag "${annotation_method}/${id}"
    label 'agat'

    input:
    tuple val(id), val(meta), val(annotation_method), path(genome), path(gff)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_canonical.gff"), emit: gff
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_proteins.fasta"), emit: fasta

    script:
    """
    # Keep the longest isoform per gene
    agat_sp_keep_longest_isoform.pl \
        --gff ${gff} \
        -o ${id}_${annotation_method}_canonical.gff

    # Extract protein sequences
    agat_sp_extract_sequences.pl \
        -g ${id}_${annotation_method}_canonical.gff \
        -f ${genome} \
        -o ${id}_${annotation_method}_proteins.fasta \
        --protein
    """

    stub:
    """
    command -v agat_sp_keep_longest_isoform.pl >/dev/null
    command -v agat_sp_extract_sequences.pl >/dev/null

    touch ${id}_${annotation_method}_canonical.gff
    touch ${id}_${annotation_method}_proteins.fasta
    """
}