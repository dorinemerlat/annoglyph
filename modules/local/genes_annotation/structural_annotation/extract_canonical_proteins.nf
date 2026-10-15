process EXTRACT_CANONICAL_PROTEINS {
    tag "${meta.annotation_method}/${meta.id}"
    label 'agat'

    input:
    tuple val(meta), path(genome), path(gff)

    output:
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_canonical.gff"), emit: gff
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_proteins.fasta"), emit: fasta

    script:
    """
    # Keep the longest isoform per gene
    agat_sp_keep_longest_isoform.pl \
        --gff ${gff} \
        -o ${meta.id}_${meta.annotation_method}_canonical.gff

    # Extract protein sequences
    agat_sp_extract_sequences.pl \
        -g ${meta.id}_${meta.annotation_method}_canonical.gff \
        -f ${genome} \
        -o ${meta.id}_${meta.annotation_method}_proteins.fasta \
        --protein
    """

    stub:
    """
    command -v agat_sp_keep_longest_isoform.pl >/dev/null
    command -v agat_sp_extract_sequences.pl >/dev/null

    touch ${meta.id}_${meta.annotation_method}_canonical.gff
    touch ${meta.id}_${meta.annotation_method}_proteins.fasta
    """
}