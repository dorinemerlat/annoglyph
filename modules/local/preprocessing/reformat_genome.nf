process REFORMAT_GENOME {
    tag "${meta.id}"
    label 'annoglyph_python'

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("genome_${meta.id}.fa"), emit: fasta
    tuple val(meta), path("genome_${meta.id}.fa.tsv"), emit: mapping

    script:
    def unmodified_genome = "${meta.id}_unreformated.fa"
    def genome_code = meta.mnemonic.annoglyph

    """
    mv ${genome} ${unmodified_genome}

    reformat_genome.py \
        --fasta ${unmodified_genome} \
        --prefix ${genome_code} \
        --output genome_${meta.id}.fa \
        --mapping genome_${meta.id}.fa.tsv
    """

    stub:
    def genome_code = meta.mnemonic.annoglyph
    """
    # Check that the script is available
    command -v reformat_genome.py >/dev/null
    
    cat > genome_${meta.id}.fa << EOF
>${genome_code}_1
ATGC
EOF

    cat > genome_${meta.id}.fa.tsv << EOF
previous_IDs	AnnoGlyph_IDs
stub_sequence	${genome_code}_1
EOF
    """
}