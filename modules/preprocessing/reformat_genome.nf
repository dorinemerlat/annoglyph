process REFORMAT_GENOME {
    tag "${id}"
    label 'exogap_python'

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path("genome_${id}.fa"), emit: fasta
    tuple val(id), val(meta), path("genome_${id}.fa.tsv"), emit: mapping

    script:
    def unmodified_genome = "${id}_unreformated.fa"
    def genome_code = meta.mnemonic.exogap

    """
    mv ${genome} ${unmodified_genome}

    reformat_genome.py \
        --fasta ${unmodified_genome} \
        --prefix ${genome_code} \
        --output genome_${id}.fa \
        --mapping genome_${id}.fa.tsv
    """

    stub:
    """
    # Check that the script is available
    command -v reformat_genome.py >/dev/null
    
    cat > genome_${id}.fa << EOF
>${genome_code}_1
ATGC
EOF

    cat > genome_${id}.fa.tsv << EOF
previous_IDs	EXOGAP_IDs
stub_sequence	${genome_code}_1
EOF
    """
}