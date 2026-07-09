process TRNASCANSE {
    tag "${id}"
    label 'trnascan'
    cpus 8
    scratch false
    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path("${id}_trna.gff"), emit: gff
    tuple val(id), val(meta), path("${id}_trna.fa"), emit: fasta
    tuple val(id), val(meta), path("${id}_trna.out"), emit: out
    tuple val(id), val(meta), path("${id}_trna.struct"), emit: struct
    tuple val(id), val(meta), path("${id}_trna.stats"), emit: stats

    script:
    """
    export TMPDIR=\$PWD/tmp
    mkdir -p \$TMPDIR

    tRNAscan-SE -E \
        --thread ${task.cpus} \
        --output ${id}_trna.out \
        --struct ${id}_trna.struct \
        --gff ${id}_trna.gff \
        --fasta ${id}_trna.fa \
        --stats ${id}_trna.stats \
        --progress \
        ${genome}
    # --forceow
    """

    stub:
    """
    touch ${id}_trna.gff ${id}_trna.fa ${id}_trna.out ${id}_trna.struct ${id}_trna.stats
    """
}
