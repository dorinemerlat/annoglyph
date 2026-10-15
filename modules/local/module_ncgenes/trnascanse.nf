process TRNASCANSE {
    tag "${meta.id}"
    label 'trnascan'
    cpus 8
    scratch false
    label 'program_trnascan'
    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${meta.id}_trna.gff"), emit: gff
    tuple val(meta), path("${meta.id}_trna.fa"), emit: fasta
    tuple val(meta), path("${meta.id}_trna.out"), emit: out
    tuple val(meta), path("${meta.id}_trna.struct"), emit: struct
    tuple val(meta), path("${meta.id}_trna.stats"), emit: stats

    script:
    """
    export TMPDIR=\$PWD/tmp
    mkdir -p \$TMPDIR

    tRNAscan-SE -E \
        --thread ${task.cpus} \
        --output ${meta.id}_trna.out \
        --struct ${meta.id}_trna.struct \
        --gff ${meta.id}_trna.gff \
        --fasta ${meta.id}_trna.fa \
        --stats ${meta.id}_trna.stats \
        --progress \
        ${params.program_args.trnascan ?: ''} ${genome}
    # --forceow
    """

    stub:
    """
    touch ${meta.id}_trna.gff ${meta.id}_trna.fa ${meta.id}_trna.out ${meta.id}_trna.struct ${meta.id}_trna.stats
    """
}
