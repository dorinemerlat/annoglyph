process REPEATMASKER {
    tag "${meta.id}"
    label 'repeatmasker'

    input:
    tuple val(meta), path(genome), path(library)

    output:
    tuple val(meta), path("genome_${meta.id}.fa.masked"), emit: masked
    tuple val(meta), path("genome_${meta.id}.fa.cat*"), emit: cat
    tuple val(meta), path("genome_${meta.id}.fa.out"), emit: out
    tuple val(meta), path("genome_${meta.id}.fa.out.gff"), emit: gff
    tuple val(meta), path("genome_${meta.id}.fa.tbl"), emit: tbl
    tuple val(meta), path("genome_${meta.id}.fa.align"), emit: align
    tuple val(meta), path("${meta.id}_repeatmasker.log"), emit: log

    script:
    def staged_genome = "genome_${meta.id}.fa"

    """
    set -o pipefail

    # Stage the genome file to a new name to avoid issues with RepeatMasker output files
    ln -s ${genome} ${staged_genome}

    RepeatMasker \
        -pa ${task.cpus} \
        -e ncbi \
        -lib ${library} \
        -a \
        -gccalc \
        -norna \
        -excln \
        -s \
        -xsmall \
        -gff \
        ${staged_genome} \
        2>&1 | tee -a ${meta.id}_repeatmasker.log
    """

    stub:
    def staged_genome = "genome_${meta.id}.fa"
    """
    touch ${staged_genome}.masked
    touch ${staged_genome}.cat.gz
    touch ${staged_genome}.out
    touch ${staged_genome}.out.gff
    touch ${staged_genome}.tbl
    touch ${staged_genome}.align
    """
}
