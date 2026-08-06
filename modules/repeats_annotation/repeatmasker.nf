process REPEATMASKER {
    tag "${id}"
    label 'repeatmasker'
    scratch false
    input:
    tuple val(id), val(meta), path(genome), path(library)

    output:
    tuple val(id), val(meta), path("genome_${id}.fa.masked"), emit: masked
    tuple val(id), val(meta), path("genome_${id}.fa.cat*"), emit: cat
    tuple val(id), val(meta), path("genome_${id}.fa.out"), emit: out
    tuple val(id), val(meta), path("genome_${id}.fa.out.gff"), emit: gff
    tuple val(id), val(meta), path("genome_${id}.fa.tbl"), emit: tbl
    tuple val(id), val(meta), path("genome_${id}.fa.align"), emit: align
    tuple val(id), val(meta), path("${id}_repeatmasker.log"), emit: log

    script:
    def staged_genome = "genome_${id}.fa"

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
        2>&1 | tee -a ${id}_repeatmasker.log
    """

    stub:
    """
    touch ${staged_genome}.masked
    touch ${staged_genome}.cat.gz
    touch ${staged_genome}.out
    touch ${staged_genome}.out.gff
    touch ${staged_genome}.tbl
    touch ${staged_genome}.align
    """
}