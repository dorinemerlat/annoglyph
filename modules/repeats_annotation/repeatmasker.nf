process REPEATMASKER {
    tag "${id}"
    label 'repeatmasker'

    input:
    tuple val(id), val(meta), path(genome), path(library)

    output:
    tuple val(id), val(meta), path("genome_${id}.fa.masked"), emit: masked
    tuple val(id), val(meta), path("genome_${id}.fa.cat.gz"), emit: cat
    tuple val(id), val(meta), path("genome_${id}.fa.out"), emit: out
    tuple val(id), val(meta), path("genome_${id}.fa.out.gff"), emit: gff
    tuple val(id), val(meta), path("genome_${id}.fa.tbl"), emit: tbl
    tuple val(id), val(meta), path("genome_${id}.fa.align"), emit: align

    script:
    def staged_genome = "genome_${id}.fa"

    """
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
        ${staged_genome}
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