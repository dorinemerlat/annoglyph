process REPEATMASKER {
    tag "$id"
    cpus 30
    time '20d'
    label 'repeatmasker'
    errorStrategy 'ignore'
    scratch false

    input:
    tuple val(id), val(meta), path(genome), path(library)

    output:
    tuple val(id), val(meta), path("genome_${id}.fa.masked"),  emit: masked
    tuple val(id), val(meta), path("genome_${id}.fa.cat.gz"),  emit: cat
    tuple val(id), val(meta), path("genome_${id}.fa.out"),     emit: out
    tuple val(id), val(meta), path("genome_${id}.fa.out.gff"), emit: gff
    tuple val(id), val(meta), path("genome_${id}.fa.tbl"),     emit: tbl
    tuple val(id), val(meta), path("genome_${id}.fa.align"),   emit: align

    script:
    def staged_genome = "genome_${id}.fa"
    """
    ln -s ${genome} ${staged_genome}
    RepeatMasker -pa ${task.cpus} -e ncbi -lib $library -a -gccalc -norna -excln  -s -xsmall -gff ${staged_genome}
    """

    stub:
    """
    touch ${staged_genome}.masked ${staged_genome}.cat ${staged_genome}.out ${staged_genome}.out.gff ${staged_genome}.tbl ${staged_genome}.align
    """
}
