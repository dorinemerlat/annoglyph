process REPEATMODELER {
    tag "${id}"
    cpus 50
    time '20d'
    label 'repeatmodeler'

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path("${id}-families.fa")

    script:
    """
    BuildDatabase -name $id -engine ncbi $genome

    RepeatModeler -threads ${task.cpus} -engine ncbi -database $id -LTRStruct
    """

    stub:
    """
    touch ${id}_families.fa
    """
}
