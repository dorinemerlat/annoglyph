process BARRNAP {

    tag "${kingdom}/${meta.id}"
    label 'exogap_tools'
    cpus 30

    label 'program_barrnap'
    input:
    tuple val(meta), path(genome), val(kingdom)

    output:
    tuple val(meta), path("${meta.id}_barrnap_${kingdom}.gff"), emit: gff

    script:
    """
    barrnap --quiet --kingdom ${kingdom} --threads ${task.cpus} ${params.program_args.barrnap ?: ''} "${genome}"  > barrnap.out

    cp barrnap.out "${meta.id}_barrnap_${kingdom}.gff"
    """

    stub:
    """
    touch ${meta.id}_barrnap_${kingdom}.gff
    """
}