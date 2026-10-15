process AGAT_SPKEEPLONGESTISOFORM {
    tag "${meta.annotation_method}/${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.6.1--pl5321hdfd78af_1' :
        'quay.io/biocontainers/agat:1.6.1--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(gxf)
    path config

    output:
    tuple val(meta), path("${output}"), emit: gff
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def config_param = config ? "--config ${config}" : ""
    output           = "${prefix}.gff"
    """
    agat_sp_keep_longest_isoform.pl \\
        --gff ${gxf} \\
        ${config_param} \\
        --out ${output} \\
        ${args}

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        agat: \$(agat --version | sed 's/v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    output     = "${prefix}.gff"
    """
    touch ${output}

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        agat: \$(agat --version | sed 's/v//')
    END_VERSIONS
    """
}
