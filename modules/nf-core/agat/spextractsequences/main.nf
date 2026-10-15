process AGAT_SPEXTRACTSEQUENCES {
    tag "${meta.annotation_method}/${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.6.1--pl5321hdfd78af_1' :
        'quay.io/biocontainers/agat:1.6.1--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(gxf), path(fasta)
    path config

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def config_arg  = config ? "-c ${config}" : ''
    if( "$fasta" == "${prefix}.fasta" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    agat_sp_extract_sequences.pl \\
        ${args} \\
        -g ${gxf} \\
        -f ${fasta} \\
        ${config_arg} \\
        -o ${prefix}.fasta

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        agat: \$(agat --version | sed 's/v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "${fasta}" == "${prefix}.fasta" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.fasta

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        agat: \$(agat --version | sed 's/v//')
    END_VERSIONS
    """
}
