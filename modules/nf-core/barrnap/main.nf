process BARRNAP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/46925ba4c44fde31056890ad1af51140c91616055f8b7ac073f51db21d321e32/data'
        : 'community.wave.seqera.io/library/barrnap:0.9--847dcfa21d3169b7' }"

    input:
    tuple val(meta), path(fasta), val(dbname)

    output:
    path "versions.yml", emit: versions
    tuple val(meta), path("*.gff"), emit: gff

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    db         = dbname ? "${dbname}" : 'bac'
    input    = fasta =~ /\.gz$/ ? fasta.name.take(fasta.name.lastIndexOf('.')) : fasta
    gunzip   = fasta =~ /\.gz$/ ? "gunzip -c ${fasta} > ${input}" : ""

    """
    $gunzip

    barrnap \\
        $args \\
        --threads $task.cpus \\
        --kingdom $db \\
        $input \\
        > ${prefix}_${db}.gff


    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        barrnap: \$( barrnap --version 2>&1 | sed "s/barrnap //g" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    db = dbname ? "${dbname}" : 'bac'
    """
    touch ${prefix}_${db}.gff

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        barrnap: \$( barrnap --version 2>&1 | sed "s/barrnap //g" )
    END_VERSIONS
    """
}
