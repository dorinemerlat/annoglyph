process INTERPROSCAN {
    tag "${meta.annotation_method}/${meta.id}"
    // will throw NullPointer exceptions and crush with more than 1 cpu
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/interproscan:5.59_91.0--hec16e2b_1' :
        'quay.io/biocontainers/interproscan:5.59_91.0--hec16e2b_1' }"

    input:
    tuple val(meta), path(fasta)
    path(interproscan_database, stageAs: 'data')

    output:
    tuple val(meta), path('*.tsv') , optional: true, emit: tsv
    tuple val(meta), path('*.xml') , optional: true, emit: xml
    tuple val(meta), path('*.gff3'), optional: true, emit: gff3
    tuple val(meta), path('*.json'), optional: true, emit: json
    path "versions.yml", emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args             = task.ext.args ?: ''
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def is_compressed    = fasta.getExtension() == "gz"
    def fasta_name       = is_compressed ? fasta.getBaseName() : fasta
    def uncompress_input = is_compressed ? "gzip -c -d ${fasta} > ${fasta_name}" : ''
    """
    $uncompress_input

    if [ -d 'data' ]; then
        # Find interproscan.properties to link data/ from work directory
        INTERPROSCAN_DIR="\$( dirname "\$( dirname "\$( which interproscan.sh )" )" )"
        INTERPROSCAN_PROPERTIES="\$( find "\$INTERPROSCAN_DIR/share" -name "interproscan.properties" )"
        cp "\$INTERPROSCAN_PROPERTIES" .
        sed -i \\
            -e "s|^bin\\.directory=.*|bin.directory=\$( dirname "\$INTERPROSCAN_PROPERTIES" )/bin|" \\
            -e "s|^data\\.directory=.*|data.directory=\$PWD/data|" \\
            interproscan.properties
        # must be absolute: interproscan.sh cd's into its install dir before reading it
        export INTERPROSCAN_CONF=\$PWD/interproscan.properties
    elif [ '${task.ext.use_container_data ?: false}' != 'true' ]; then
        echo 'Provide --interproscan_db with the complete InterProScan 5.59-91.0 data directory; the image only contains test data.' >&2
        exit 1
    fi

    # The retained image uses BusyBox ps, lacking the procps syntax used by IPS.
    if ! ps -o pid,vsize,rss,command -p \$\$ >/dev/null 2>&1; then
        mkdir -p interproscan_compat
        cp ${moduleDir}/resources/ps interproscan_compat/ps
        chmod +x interproscan_compat/ps
        export PATH="\$PWD/interproscan_compat:\$PATH"
    fi

    # AGAT protein FASTAs can contain stop symbols, which InterProScan rejects.
    sed '/^>/!s/\\*//g' ${fasta_name} > ${prefix}.input.fasta

    interproscan.sh \\
        --cpu ${task.cpus} \\
        --input ${prefix}.input.fasta \\
        ${args} \\
        --output-file-base ${prefix}

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        interproscan: \$(interproscan.sh --version | sed '1!d; s/.*version //')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.{tsv,xml,json,gff3}

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        interproscan: \$(interproscan.sh --version | sed '1!d; s/.*version //')
    END_VERSIONS
    """
}
