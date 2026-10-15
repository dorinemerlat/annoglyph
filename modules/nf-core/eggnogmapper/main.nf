process EGGNOGMAPPER {
    tag "${meta.annotation_method}/${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.13--pyhdfd78af_2':
        'quay.io/biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(search_mode), path(db)
    path(eggnog_data_dir)

    output:
    tuple val(meta), path("${prefix}/*.emapper.annotations")   , emit: annotations
    tuple val(meta), path("${prefix}/*.emapper.seed_orthologs"), emit: seed_orthologs, optional: true
    tuple val(meta), path("${prefix}/*.emapper.hits")          , emit: hits     , optional: true
    tuple val(meta), path("${prefix}/*.emapper.orthologs"), emit: orthologs
    // emapper.py --version reports the enclosing repository's git tag when the package sits inside one, as it does under conda.
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ''
    prefix            = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.extension == '.gz'
    def fasta_name    = is_compressed ? fasta.baseName : "$fasta"
    def db_flags = ['diamond': '--dmnd_db', 'novel_fams': '--dmnd_db', 'mmseqs': '--mmseqs_db', 'hmmer': '--database', 'no_search': '--annotate_hits_table', 'cache': '--cache']
    def db_path  = (db instanceof Path && db.isDirectory()) ? "${db}/${db.name}"                    : "$db"
    def db_arg   = db && db_flags[search_mode]              ? "${db_flags[search_mode]} ${db_path}" : ''
    def dbmem    = task.ext.dbmem ? '--dbmem' : ''
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    mkdir -p ${prefix}
    emapper.py \\
        $args \\
        --cpu ${task.cpus} \\
        -i ${fasta_name} \\
        --data_dir ${eggnog_data_dir} \\
        -m ${search_mode} \\
        $db_arg \\
        ${dbmem} \\
        --output_dir ${prefix} \\
        --output ${prefix}

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        eggnog-mapper: \$(python -c "from importlib.metadata import version; print(version('eggnog-mapper'))")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    touch ${prefix}/${prefix}.emapper.annotations
    touch ${prefix}/${prefix}.emapper.seed_orthologs
    touch ${prefix}/${prefix}.emapper.hits
    touch ${prefix}/${prefix}.emapper.orthologs

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        eggnog-mapper: \$(python -c "from importlib.metadata import version; print(version('eggnog-mapper'))")
    END_VERSIONS
    """
}
