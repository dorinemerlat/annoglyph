process COLLECT_PROVENANCE {
    tag 'run_provenance'
    container null
    storeDir null
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode, overwrite: true
    input:
    val snapshot
    tuple val(input_names), val(original_paths), path(input_files, stageAs: 'inputs??/*')
    path version_files, stageAs: 'versions??/*'
    tuple val(database_names), path(database_artifacts, stageAs: 'databases??/*')
    path lock
    path provenance_script
    val completion
    output:
    path 'provenance', emit: report
    script:
    def normalize = { value -> value instanceof List ? value : (value ? [value] : []) }
    def inputs = normalize(input_files).withIndex().collect { path, index -> [name: input_names[index], original: original_paths[index], staged: path.toString()] }
    def config = snapshot + [used_database_artifacts: normalize(database_artifacts).withIndex().collect { path, index -> [name: database_names[index], staged: path.toString()] }]
    def snapshotJson = groovy.json.JsonOutput.toJson(config)
    def inputsJson = groovy.json.JsonOutput.toJson(inputs)
    def versionsJson = groovy.json.JsonOutput.toJson(normalize(version_files).collect { it.toString() })
    def lockArg = lock ? "--database-lock ${lock}" : ''
    """
    cat > snapshot.json <<'SNAPSHOT'
    ${snapshotJson}
    SNAPSHOT
    cat > inputs.json <<'INPUTS'
    ${inputsJson}
    INPUTS
    cat > versions.json <<'VERSIONS'
    ${versionsJson}
    VERSIONS
    python3 ${provenance_script} --snapshot snapshot.json --inputs inputs.json --versions versions.json ${lockArg} --output provenance
    """
    stub:
    """
    mkdir provenance
    echo '{}' > provenance/parameters.json
    """
}
