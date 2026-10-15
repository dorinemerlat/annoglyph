class ProvenanceUtils {
    static def plain(value) {
        if (value instanceof Map) return value.collectEntries { key, item ->
            [(key.toString()): key.toString().toLowerCase() ==~ /.*(password|secret|token|credential).*/ ? '[REDACTED]' : plain(item)]
        }
        if (value instanceof List) return value.collect { plain(it) }
        if (value == null || value instanceof Number || value instanceof Boolean || value instanceof String) return value
        return value.toString()
    }

    static Map snapshot(params, workflow, nextflowVersion) {
        def blast = params.blast_db ?: ''
        def blastPath = new File(blast.toString())
        return [params: plain(params),
            workflow: [run_name: workflow.runName, session_id: workflow.sessionId.toString(),
                nextflow_version: nextflowVersion.toString(), revision: workflow.revision,
                profile: workflow.profile, project_dir: workflow.projectDir.toString()],
            singularity_cache: System.getenv('NXF_SINGULARITY_CACHEDIR') ?: "${workflow.projectDir}/.singularity",
            databases: [swissprot_blast: blastPath.isAbsolute() ? blastPath.parent : null,
                eggnog: params.eggnog_db, interproscan_data: params.interproscan_db,
                interproscan_container: params.interproscan_container]]
    }
}
