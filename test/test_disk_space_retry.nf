// Run from a separate launch directory:
// nextflow -c /path/to/annoglyph/conf/modules.config run /path/to/annoglyph/test/test_disk_space_retry.nf
nextflow.enable.dsl = 2

params.program_args = [:]
// Other validation modes: ordinary (no retry), persistent (exhaust maxRetries).
params.failure_mode = 'recover'

process WRAPPER_DISK_FAILURE {
    // Fail before the tool starts: the diagnostic is only in .command.log.
    beforeScript { task.attempt == 1 ? 'echo "mktemp: No space left on device" >&2; exit 1' : ':' }
    output:
    path 'wrapper.txt'
    script:
    """
    echo '${task.attempt}' > wrapper.txt
    """
}

process TOOL_DISK_FAILURE {
    input:
    val stream
    output:
    path "${stream}.txt"
    script:
    def redirect = stream == 'stderr' ? '>&2' : ''
    def message = params.failure_mode == 'ordinary' ? 'invalid input data' : 'write failed: No space left on device'
    """
    if [ '${task.attempt}' = 1 ] || [ '${params.failure_mode}' != recover ]; then
        echo '${message}' ${redirect}
        exit 1
    fi
    echo '${task.attempt}' > ${stream}.txt
    """
}

workflow {
    WRAPPER_DISK_FAILURE()
    TOOL_DISK_FAILURE(Channel.of('stdout', 'stderr'))
    WRAPPER_DISK_FAILURE.out.view { assert it.text.trim() == '2'; 'Wrapper failure recovered' }
    TOOL_DISK_FAILURE.out.view { assert it.text.trim() == '2'; "Tool failure recovered: ${it.name}" }
}
