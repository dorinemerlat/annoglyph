// Simulate STAR's exact error without allocating a large BAM or consuming RAM.
nextflow.enable.dsl = 2
params.program_args = [:]
params.star_failure = 'reported'

process STAR_MEMORY_FAILURE {
    output:
    path 'passed.txt'
    script:
    long memory = task.memory.toBytes()
    if (task.attempt > 1) {
        long expected = params.star_failure in ['reported', 'mixed'] ? 7L * 1024 * 1024 * 1024 : 4L * 1024 * 1024 * 1024
        assert memory == expected : "Expected ${expected} bytes on retry, got ${memory}"
    }
    def diagnostic = params.star_failure == 'ordinary' ? 'invalid STAR argument' :
        'EXITING because of fatal ERROR: not enough memory for BAM sorting:'
    def solution = params.star_failure in ['reported', 'mixed'] ? 'SOLUTION: re-run STAR with at least --limitBAMsortRAM 5368709120' : ''
    """
    if [ '${task.attempt}' = 1 ]; then
        echo '${diagnostic}' >&2
        echo '${solution}' >&2
        exit 102
    fi
    if [ '${params.star_failure}' = mixed ] && [ '${task.attempt}' = 2 ]; then
        echo 'Disk quota exceeded' >&2
        exit 1
    fi
    echo 'STAR_MEMORY_RETRY_OK: attempt ${task.attempt}, ${memory} bytes' > passed.txt
    """
}

workflow {
    STAR_MEMORY_FAILURE()
    STAR_MEMORY_FAILURE.out.view { it.text }
}
