include { SAMTOOLS_MERGE } from './modules/nf-core/samtools/merge/main.nf'
workflow {
    inputs = Channel.of([[id:'merge_smoke'], [file("${projectDir}/test/validation/genes_merge_a.bam"), file("${projectDir}/test/validation/genes_merge_b.bam")], []])
    SAMTOOLS_MERGE(inputs, Channel.value([[:], [], [], []]), '')
    SAMTOOLS_MERGE.out.bam.view { meta, bam ->
        assert bam.size() > 0
        "MERGE_REAL_OK: ${bam.name}"
    }
}
