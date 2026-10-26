// Real tiny-data validation: fastp -> STAR -> samtools sort -> merge.
include { CLEAN_RNASEQ_READS } from './modules/local/genes_annotation/rnaseq_assembly/clean_reads.nf'
include { STAR_GENOMEGENERATE } from './modules/nf-core/star/genomegenerate/main.nf'
include { STAR_ALIGN as STAR_SINGLE; STAR_ALIGN as STAR_PAIRED } from './modules/nf-core/star/align/main.nf'
include { SORT_RNASEQ_BAM } from './modules/local/genes_annotation/rnaseq_assembly/sort_bam.nf'
include { SAMTOOLS_MERGE } from './modules/nf-core/samtools/merge/main.nf'

process CHECK_RNASEQ_BAM {
    container 'docker://dorinemerlat/varus:v1.1'
    input:
    path bam
    output:
    path 'passed.txt'
    script:
    """
    samtools quickcheck '${bam}'
    test "\$(samtools view -c '${bam}')" -eq 90
    samtools view -H '${bam}' | grep 'SO:coordinate'
    samtools view '${bam}' | awk 'BEGIN { prev=0 } { if (\$4 < prev) exit 1; prev=\$4 }'
    echo 'RNASEQ_REAL_TOOLS_OK: 90 mapped reads, coordinate-sorted merged BAM' > passed.txt
    """
}

workflow GENES_ANNOTATION {
    main:
    def root = "${projectDir}/test/data/rnaseq"
    def genomeMeta = [id:'tiny_genome', assembly_size:30000]
    def se = [file("${root}/single.fastq")]
    def pe = [file("${root}/one/R1.fastq"), file("${root}/one/R2.fastq")]
    CLEAN_RNASEQ_READS(Channel.of([RnaseqUtils.readset(genomeMeta, se, 'local_se'), se],
                                  [RnaseqUtils.readset(genomeMeta, pe, 'local_pe'), pe]))
    STAR_GENOMEGENERATE(Channel.of([genomeMeta, file("${root}/genome.fa")]), Channel.value([[:], []]))
    indexed = CLEAN_RNASEQ_READS.out.reads.map { meta, reads, report ->
        def stats = RnaseqUtils.readStats(report)
        assert stats.rnaseq_read_count == (meta.single_end ? 30 : 60)
        [meta + stats, (reads instanceof List ? reads : [reads]).sort { it.name }]
    }.combine(STAR_GENOMEGENERATE.out.index).map { meta, reads, indexMeta, index ->
        [meta + [star_memory_gb:RnaseqUtils.starMemoryGb(RnaseqUtils.indexBytes(index), meta.rnaseq_read_count as long, params)], reads, index]
    }.branch {
        single: it[0].single_end
        paired: !it[0].single_end
    }
    singles = indexed.single.multiMap { meta, reads, index ->
        reads: [meta, reads]
        index: [[id:meta.genome_id], index]
    }
    pairs = indexed.paired.multiMap { meta, reads, index ->
        reads: [meta, reads]
        index: [[id:meta.genome_id], index]
    }
    STAR_SINGLE(singles.reads, singles.index, Channel.value([[:], []]), true)
    STAR_PAIRED(pairs.reads, pairs.index, Channel.value([[:], []]), true)
    SORT_RNASEQ_BAM(STAR_SINGLE.out.bam.concat(STAR_PAIRED.out.bam))
    merge = SORT_RNASEQ_BAM.out.bam.collect(flat:false).map { rows ->
        [RnaseqUtils.mergeMeta(rows.collect { it[0] }), rows.collect { it[1] }.sort { it.name }, []]
    }
    SAMTOOLS_MERGE(merge, Channel.value([[:], [], [], []]), '')
    CHECK_RNASEQ_BAM(SAMTOOLS_MERGE.out.bam.map { meta, bam -> bam })
    CHECK_RNASEQ_BAM.out.view { it.text }
}

workflow AnnoGlyph {
    main:
    GENES_ANNOTATION()
}

workflow {
    AnnoGlyph()
}
