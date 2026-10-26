// Portable channel regression: run with -stub-run and test/validation/rnaseq_runs.config.
include { RNASEQ_ASSEMBLY } from './subworkflows/genes_annotation/rnaseq_assembly.nf'

workflow GENES_ANNOTATION {
    take:
    genomes
    main:
    RNASEQ_ASSEMBLY(genomes)
    RNASEQ_ASSEMBLY.out.qc_reports.collect(flat: false).view { rows ->
        def count = params.test_local_libraries ? 8 : 2
        assert rows.size() == count // shared SRA run plus optional local libraries
        assert rows.collect { it[0].id }.toSet().size() == count
        assert rows.count { it[0].run_id == 'SRR000001' } == 2
        "RNASEQ_QC_OK: ${rows.size()} independently cached libraries"
    }
    RNASEQ_ASSEMBLY.out.run_bams.collect(flat: false).ifEmpty { [] }.view { rows ->
        assert rows.size() == (params.test_clean_read_count == 0 ? 0 : params.test_local_libraries ? 8 : 2)
        assert rows.every { it[1].name == "${it[0].id}.sorted.bam" }
        assert rows.groupBy { it[0].genome_id }.values().every { it.size() == (params.test_local_libraries ? 4 : 1) }
        "RNASEQ_SORT_OK: ${rows.size()} independent sorted BAMs"
    }
    RNASEQ_ASSEMBLY.out.bams.collect(flat: false).ifEmpty { [] }.view { rows ->
        assert rows.size() == (params.test_clean_read_count == 0 ? 0 : 2)
        assert rows.collect { it[0].id }.toSet() == (params.test_clean_read_count == 0 ? [] : ['sample_a', 'sample_b']).toSet()
        assert rows.every { !it[0].containsKey('genome_meta') && !it[0].containsKey('run_id') && it[0].rnaseq_aligner == 'star' }
        assert rows.every { params.test_local_libraries ? it[1].name == "${it[0].id}.bam" : it[1].name.endsWith('.sorted.bam') }
        "RNASEQ_MERGE_OK: ${rows.size()} genome BAMs, original metadata restored"
    }
}

workflow AnnoGlyph {
    main:
    def root = "${projectDir}/test/data/rnaseq"
    def genome = file("${root}/genome.fa", checkIfExists: true)
    genomes = Channel.of('sample_a', 'sample_b').map { id ->
        [[id:id, name:'Saccharomyces cerevisiae', assembly_size:30000,
          lineage:[species:[name:'Saccharomyces cerevisiae']],
          'RNASeq-file':params.test_local_libraries ? [single_end:[file("${root}/single.fastq")],
                         paired_end:[[file("${root}/one/R1.fastq"), file("${root}/one/R2.fastq")],
                                     [file("${root}/two/R1.fastq"), file("${root}/two/R2.fastq")]]] : [:]], genome]
    }
    GENES_ANNOTATION(genomes)
}

workflow {
    AnnoGlyph()
}
