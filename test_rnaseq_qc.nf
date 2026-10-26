// Real fastp checks. See test/README.md for default, disabled and stricter runs.
include { CLEAN_RNASEQ_READS } from './modules/local/genes_annotation/rnaseq_assembly/clean_reads.nf'
params.test_qc_expected_good = 30
params.test_qc_expected_mixed = 1

workflow GENES_ANNOTATION {
    main:
    def root = "${projectDir}/test/data/rnaseq"
    def genome = [id:'qc_genome']
    libraries = Channel.of(
        ['good_se', [file("${root}/single.fastq")]],
        ['good_pe', [file("${root}/one/R1.fastq"), file("${root}/one/R2.fastq")]],
        ['mixed_se', [file("${root}/mixed/single.fastq")]],
        ['mixed_pe', [file("${root}/mixed/R1.fastq"), file("${root}/mixed/R2.fastq")]]
    ).map { name, reads -> [RnaseqUtils.readset(genome, reads, 'fixture', name), reads] }
    CLEAN_RNASEQ_READS(libraries)
    CLEAN_RNASEQ_READS.out.reads.collect(flat:false).view { rows ->
        assert rows.size() == 4
        rows.each { meta, reads, report ->
            def files = (reads instanceof List ? reads : [reads]).sort { it.name }
            def expected = (meta.run_id.startsWith('mixed') ? params.test_qc_expected_mixed : params.test_qc_expected_good) as int
            assert RnaseqUtils.readStats(report).rnaseq_read_count == expected * files.size()
            def ids = files.collect { path ->
                def stream = new java.util.zip.GZIPInputStream(new FileInputStream(path.toString()))
                def lines = stream.withReader { it.readLines() }
                assert lines.size() == expected * 4
                lines.collate(4).collect { record ->
                    assert record.size() == 4 && record[0].startsWith('@') && record[2].startsWith('+')
                    assert record[1].size() == record[3].size()
                    record[0].split(/\s+/)[0].replaceAll(/\/[12]$/, '')
                }
            }
            if (!meta.single_end) assert ids[0] == ids[1]
        }
        'RNASEQ_REAL_QC_OK: filtering, empty libraries and paired synchronization verified'
    }
}

workflow AnnoGlyph {
    main:
    GENES_ANNOTATION()
}

workflow {
    AnnoGlyph()
}
