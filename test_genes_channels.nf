include { RNASEQ_ASSEMBLY } from './subworkflows/genes_annotation/rnaseq_assembly.nf'
include { STRUCTURAL_ANNOTATION } from './subworkflows/genes_annotation/structural_annotation.nf'
include { FUNCTIONAL_ANNOTATION } from './subworkflows/genes_annotation/functional_annotation.nf'

workflow {
    def genome = file("${projectDir}/test/test_full2/cache/preprocessing/reformat_genome/saccharomyces-cerevisiae-s288c/genome_saccharomyces-cerevisiae-s288c.fa", checkIfExists: true)
    def reads = file("${projectDir}/test/validation/genes_reads.fastq", checkIfExists: true)
    genomes = Channel.of(
        [[id:'sample_a', name:'Saccharomyces cerevisiae', mnemonic:[annoglyph:'AAAA'], lineage:[species:[name:'Saccharomyces cerevisiae']], 'RNASeq-file':[single_end:[reads], paired_end:[[reads, reads]]]], genome],
        [[id:'sample_b', name:'Saccharomyces cerevisiae', mnemonic:[annoglyph:'BBBB'], lineage:[species:[name:'Saccharomyces cerevisiae']], 'RNASeq-file':[single_end:[reads], paired_end:[[reads, reads]]]], genome]
    )
    RNASEQ_ASSEMBLY(genomes)
    STRUCTURAL_ANNOTATION(genomes, RNASEQ_ASSEMBLY.out.bams, Channel.value(reads), Channel.value('fungi_odb12'))
    FUNCTIONAL_ANNOTATION(STRUCTURAL_ANNOTATION.out.gff, STRUCTURAL_ANNOTATION.out.proteins)
    FUNCTIONAL_ANNOTATION.out.gff.collect(flat: false).view { rows ->
        def expectedMethods = params.test_read_count <= 1000000 ? ['braker2','galba'] : params.rna_seq_aligner == 'both' ? ['braker2','galba','braker-with-star','braker-with-varus'] : ['braker2','galba',"braker-with-${params.rna_seq_aligner}"]
        assert rows.size() == 2 * expectedMethods.size()
        assert rows.collect { it[0].id }.toSet() == ['sample_a','sample_b'].toSet()
        assert rows.collect { it[0].annotation_method }.toSet() == expectedMethods.toSet()
        rows.each { meta, gff -> assert gff.name == "${meta.id}_${meta.annotation_method}_with_functional_annotations.gff" }
        "GENE_CHANNELS_OK: ${rows.size()} annotations with correct sample/method associations"
    }
}
