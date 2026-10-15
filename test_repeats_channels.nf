nextflow.enable.dsl = 2
include { REPEATS_ANNOTATION } from './subworkflows/repeats_annotation'
workflow {
    genomes = Channel.of(
        [[id: 'test-cerevisiae', mnemonic: [exogap: 'SACCE'], assembly_size: '12157105', busco_dataset: 'saccharomycetaceae_odb12'], file("${projectDir}/test/test_full2/genomes/GCF_000146045.2/genome.fa")],
        [[id: 'test-paradoxus', mnemonic: [exogap: 'SACPA'], assembly_size: '11500000', busco_dataset: 'saccharomycetaceae_odb12'], file("${projectDir}/test/test_full2/genomes/GCF_002079055.1/genome.fa")]
    )
    datasets = Channel.of([[id: 'saccharomycetaceae_odb12'], file("${projectDir}/test/test_full2/cache/preprocessing/busco_download/saccharomycetaceae_odb12/busco_downloads/lineages/saccharomycetaceae_odb12")])
    REPEATS_ANNOTATION(genomes, Channel.empty(), Channel.empty(), datasets)
    REPEATS_ANNOTATION.out.masked_genomes.toList().map { results ->
        assert results.size() == 2 : 'Both genomes must produce a masked assembly'
        assert results.collect { it[0].id }.toSet() == ['test-cerevisiae', 'test-paradoxus'].toSet()
        results.each { meta, fasta ->
            assert meta instanceof Map
            assert fasta.name == "genome_${meta.id}.fa.masked"
        }
        "Validated two masked assemblies (group_consensus_sequences=${params.group_consensus_sequences})"
    }.view()
    REPEATS_ANNOTATION.out.repeat_gff.toList().map { results ->
        assert results.size() == 2 : 'Both genomes must produce repeat annotations'
        results.each { meta, gff ->
            assert meta instanceof Map
            assert gff.name == "${meta.id}_repeats.gff"
        }
        'Validated two repeat annotation files'
    }.view()
}
