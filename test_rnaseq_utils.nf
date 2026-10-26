workflow {
    def root = "${projectDir}/test/data/rnaseq"
    def genome = [id:'genome_a', name:'Original genome']
    def pair = [file("${root}/one/R1.fastq"), file("${root}/one/R2.fastq")]
    def first = RnaseqUtils.readset(genome, pair, 'sra', 'SRR123')
    assert first.run_id == 'SRR123' && first.genome_id == 'genome_a' && !first.single_end
    assert first == RnaseqUtils.readset(genome, pair, 'sra', 'SRR123')
    assert !genome.containsKey('run_id') && first.genome_meta == genome
    def other = RnaseqUtils.readset(genome, [file("${root}/two/R1.fastq"), file("${root}/two/R2.fastq")], 'local_pe')
    assert first.id != other.id
    assert RnaseqUtils.mergeMeta([first, other]) == RnaseqUtils.mergeMeta([other, first])
    assert RnaseqUtils.mergeMeta([first]).rnaseq_run_signature != RnaseqUtils.mergeMeta([first, other]).rnaseq_run_signature
    def foreign = RnaseqUtils.readset([id:'genome_b'], pair, 'sra', 'SRR123')
    assert foreign.id != first.id
    try {
        RnaseqUtils.mergeMeta([first, foreign])
        assert false: 'Mixed genomes must not merge'
    } catch (IllegalArgumentException expected) {}
    try {
        RnaseqUtils.readset(genome, pair + [pair[0]], 'local')
        assert false: 'Invalid library arity must be rejected'
    } catch (IllegalArgumentException expected) {}
    assert RnaseqUtils.starMemoryGb(0L, 0L, params) == params.rnaseq_star_min_memory_gb
    assert RnaseqUtils.starMemoryGb(10L * 1024 * 1024 * 1024, 100_000_000L, params) > RnaseqUtils.starMemoryGb(0L, 100L, params)
    assert RnaseqUtils.indexBytes(file(root)) > 30000
    def json = File.createTempFile('rnaseq-count-', '.json')
    try {
        json.text = '{"summary":{"after_filtering":{"total_reads":5000000000,"total_bases":750000000000}}}'
        assert RnaseqUtils.readStats(json).rnaseq_read_count == 5_000_000_000L
        json.text = '{"summary":{"after_filtering":{"total_reads":0,"total_bases":0}}}'
        assert RnaseqUtils.readStats(json).rnaseq_read_count == 0
    } finally { json.delete() }
    println 'RNASEQ_UTILS_OK: stable identities, genome isolation, cache keys and resource estimates'
}
