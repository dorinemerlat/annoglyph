import groovy.io.FileType
import groovy.json.JsonOutput
import groovy.json.JsonSlurperClassic
import java.security.MessageDigest

/** Keep biological sample identity separate from independently processed libraries. */
class RnaseqUtils {
    static String digest(def value) {
        MessageDigest.getInstance('SHA-256').digest(JsonOutput.toJson(value).getBytes('UTF-8')).encodeHex().toString().take(16)
    }

    static Map readset(Map genome, def reads, String source, String accession = '') {
        def files = reads instanceof List ? reads : [reads]
        if (!(files.size() in [1, 2])) throw new IllegalArgumentException('Expected one FASTQ or an ordered R1/R2 pair')
        def identity = files.collect { path ->
            def f = new File(path.toString())
            [f.canonicalPath, f.length(), f.lastModified()]
        }
        def runId = accession ?: "${source}_${digest(identity.collect { it[0] })}"
        def signature = digest(identity)
        def key = "${runId}_${signature}"
        return genome + [id: "${genome.id}__${key}".toString(), genome_id: genome.id,
                         genome_meta: new LinkedHashMap(genome), readset_id: key.toString(),
                         run_id: runId.toString(), readset_signature: signature,
                         single_end: files.size() == 1]
    }

    static Map readStats(def path) {
        def report = new JsonSlurperClassic().parse(new File(path.toString()))
        def after = report.summary.after_filtering
        return [rnaseq_read_count: after.total_reads as long, rnaseq_total_bases: after.total_bases as long]
    }

    static long indexBytes(def directory) {
        long bytes = 0
        new File(directory.toString()).eachFileRecurse(FileType.FILES) { bytes += it.length() }
        return bytes
    }

    static int starMemoryGb(long indexBytes, long reads, def params) {
        // Initial scheduling estimate, not an exact STAR peak-memory prediction.
        double gib = indexBytes / (1024d * 1024 * 1024)
        double estimate = gib * (params.rnaseq_star_index_memory_factor as double) + 4d +
            reads / 1_000_000d * (params.rnaseq_star_memory_per_million_reads_gb as double)
        return Math.ceil(Math.max(params.rnaseq_star_min_memory_gb as double, estimate)) as int
    }

    static Map mergeMeta(def metas) {
        def genome = metas.first().genome_meta
        if (!metas.every { it.genome_meta == genome }) throw new IllegalArgumentException('Cannot merge libraries from different genomes')
        return genome + [genome_meta: genome,
                         rnaseq_run_signature: digest(metas.collect { it.readset_id }.sort())]
    }
}
