/** Explicitly validated imports from an older cache; never infer absent outputs. */
class LegacyCacheUtils {
    static Map load(def manifestPath, def cacheDir) {
        if (!manifestPath) return [:]
        def source = new File(manifestPath.toString())
        def manifest = new groovy.json.JsonSlurperClassic().parse(source)
        if (manifest.schema_version != 1 || !(manifest.genomes instanceof Map)) {
            throw new IllegalArgumentException("Unsupported cache import manifest: ${source}")
        }
        if (new File(manifest.cache_dir.toString()).canonicalPath != new File(cacheDir.toString()).canonicalPath) {
            throw new IllegalArgumentException('Cache import manifest belongs to a different cacheDir')
        }
        manifest.genomes.each { id, entry ->
            [entry.reference_identity, entry.mapping, entry.masked_source, entry.repeats_source].each { identity ->
                def original = new File(identity.path.toString())
                if (!original.isFile() || original.length() != identity.bytes ||
                    original.lastModified() != (identity.mtime_ns as BigInteger).intdiv(1000000L)) {
                    throw new IllegalArgumentException("Imported cache source changed; revalidate before reuse: ${original}")
                }
            }
            def paths = [entry.reference, entry.masked, entry.repeats_gff, entry.family_table, entry.families, entry.landscape] + entry.predictions.values()
            paths.each { path ->
                if (!new File(path.toString()).isFile()) {
                    throw new IllegalArgumentException("Missing imported cache artifact for ${id}: ${path}")
                }
            }
        }
        return manifest
    }

    static Map genome(Map manifest, String id) {
        if (!manifest) return [:]
        def entry = manifest.genomes[id]
        if (!entry) throw new IllegalArgumentException("Genome absent from cache import manifest: ${id}")
        return entry
    }

    static boolean hasPrediction(Map manifest, String id, String method) {
        return (genome(manifest, id).predictions ?: [:]).containsKey(method)
    }

    static boolean needsRnaseq(Map manifest, String id, String aligner) {
        if (!manifest) return true
        def enabled = aligner == 'both' ? ['star', 'varus'] : [aligner]
        return enabled.any { !hasPrediction(manifest, id, "braker-with-${it}".toString()) }
    }
}
