import nextflow.Channel

/*
 * Utility methods for parsing the EXOGAP input samplesheet.
 *
 * This class replaces nf-schema/samplesheetToList for the samplesheet parsing step.
 * The goal is to avoid plugin-dependent parsing issues and keep the input logic
 * explicit, maintainable and compatible with recent Nextflow versions.
 *
 * Expected CSV columns:
 *   - name        : species/genome name
 *   - taxid       : NCBI taxonomic identifier
 *   - fasta       : genome FASTA path
 *   - RNASeq-dir  : optional directory containing RNA-seq FASTQ files
 *   - SRA         : optional SRA/ENA/DDBJ accession(s), separated by ';', ',' or spaces
 *
 * Output channel:
 *   tuple(genome_id, meta, fasta)
 */
class SamplesheetUtils {

    static def buildGenomesChannel(params, log, fileFunc) {

        if (!params.input) {
            log.error "No samplesheet provided. You must provide a samplesheet with --input or in nextflow.config."
            System.exit(1)
        }

        def samplesheetFile = fileFunc(params.input, [checkIfExists: true])

        return Channel
            .fromPath(samplesheetFile)
            .splitCsv(header: true, sep: ',')
            .map { row ->

                def meta = row.collectEntries { k, v ->
                    [(k): (v == null || v.toString().trim() == '' ? null : v)]
                }

                SamplesheetUtils.checkRequiredFields(meta, row, log)

                def fasta = fileFunc(meta.fasta, [checkIfExists: true])
                def genomeId = SamplesheetUtils.formatGenomeId(meta.name)

                meta.taxid = meta.taxid.toString()

                SamplesheetUtils.parseRnaseqDirectory(meta, log, fileFunc)
                SamplesheetUtils.parseSraAccessions(meta)

                return [genomeId, meta, fasta]
            }
    }

    /*
     * Check mandatory columns before the workflow starts.
     */
    static void checkRequiredFields(meta, row, log) {
        if (!meta.name) {
            log.error "Missing name in samplesheet row: ${row}"
            System.exit(1)
        }

        if (!meta.taxid) {
            log.error "Missing taxid for sample: ${meta.name}"
            System.exit(1)
        }

        if (!meta.fasta) {
            log.error "Missing fasta for sample: ${meta.name}"
            System.exit(1)
        }
    }

    /*
     * Convert a species/genome name into a stable lowercase identifier.
     */
    static String formatGenomeId(String name) {
        return name
            .toLowerCase()
            .replaceAll(/\s+/, '-')
    }

    /*
     * Parse an optional RNASeq-dir column.
     *
     * FASTQ files are classified as:
     *   - paired-end if filenames end with _1/_2, _R1/_R2, _r1/_r2
     *   - single-end otherwise
     *
     * The resulting structure is stored in:
     *   meta['RNASeq-file'] = [
     *       single_end: [...],
     *       paired_end: [[r1, r2], ...]
     *   ]
     *
     * The original RNASeq-dir key is removed afterwards because downstream
     * modules expect RNASeq-file, not RNASeq-dir.
     */
    static void parseRnaseqDirectory(meta, log, fileFunc) {

        def rnaseqFiles = [
            single_end: [],
            paired_end: []
        ]

        def rnaseqDirPath = meta['RNASeq-dir']

        if (rnaseqDirPath) {
            def rnaseqDir = fileFunc(rnaseqDirPath, [checkIfExists: true])

            if (!rnaseqDir.isDirectory()) {
                log.error "RNASeq-dir for sample '${meta.name}' must be a directory: ${rnaseqDirPath}"
                System.exit(1)
            }

            def allowedExtensionsRegex = ~/(?i)^.*\.(?:fq|fastq)(?:\.gz)?$/

            def readFiles = (new File(rnaseqDir.toString()).listFiles() ?: [])
                .findAll { it.isFile() && it.name ==~ allowedExtensionsRegex }
                .sort { a, b -> a.name <=> b.name }

            def pairedFiles = [:].withDefault { [r1: null, r2: null] }
            def singleFiles = [:]

            readFiles.each { readFile ->
                def readStem = readFile.name.replaceFirst(/(?i)\.(?:fq|fastq)(?:\.gz)?$/, '')
                def readMatch = readStem =~ /(?i)^(.*)_r?([12])$/

                if (readMatch.matches()) {
                    def runId = readMatch[0][1]
                    def mateId = readMatch[0][2]

                    if (mateId == '1') {
                        pairedFiles[runId].r1 = fileFunc(readFile.toString())
                    } else {
                        pairedFiles[runId].r2 = fileFunc(readFile.toString())
                    }
                } else {
                    singleFiles[readStem] = fileFunc(readFile.toString())
                }
            }

            pairedFiles.keySet().sort().each { runId ->
                def mates = pairedFiles[runId]

                if (mates.r1 && mates.r2) {
                    rnaseqFiles.paired_end << [mates.r1, mates.r2]
                } else {
                    log.warn "Incomplete paired-end RNA-seq run for sample '${meta.name}', run '${runId}'"
                }
            }

            singleFiles.keySet().sort().each { runId ->
                rnaseqFiles.single_end << singleFiles[runId]
            }
        }

        meta['RNASeq-file'] = (rnaseqFiles.single_end || rnaseqFiles.paired_end) ? rnaseqFiles : null
        meta.remove('RNASeq-dir')
    }

    /*
     * Parse optional SRA accessions.
     *
     * Accepts one or several accessions separated by:
     *   - semicolon
     *   - comma
     *   - spaces
     *   - tabs
     */
    static void parseSraAccessions(meta) {
        def sraIds = []

        if (meta['SRA']) {
            sraIds = meta['SRA']
                .toString()
                .split(/[;, \t]+/)
                .findAll { it }
        }

        meta['SRA'] = sraIds ?: null
    }
}