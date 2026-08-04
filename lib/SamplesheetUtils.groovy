import nextflow.Channel


/*
 * Utility methods for formatting the EXOGAP input samplesheet.
 */
class SamplesheetUtils {

    /*
     * Build the genome input channel from validated samplesheet entries.
     *
     * Output:
     *   [genome_id, meta, fasta]
     */
    static def buildGenomesChannel(samplesheet, log, fileFunc) {
        return Channel
            .fromList(samplesheet)
            .map { unformattedMeta, fastaPath ->
                def meta = normalizeMetadata(unformattedMeta)
                def genomeId = formatGenomeId(meta.name)
                def fasta = fileFunc(fastaPath, [checkIfExists: true])

                meta.taxid = meta.taxid.toString()

                parseRnaseqDirectory(meta, log, fileFunc)
                parseSraAccessions(meta)

                [genomeId, meta, fasta]
            }
    }


    /*
     * Replace empty values with null.
     */
    static Map normalizeMetadata(Map metadata) {
        return metadata.collectEntries { key, value ->
            def normalizedValue = value instanceof List && value.isEmpty()
                ? null
                : value

            [(key): normalizedValue]
        }
    }


    /*
     * Convert a genome name into a lowercase identifier.
     */
    static String formatGenomeId(String name) {
        return name
            .toLowerCase()
            .replaceAll(/\s+/, '-')
    }


    /*
     * Parse RNA-seq files from the optional RNASeq-dir column.
     */
    static void parseRnaseqDirectory(Map meta, log, fileFunc) {
        def rnaseqFiles = [
            single_end: [],
            paired_end: []
        ]

        if (meta['RNASeq-dir']) {
            def rnaseqDir = fileFunc(meta['RNASeq-dir'], [checkIfExists: true])

            def readFiles = (new File(rnaseqDir.toString()).listFiles() ?: [])
                .findAll { readFile ->
                    readFile.isFile() &&
                        readFile.name ==~ /(?i)^.*\.(?:fq|fastq)(?:\.gz)?$/
                }
                .sort { first, second ->
                    first.name <=> second.name
                }

            def pairedFiles = [:].withDefault {
                [r1: null, r2: null]
            }

            def singleFiles = [:]

            readFiles.each { readFile ->
                def readStem = readFile.name.replaceFirst(
                    /(?i)\.(?:fq|fastq)(?:\.gz)?$/,
                    ''
                )

                def readMatch = readStem =~ /(?i)^(.*)_r?([12])$/

                if (readMatch.matches()) {
                    def runId = readMatch[0][1]
                    def mateId = readMatch[0][2]
                    def mateKey = mateId == '1' ? 'r1' : 'r2'

                    pairedFiles[runId][mateKey] = fileFunc(readFile.toString())
                } else {
                    singleFiles[readStem] = fileFunc(readFile.toString())
                }
            }

            pairedFiles.keySet().sort().each { runId ->
                def mates = pairedFiles[runId]

                if (mates.r1 && mates.r2) {
                    rnaseqFiles.paired_end << [mates.r1, mates.r2]
                } else {
                    log.warn "Incomplete paired-end RNA-seq run for '${meta.name}': ${runId}"
                }
            }

            singleFiles.keySet().sort().each { runId ->
                rnaseqFiles.single_end << singleFiles[runId]
            }
        }

        meta['RNASeq-file'] = rnaseqFiles.single_end || rnaseqFiles.paired_end
            ? rnaseqFiles
            : null

        meta.remove('RNASeq-dir')
    }


    /*
     * Parse optional SRA accessions.
     */
    static void parseSraAccessions(Map meta) {
        def sraIds = meta.SRA instanceof List
            ? meta.SRA.findAll()
            : meta.SRA
                ?.toString()
                ?.split(/[;, \t]+/)
                ?.findAll()

        meta.SRA = sraIds ?: null
    }
}