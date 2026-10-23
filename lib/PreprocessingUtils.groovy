/*
 * Utility methods used by the AnnoGlyph preprocessing workflow.
 */
class PreprocessingUtils {

    /*
     * Parse an NCBI taxonomy lineage file.
     */
    static Map parseLineage(lineageFile) {
        def lineage = [:]

        lineageFile.text.readLines().each { line ->
            line = line.trim()

            if (!line || line.startsWith('#')) {
                return
            }

            def columns = line.split('\t', -1)

            if (columns[0] == 'rank' || columns.size() < 3) {
                return
            }

            def rank = columns[0]
            def taxid = columns[1]
            def name = columns[2]

            lineage[rank] = [
                name: name,
                taxid: taxid
            ]
        }

        return lineage
    }


    /*
     * Read the first non-empty line from a file.
     */
    static String readFirstLine(inputFile) {
        return inputFile.text
            .readLines()
            .find { it.trim() }
            ?.trim()
    }


    /*
     * Generate a five-character mnemonic from a genome identifier.
     */
    static String generateFallbackMnemonic(String genomeId) {
        def nameParts = genomeId.tokenize('-')

        if (nameParts.size() >= 2) {
            return (
                nameParts[0].take(3) +
                nameParts[1].take(2)
            ).toUpperCase()
        }

        return genomeId
            .replaceAll(/[^A-Za-z0-9]/, '')
            .take(5)
            .toUpperCase()
    }


    /*
     * Add numeric suffixes when several genomes share a mnemonic.
     */
    static List makeUniqueMnemonics(List genomeIds, String mnemonic) {
        if (genomeIds.size() == 1) {
            return [[genomeIds[0], mnemonic]]
        }

        return genomeIds.withIndex().collect { genomeId, index ->
            [genomeId, "${mnemonic}${index + 1}"]
        }
    }


    /*
     * Check whether the user provided a mnemonic.
     */
    static boolean hasProvidedMnemonic(Map meta) {
        return meta.mnemonic &&
            meta.mnemonic.toString().trim() &&
            meta.mnemonic.toString() != 'NA'
    }


    /*
     * Select user-provided or automatically generated mnemonics.
     */
    static Map selectMnemonics(Map meta, String officialMnemonic, String annoglyphMnemonic) {
        if (hasProvidedMnemonic(meta)) {
            return [
                official: meta.mnemonic,
                annoglyph: meta.mnemonic
            ]
        }

        return [
            official: officialMnemonic,
            annoglyph: annoglyphMnemonic
        ]
    }


    /*
     * Extract the total assembly length from SeqKit stats.
     */
    static long parseAssemblySize(seqkitFile) {
        def lines = seqkitFile.text.readLines()

        if (lines.size() < 2) {
            throw new IllegalArgumentException(
                "SeqKit statistics file does not contain a data row: ${seqkitFile}"
            )
        }

        def columns = lines[1].split('\t', -1)

        if (columns.size() < 5) {
            throw new IllegalArgumentException(
                "SeqKit statistics file has an unexpected format: ${seqkitFile}"
            )
        }

        return columns[4].replace(',', '').trim().toLong()
    }
}