/*
 * Utility methods used by the repeat annotation workflow.
 */
class RepeatAnnotationUtils {

    /*
     * Group per-genome summary tables and define the merged output name.
     */
    static def gatherSummaryTables(channel, String prefix) {
        return channel
            .map { id, meta, table ->
                [prefix, table]
            }
            .groupTuple(by: 0)
            .map { group, tables ->
                [group, tables, "${prefix}.tsv"]
            }
    }
}