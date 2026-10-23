/*
 * Utility methods for AnnoGlyph logging.
 */
class LogUtils {

    /*
     * Print the execution parameters before starting the workflow.
     */
    static void printExecutionSummary(workflow, params, samplesheet, log) {
        def genomeCount = samplesheet.size()

        def genomesWithRnaseq = samplesheet.count { entry ->
            def meta = entry[0]
            meta['RNASeq-dir']
        }

        def genomesWithSra = samplesheet.count { entry ->
            def meta = entry[0]
            meta.SRA
        }

        log.info """
        ============================================================
                             AnnoGlyph
        ============================================================

        Core Nextflow options
          Run name           : ${workflow.runName}
          Run group          : ${params.run_group ?: "none"}
          Profile            : ${workflow.profile}
          Container engine   : ${workflow.containerEngine}
          Launch directory   : ${workflow.launchDir}
          Work directory     : ${workflow.workDir}
          Project directory  : ${workflow.projectDir}

        Input/output
          Samplesheet        : ${params.input}
          Genomes            : ${genomeCount}
          Genomes with RNA-seq: ${genomesWithRnaseq}
          Genomes with SRA   : ${genomesWithSra}
          Output directory   : ${params.outdir}
          Publish mode       : ${params.publish_dir_mode}

        Annotation modules
          Repeats            : ${formatBoolean(params.module_repeats)}
          Protein-coding genes: ${formatBoolean(params.module_genes)}
          Non-coding genes   : ${formatBoolean(params.module_ncgenes)}
          Comparative genomics: ${formatBoolean(params.comparative_genomics)}

        Shared annotation
          BUSCO lineage      : ${params.busco_lineage ?: 'automatic'}

        Repeat annotation
          Reference library  : ${formatPath(params.reference_library)}
          Group consensus    : ${formatBoolean(params.group_consensus_sequences)}

        Protein-coding gene annotation
          Protein set        : ${formatPath(params.protein_set)}
          BLAST database     : ${formatPath(params.blast_db)}
          RNA-seq aligner    : ${params.rna_seq_aligner}
          Protein-only tool  : ${params.tool_only_protein_evidence}
          OMArk threshold    : ${params.omark_contamination_threshold}

        ============================================================
        """.stripIndent()
    }


    /*
     * Print the workflow completion status.
     */
    static void printCompletionSummary(workflow, log) {
        if (workflow.success) {
            log.info "[AnnoGlyph] Pipeline completed successfully."
        } else {
            log.error "[AnnoGlyph] Pipeline completed with errors."
        }
    }


    /*
     * Format boolean values for display.
     */
    static String formatBoolean(value) {
        return value ? 'yes' : 'no'
    }


    /*
     * Return the final component of a path.
     */
    static String formatPath(value) {
        return value ? new File(value.toString()).name : 'none'
    }
}