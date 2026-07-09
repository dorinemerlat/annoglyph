#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/exogap
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/exogap
    Website: https://nf-co.re/exogap
    Slack  : https://nfcore.slack.com/channels/exogap
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// Print help message if needed
if (params.help) {
    def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
    def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
    def String command = "nextflow run ${workflow.manifest.name} -c nextflow.config -profile singularity"
    log.info logo + paramsHelp(command) + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)
    System.exit(0)
}

// // Validate input parameters
// if (params.validate_params) {
//     validateParameters()
// }
// WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check is the user provided a samplesheet
if (params.input) {
    // Check if the samplesheet exists
    def samplesheet_file = file(params.input, checkIfExists: true)
    // Convert samplesheet to list of maps
    samplesheet = samplesheetToList(samplesheet_file, "$projectDir/assets/schema_input.json")
    log.info "Found ${samplesheet.size()} genomes in the samplesheet: ${params.input}"
} else {
    log.error "No samplesheet provided. You must provide a samplesheet with --input or in nextflow.config file to run the workflow."
    System.exit(1)
}

// check if at least one module is activated
if (params.module_repeats == false && params.module_genes == false && params.module_ncgenes == false) {
    log.error "You have deactivated all annotation modules. Please activate at least one module to run the workflow."
    System.exit(1)
}

// if module_repeats is activated, check its parameters
if (params.module_repeats == true) {
    // check group_consensus_sequences
    if (params.group_consensus_sequences != true && params.group_consensus_sequences != false) {
        log.error "Parameter error: --group_consensus_sequences must be true or false."
        System.exit(1)
    }

    // Check reference_library parameters
    if (params.reference_library) {
        def ref_lib_file = file(params.reference_library)
        if (!ref_lib_file.exists()) {
            log.warn "The file provided for --reference_library does not exist: ${params.reference_library}"
        } else {
            params.reference_library = ref_lib_file
            log.info "Reference library file provided for repeat annotation: ${params.reference_library}"
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EXTRACT GENOME INFORMATION FROM SAMPLESHEET
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

Channel.fromList(samplesheet)
    .map { unformated_meta, fasta_path ->
        // format meta
        def meta = unformated_meta.collectEntries { key, value ->
            if (value instanceof List && value.isEmpty()) {
                [(key): null]
            } else {
                [(key): value]
            }
        }
        
        // create a unique ID for each genome based on its name, in lowercase and with spaces replaced by '-'
        def genome_id = meta.name.toLowerCase().replace(' ', '-') 
        
        // format fasta file
        def fasta = file(fasta_path, checkIfExists: true)

        // Make sure taxid is a string
        meta = Utils.updateTopLevelKey(
            meta,
            'taxid',
            meta.taxid.toString()
        )

        // RNA-seq directory from samplesheet column RNASeq-dir
        def rnaseq_files = [
            single_end: [],
            paired_end: []
        ]

        def rnaseq_dir_path = meta['RNASeq-dir']

        if (rnaseq_dir_path) {
            def rnaseq_dir = file(rnaseq_dir_path, checkIfExists: true)

            if (!rnaseq_dir.isDirectory()) {
                log.error "RNASeq-dir for sample '${meta.name}' must be a directory: ${rnaseq_dir_path}"
                System.exit(1)
            }

            def allowed_extensions_regex = ~/(?i)^.*\.(?:fq|fastq)(?:\.gz)?$/

            def read_files = (new File(rnaseq_dir.toString()).listFiles() ?: [])
                .findAll { it.isFile() && it.name ==~ allowed_extensions_regex }
                .sort { a, b -> a.name <=> b.name }

            def paired_files = [:].withDefault { [r1: null, r2: null] }
            def single_files = [:]

            read_files.each { read_file ->
                def read_stem = read_file.name.replaceFirst(/(?i)\.(?:fq|fastq)(?:\.gz)?$/, '')
                def read_match = read_stem =~ /(?i)^(.*)_r?([12])$/

                if (read_match.matches()) {
                    def run_id = read_match[0][1]
                    def mate_id = read_match[0][2]

                    if (mate_id == '1') {
                        paired_files[run_id].r1 = file(read_file)
                    } else {
                        paired_files[run_id].r2 = file(read_file)
                    }
                } else {
                    if (single_files.containsKey(read_stem)) {
                        log.warn "Duplicate single-end RNA-seq run id '${read_stem}' for sample '${meta.name}'. Keeping latest file: ${read_file}"
                    }
                    single_files[read_stem] = file(read_file)
                }
            }

            paired_files.keySet().sort().each { run_id ->
                def mates = paired_files[run_id]

                if (mates.r1 && mates.r2) {
                    rnaseq_files.paired_end << [mates.r1, mates.r2]
                } else {
                    log.warn "Incomplete paired-end RNA-seq run for sample '${meta.name}', run '${run_id}' (r1=${mates.r1}, r2=${mates.r2})"
                }
            }

            single_files.keySet().sort().each { run_id ->
                rnaseq_files.single_end << single_files[run_id]
            }
        }

        // update meta field for RNASeq-file
        if (rnaseq_files.single_end || rnaseq_files.paired_end) {
            meta['RNASeq-file'] = rnaseq_files
        } else {
            meta['RNASeq-file'] = null
        }

        meta.remove('RNASeq-dir') 

        // SRA id from samplesheet column SRA
        def sra_ids = []
        if (meta['SRA']) {
            if (meta['SRA'] instanceof List) {
                sra_ids = meta['SRA'].findAll { it }
            } else {
                sra_ids = meta['SRA']
                    .toString()
                    .split(/[;, \t]+/)
                    .findAll { it }
            }
        }

        // update meta field for SRA
        meta['SRA'] = sra_ids ?: null


        [ genome_id, meta, fasta ] 
    }
    .set { genomes }
    

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXOGAP } from './workflows/exogap'


// WORKFLOW: Run main nf-core/exogap analysis pipeline

workflow {
    log.info "Starting workflow EXOGAP"

    EXOGAP(genomes)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
