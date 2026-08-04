#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    validateParameters
    samplesheetToList
} from 'plugin/nf-schema'

// include { EXOGAP } from './workflows/exogap'

validateParameters()

def samplesheet = samplesheetToList(
    params.input,
    "$projectDir/assets/schema_input.json"
)

def genomes = SamplesheetUtils.buildGenomesChannel(
    samplesheet,
    log,
    { path, options = [:] ->
        file(path, checkIfExists: options.get('checkIfExists', false))
    }
)

LogUtils.printExecutionSummary(
    workflow,
    params,
    samplesheet,
    log
)

workflow {
    log.info 'Starting EXOGAP'

    // EXOGAP(genomes)
}