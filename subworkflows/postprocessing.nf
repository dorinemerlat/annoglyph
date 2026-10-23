include { EXPORT_ANNOTATION } from "$projectDir/modules/local/postprocessing/export_annotation.nf"

include { COMPARATIVE_REPORT } from "$projectDir/modules/local/postprocessing/comparative_report.nf"

workflow POSTPROCESSING {
    take:
    genomes
    coding_gff
    repeat_gff
    ncrna_gff
    species_assets
    comparative_assets
    taxonomy
    tree
    main:
    grouped_assets = species_assets.groupTuple(by: 0)
        .map { id, assets -> [id, assets.collect { it[0] }, assets.collect { it[1] }] }
    // Join on sample ID: selected annotations carry additional method/QC metadata.
    inputs = genomes.map { meta, genome -> [meta.id, meta, genome] }
        .join(coding_gff.map { meta, gff -> [meta.id, meta, gff] }, remainder: true, failOnDuplicate: true)
        .filter { row -> row[1] != null }
        .map { row ->
            def (id, meta, genome, coding_meta, coding) = row
            [id, meta + (coding_meta ?: [:]), genome, coding ?: []]
        }
        .join(repeat_gff.map { meta, gff -> [meta.id, gff] }, remainder: true, failOnDuplicate: true)
        .filter { row -> row[1] != null }
        .map { row ->
            def (id, meta, genome, coding, repeats) = row
            [id, meta, genome, coding, repeats ?: []]
        }
        .join(ncrna_gff.map { meta, gff -> [meta.id, gff] }, remainder: true, failOnDuplicate: true)
        .filter { row -> row[1] != null }
        .map { row ->
            def (id, meta, genome, coding, repeats, ncrna) = row
            [id, meta, genome, coding, repeats, ncrna ?: []]
        }
        .join(grouped_assets, remainder: true, failOnDuplicate: true)
        .filter { row -> row[1] != null }
        .map { row ->
            def (id, meta, genome, coding, repeats, ncrna, assets, destinations) = row
            [meta, genome, coding, repeats, ncrna, assets ?: [], destinations ?: []]
        }
    if (params.validate_submission && (!params.submission_metadata || !params.submission_template || !params.table2asn_binary || !params.table2asn_sha256)) {
        error 'Submission validation requires submission_metadata, submission_template, table2asn_binary and table2asn_sha256'
    }
    EXPORT_ANNOTATION(inputs,
        params.submission_metadata ? Channel.value(file(params.submission_metadata, checkIfExists: true)) : [],
        params.submission_template ? Channel.value(file(params.submission_template, checkIfExists: true)) : [],
        params.table2asn_binary ? Channel.value(file(params.table2asn_binary, checkIfExists: true)) : [])
    all_comparative_assets = comparative_assets.collect(flat: false)
    COMPARATIVE_REPORT(EXPORT_ANNOTATION.out.bundle.map { meta, folder -> folder }.collect(),
        taxonomy.first(), tree.first(), file("${projectDir}/assets/results_README.md", checkIfExists: true),
        all_comparative_assets.map { assets -> assets.collect { it[0] } },
        all_comparative_assets.map { assets -> assets.collect { it[1] } })
    emit:
    comparative = COMPARATIVE_REPORT.out.comparative
    bundles = EXPORT_ANNOTATION.out.bundle
}
