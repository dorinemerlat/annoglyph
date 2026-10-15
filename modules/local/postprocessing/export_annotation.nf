process EXPORT_ANNOTATION {
    tag "layout_v2/${meta.id}"
    label 'exogap_python'
    memory 8.GB
    publishDir "${params.outdir}/species", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(genome), path(coding), path(repeats), path(ncrna), path(assets, stageAs: 'assets??/*'), val(destinations)

    path submission_metadata
    path submission_template
    path validator

    output:
    tuple val(meta), path("${meta.id}"), emit: bundle

    script:
    def assetFiles = assets instanceof List ? assets : (assets ? [assets] : [])
    def assetManifest = groovy.json.JsonOutput.toJson(assetFiles.withIndex().collect { path, index ->
        [source: path.toString(), destination: destinations[index]]
    })
    def codingArg = coding ? "--coding ${coding}" : ''
    def repeatsArg = repeats ? "--repeats ${repeats}" : ''
    def ncrnaArg = ncrna ? "--ncrna ${ncrna}" : ''
    def locusPrefix = meta.locus_tag_prefix ?: params.locus_tag_prefix ?: meta.mnemonic.exogap
    def organism = meta.lineage?.species?.name ?: meta.name
    def code = meta.translation_table ?: params.translation_table
    def metadataArg = submission_metadata ? "--submission-metadata ${submission_metadata}" : ''
    def validationCommand = params.validate_submission ? "validate_submission.py --bundle ${meta.id} --template ${submission_template} --binary ${validator} --sha256 ${params.table2asn_sha256} --fail-on-errors ${params.submission_fail_on_errors}" : ''
    def selectedMethod = meta.annotation_method ?: 'none'
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib
    cat > assets.json <<'ASSET_MANIFEST'
    ${assetManifest}
    ASSET_MANIFEST
    export_annotation.py --genome ${genome} ${codingArg} ${repeatsArg} ${ncrnaArg} \
        --sample '${meta.id}' --organism '${organism}' --taxid '${meta.taxid}' \
        --method '${selectedMethod}' --locus-prefix '${locusPrefix}' \
        --translation-table ${code} ${metadataArg} --assets assets.json --output ${meta.id}
    report_species.py --bundle ${meta.id} --assets assets.json
    ${validationCommand}
    checksum_bundle.py --bundle ${meta.id}
    """

    stub:
    """
    mkdir -p ${meta.id}/genome/sequences ${meta.id}/genome/stats ${meta.id}/genome/plots
    mkdir -p ${meta.id}/protein_coding/sequences ${meta.id}/protein_coding/selection ${meta.id}/protein_coding/stats ${meta.id}/protein_coding/plots
    mkdir -p ${meta.id}/noncoding/sequences ${meta.id}/noncoding/stats ${meta.id}/noncoding/plots
    mkdir -p ${meta.id}/repeats/sequences ${meta.id}/repeats/stats ${meta.id}/repeats/plots ${meta.id}/submission
    echo '##gff-version 3' > ${meta.id}/genome/annotation.gff3
    echo 'genome_id\torganism\tgenome_size_bp\tcontigs\tN50_bp\tL50\tGC_percent\tN_bp\tsoftmasked_bp\trepeat_count\trepeat_covered_bp\trepeat_coverage_percent\tprotein_count' > ${meta.id}/genome/stats/overview.tsv
    touch ${meta.id}/genome/sequences/genome.softmasked.fna ${meta.id}/genome/sequences/genome.fna
    touch ${meta.id}/protein_coding/sequences/cds.fna ${meta.id}/protein_coding/sequences/proteins.faa ${meta.id}/noncoding/sequences/ncrna.fna
    """
}
