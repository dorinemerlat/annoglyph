process COMPARATIVE_REPORT {
    tag 'results_layout_v2'
    label 'annoglyph_python'
    memory 8.GB
    publishDir "${params.outdir}", mode: params.publish_dir_mode, overwrite: true
    input:
    path bundles, stageAs: 'bundles/*'
    path taxonomy
    path tree
    path readme
    path assets, stageAs: 'assets??/*'
    val destinations
    output:
    path 'comparative', emit: comparative
    path 'README.md', emit: readme
    script:
    def assetFiles = assets instanceof List ? assets : (assets ? [assets] : [])
    def assetManifest = groovy.json.JsonOutput.toJson(assetFiles.withIndex().collect { path, index ->
        [source: path.toString(), destination: destinations[index]]
    })
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib
    cat > assets.json <<'ASSET_MANIFEST'
    ${assetManifest}
    ASSET_MANIFEST
    report_comparative.py --bundles bundles --taxonomy ${taxonomy} --tree ${tree} --output comparative --assets assets.json
    cp ${readme} README.md
    """
    stub:
    """
    mkdir comparative
    cp ${readme} README.md
    """
}
