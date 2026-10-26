process RESOLVE_RFAM_METADATA {
    tag 'all_genomes'
    label 'annoglyph_python'
    storeDir null
    input:
    path rfam
    path tables, stageAs: 'infernal_tables/*'
    output:
    path 'rfam_resolved', emit: metadata
    path 'rfam_resolved/metadata_resolution.json', emit: audit
    script:
    def cacheRoot = file(params.cacheDir).toAbsolutePath()
    """
    python ${projectDir}/bin/resolve_rfam_metadata.py --families ${rfam}/family.txt \\
        --infernal ${tables.join(' ')} --cache-dir '${cacheRoot}/rfam_metadata/CURRENT' \\
        --output rfam_resolved
    """
}
