process NORMALIZE_INFERNAL_OUTPUTS {
    tag "${meta.id}"
    label 'annoglyph_python'
    input:
    tuple val(meta), path(table), path(report)
    output:
    tuple val(meta), path("${meta.id}_infernal.tblout"), emit: tblout
    tuple val(meta), path("${meta.id}_infernal.out"), emit: report
    script:
    """
    gzip -cd ${table} > ${meta.id}_infernal.tblout
    gzip -cd ${report} > ${meta.id}_infernal.out
    """
    stub:
    """
    echo '# cmscan fmt 2' > ${meta.id}_infernal.tblout
    touch ${meta.id}_infernal.out
    """
}
