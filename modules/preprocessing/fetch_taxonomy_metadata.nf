process FETCH_TAXONOMY_METADATA {
    tag "${id}"
    label 'exogap_python'

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path("${meta.taxid}_lineage.tsv"), emit: lineage
    tuple val(id), val(meta), path("${meta.taxid}.mnemonic"), emit: mnemonic

    script:
    """
    fetch_taxonomy_metadata.py ${meta.taxid}
    """

    stub:
    """
    # Check that the executable is available
    command -v fetch_taxonomy_metadata.py >/dev/null

    # Minimal lineage table
    cat > ${meta.taxid}_lineage.tsv << EOF
    rank	taxid	name
    EOF

    # Placeholder mnemonic
    echo "NA" > ${meta.taxid}.mnemonic
    """
}