process GET_TAXONOMY {
    tag "${id}"
    label 'exogap_python'

    input:
    tuple val(id), val(meta), val(genome)

    output:
    tuple val(id), val(meta), path("${meta.taxid}_lineage.csv"), emit: lineage
    tuple val(id), val(meta), path("${meta.taxid}.mnemonic"), emit: mnemonic

    script:
    """
    # download lineage and mnemonic
    fetch_taxonomy.py $meta.taxid
    """

    stub:
    """
    touch ${meta.taxid}_lineage.csv ${meta.taxid}.mnemonic
    """
}
