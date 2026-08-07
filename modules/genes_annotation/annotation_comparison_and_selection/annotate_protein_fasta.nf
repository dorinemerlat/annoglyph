process ANNOTATE_PROTEIN_FASTA {
    tag "${id}"
    label 'exogap_python'

    input:
    tuple val(id), val(meta), val(gff), path(proteins)

    output:
    tuple val(id), val(meta), path("${id}_with_functions.faa")

    script:
    """
    annotate_protein_fasta.py \
        --fasta ${proteins} \
        --gff ${gff} \
        --output ${id}_with_functions.faa \
        --species "${meta.lineage.species.name}" \
        --taxid "${meta.taxid}"
    """

    stub:
    """
    command -v annotate_protein_fasta.py >/dev/null

    touch ${id}_with_functions.faa
    """
}