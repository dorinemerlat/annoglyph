process ANNOTATE_PROTEIN_FASTA {
    tag "${id}"
    label 'exogap_python'

    input:
    tuple val(id), val(meta), val(gff), path(proteins)

    output:
    tuple val(id), val(meta), path("${id}_with_functions.faa")

    script:
    """
    annotate_protein_fasta.py -f ${proteins} -g ${gff} -o ${id}_with_functions.faa -s "${meta.lineage.species.name}" -t "${meta.taxid}"
    """

    stub:
    """
    touch ${id}_with_functions.faa
    """
}