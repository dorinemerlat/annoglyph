process FILTER_PROTEIN_EVIDENCE {
    tag "${meta.id}"
    label 'annoglyph_python'
    // Filtering follows Nextflow's input cache rather than a shared storeDir.
    storeDir null

    input:
    tuple val(meta), path(proteins)
    path taxonomy
    path taxonomy_db

    output:
    tuple val(meta), path('filtered_protein_evidence.fa'), emit: proteins
    tuple val(meta), path('protein_evidence_filter.json'), emit: report

    script:
    """
    filter_protein_evidence.py \\
        --proteins ${proteins} \\
        --taxonomy ${taxonomy} \\
        --target-id '${meta.id}' \\
        --taxonomy-db ${taxonomy_db} \\
        --output filtered_protein_evidence.fa \\
        --report protein_evidence_filter.json
    """
}
