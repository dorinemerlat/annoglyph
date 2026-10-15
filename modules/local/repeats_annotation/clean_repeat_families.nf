process CLEAN_REPEAT_FAMILIES {
    tag "${meta.id}"
    label 'annoglyph_python'

    input:
    tuple val(meta), path(families)
    each path(classification)

    output:
    tuple val(meta), path("${families.baseName}.cleaned_for_mchelper.fa"), emit: for_mchelper
    tuple val(meta), path("${families.baseName}.tsv"), emit: table

    script:
    """
    clean_repeat_families.py \
        --fasta ${families} \
        --classification ${classification} \
        --mnemonic ${meta.mnemonic.annoglyph}
    """

    stub:
    """
    command -v clean_repeat_families.py >/dev/null

    cat > ${families.baseName}.cleaned_for_mchelper.fa << EOF
>${meta.mnemonic.annoglyph}_TE1#Unknown
ACGTACGTACGTACGTACGTACGTACGTACGT
EOF

    cat > ${families.baseName}.tsv << EOF
id	rm2_id	rm2_length	rm_code	rm2_class	rm2_order	rm2_superfamily	rm2_family	rm2_subfamily	rm2_subfamily2
${meta.mnemonic.annoglyph}_TE1	stub_family	32	Unknown						
EOF
    """
}