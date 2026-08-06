process PREPARE_PROTEIN_EVIDENCE {
    label 'exogap_python'

    input:
    path taxonomy
    path busco_datasets

    output:
    path 'protein_evidence.fa', emit: proteins
    path '*_odb*.fa', emit: orthodb
    path '*_uniprot.fa', emit: uniprot
    path 'protein_evidence.tsv', emit: metadata

    script:
    def maxSequences = params.protein_evidence_max_sequences ?: 0

    """
    prepare_protein_evidence.py \
        --taxonomy ${taxonomy} \
        --busco-datasets ${busco_datasets} \
        --max-sequences ${maxSequences} \
        --output protein_evidence.fa \
        --metadata-output protein_evidence.tsv
    """

    stub:
    """
    command -v prepare_protein_evidence.py >/dev/null

    cat > protein_evidence.fa <<'EOF'
>test_protein
MAAA
EOF

    cp protein_evidence.fa test_odb12.fa
    cp protein_evidence.fa test_uniprot.fa

    cat > protein_evidence.tsv <<'EOF'
source	clade	taxid	version	busco_datasets	sequences	file
OrthoDB	Eukaryota	2759	12	eukaryota_odb12	1	test_odb12.fa
UniProtKB	Eukaryota	2759		eukaryota_odb12	1	test_uniprot.fa
combined	Eukaryota	2759		eukaryota_odb12	2	protein_evidence.fa
EOF
    """
}