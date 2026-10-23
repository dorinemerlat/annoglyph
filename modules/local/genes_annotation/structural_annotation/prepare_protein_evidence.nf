process PREPARE_PROTEIN_EVIDENCE {
    label 'annoglyph_python'

    input:
    path taxonomy
    val busco_dataset

    output:
    path 'protein_evidence.fa', emit: proteins
    path '*.release.json', emit: release, optional: true
    path '*_odb*.fa', emit: orthodb
    path '*_uniprot.fa', emit: uniprot
    path 'protein_evidence.tsv', emit: metadata
    path 'prepare_protein_evidence.log', emit: log

    script:
    def evaluationArg = params.evaluation_mode ? '--evaluation' : ''
    def releaseArg = params.uniprot_evidence_release ? "--uniprot-release ${params.uniprot_evidence_release}" : ''
    def minSequences = params.protein_evidence_min_sequences ?: 100000
    def maxSequences = params.protein_evidence_max_sequences ?: 0
    def buscoVersionMatch = busco_dataset =~ /_odb(\d+)/
    def orthodbVersion = buscoVersionMatch ? buscoVersionMatch[0][1] : null

    if (!orthodbVersion) {
        error "Cannot determine OrthoDB version from BUSCO dataset: ${busco_dataset}"
    }

    """
    prepare_protein_evidence.py \
        ${evaluationArg} \
        --taxonomy ${taxonomy} \
        --orthodb-version ${orthodbVersion} \
        --min-sequences ${minSequences} \
        --max-sequences ${maxSequences} \
        ${releaseArg} \
        --output protein_evidence.fa \
        --metadata-output protein_evidence.tsv \
        > prepare_protein_evidence.log 2>&1
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
source	clade	rank	taxid	version	sequences	file
OrthoDB	Eukaryota			12	1	test_odb12.fa
UniProtKB	Eukaryota	superkingdom	2759		1	test_uniprot.fa
combined					2	protein_evidence.fa
EOF

    cat > prepare_protein_evidence.log <<'EOF'
Selected UniProt evidence clade: Eukaryota
Selected OrthoDB partition: Eukaryota
Final protein evidence set: 2 sequences
EOF
    """
}