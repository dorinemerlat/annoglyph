process GENERATE_TAXONOMY_METADATA {
    label 'annoglyph_python'

    input:
    path taxonomy_input

    output:
    path 'taxonomy.tsv', emit: table
    path 'taxonomy.json', emit: json
    path 'tree.nwk', emit: newick
    path 'tree_minimal.nwk', emit: minimal_tree
    
    script:
    """
    generate_taxonomy_metadata.py --input ${taxonomy_input} --output .
    """

    stub:
    """
    command -v generate_taxonomy_metadata.py >/dev/null

    cat > taxonomy.tsv <<'EOF'
id	name	scientific_name	taxid	official_mnemonic	annoglyph_mnemonic	superkingdom	kingdom	species
test-species	Test species	Test species	12345	TEST	TEST	Eukaryota	Metazoa	Test species
EOF

    cat > taxonomy.json <<'EOF'
{
  "test-species": {
    "scientific_name": "Test species",
    "taxid": "12345",
    "mnemonic": {
      "official": "TEST",
      "annoglyph": "TEST"
    },
    "lineage": {
      "superkingdom": {
        "name": "Eukaryota",
        "taxid": "2759"
      },
      "kingdom": {
        "name": "Metazoa",
        "taxid": "33208"
      },
      "species": {
        "name": "Test species",
        "taxid": "12345"
      }
    }
  }
}
EOF

  cat > tree.nwk <<'EOF'
(Test species);
EOF
    """
}