process SELECT_BUSCO_DATASET {
    label 'busco'

    input:
    path taxonomy

    output:
    path 'busco_datasets.tsv', emit: table

    script:
    """
    select_busco_dataset.py \
        --taxonomy ${taxonomy} \
        --output busco_datasets.tsv
    """

    stub:
    """
    command -v select_busco_dataset.py >/dev/null
    command -v busco >/dev/null

    cat > busco_datasets.tsv <<'EOF'
id	busco_dataset
test-species	eukaryota_odb12
EOF
    """
}