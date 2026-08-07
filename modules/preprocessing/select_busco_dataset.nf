process SELECT_BUSCO_DATASET {
    label 'busco'

    input:
    path taxonomy

    output:
    path 'busco_dataset.txt'

    script:
    """
    select_busco_dataset.py \
        --taxonomy ${taxonomy} \
        --output busco_dataset.txt
    """

    stub:
    """
    command -v select_busco_dataset.py >/dev/null
    command -v busco >/dev/null

    cat > busco_dataset.txt <<'EOF'
eukaryota_odb12
EOF
    """
}