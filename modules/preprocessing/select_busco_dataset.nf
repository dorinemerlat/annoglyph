process SELECT_BUSCO_DATASET {
    tag "${id}"
    label 'busco'

    input:
    tuple val(id), val(meta), path(lineage)

    output:
    tuple val(id), val(meta), path("dataset.txt")

    script:
    """
    select_busco_dataset.py ${lineage} > dataset.txt
    """

    stub:
    """
    command -v select_busco_dataset.py >/dev/null
    command -v busco >/dev/null

    echo "eukaryota_odb10" > dataset.txt
    """
}