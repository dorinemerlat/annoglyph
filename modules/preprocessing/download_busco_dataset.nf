process DOWNLOAD_BUSCO_DATASET {
    tag "${busco_dataset}"
    label 'busco'

    input:
    val busco_dataset

    output:
    tuple val(busco_dataset), path("lineages/${busco_dataset}"), emit: dataset

    script:
    """
    busco \
        --download ${busco_dataset} \
        --download_path .
    """

    stub:
    """
    command -v busco >/dev/null
    mkdir -p lineages/${busco_dataset}
    """
}