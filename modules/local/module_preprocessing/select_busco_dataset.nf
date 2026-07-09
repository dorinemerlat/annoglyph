process SELECT_BUSCO_DATASET {
    tag "$id"
    label 'busco'

    input:
    tuple  val(id), val(meta), path(lineage)

    output:
    tuple val(id), val(meta), path("datasets.txt")

    script:
    """
    select_busco_datasets.sh $lineage > datasets.txt
    """

    stub:
    """
    touch datasets.txt
    """
}
