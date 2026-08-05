process PREPARE_BUSCO_HMM {
    tag "${busco_dataset}"
    label 'mchelper'

    input:
    tuple val(busco_dataset), path(busco_dataset_dir)

    output:
    tuple val(busco_dataset), path("${busco_dataset}.hmm*")

    script:
    """
    cat ${busco_dataset_dir}/hmms/*.hmm > ${busco_dataset}.hmm
    hmmpress ${busco_dataset}.hmm
    """

    stub:
    """
    touch ${busco_dataset}.hmm
    touch ${busco_dataset}.hmm.h3f
    touch ${busco_dataset}.hmm.h3i
    touch ${busco_dataset}.hmm.h3m
    touch ${busco_dataset}.hmm.h3p
    """
}