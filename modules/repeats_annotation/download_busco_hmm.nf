process DOWNLOAD_BUSCO_HMM {
    tag "${busco_dataset}"
    label 'mchelper'

    input:
    tuple val(busco_dataset), val(species_list)

    output:
    tuple val(busco_dataset), val(species_list), path("${busco_dataset}.hmm*")

    script:
    """
    # Download the matching BUSCO lineage
    hmm_list=\$(curl -s https://busco-data.ezlab.org/v4/data/lineages/ | grep -oP 'href="\\K[^"]+\\.tar\\.gz')
    matched_dataset=\$(grep ${busco_dataset} <<< "\$hmm_list")

    wget https://busco-data.ezlab.org/v4/data/lineages/\${matched_dataset}
    tar xvf \${matched_dataset}

    # Merge all HMM profiles into a single database
    cat ${busco_dataset}*/hmms/*.hmm > ${busco_dataset}.hmm

    # Index the HMM database
    conda run -p /opt/conda/envs/MCHelper/ hmmpress ${busco_dataset}.hmm
    """

    stub:
    """
    touch ${busco_dataset}.hmm \
        ${busco_dataset}.hmm.h3f \
        ${busco_dataset}.hmm.h3i \
        ${busco_dataset}.hmm.h3m \
        ${busco_dataset}.hmm.h3p
    """
}