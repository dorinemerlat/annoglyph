process DOWNLOAD_BUSCO_HMM{
    tag "${busco_dataset}"
    cpus 1
    label 'mchelper'
    time '2h'

    input: 
    tuple val(busco_dataset), val(species_list)
    
    output:
    tuple val(busco_dataset), val(species_list), path("${busco_dataset}.hmm*")

    script:
    """
    hmm_list=\$(curl -s https://busco-data.ezlab.org/v4/data/lineages/ | grep -oP 'href="\\K[^"]+\\.tar\\.gz' )
    matched_dataset=\$(grep $busco_dataset <<< "\$hmm_list")

    wget https://busco-data.ezlab.org/v4/data/lineages/\${matched_dataset}
    tar xvf \${matched_dataset}
    cat ${busco_dataset}*/hmms/*.hmm > ${busco_dataset}.hmm

    conda run -p /opt/conda/envs/MCHelper/ hmmpress ${busco_dataset}.hmm
    """
}
