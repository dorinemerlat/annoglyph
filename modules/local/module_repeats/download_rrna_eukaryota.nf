process DOWNLOAD_RRNA_EUKARYOTA {
    tag "rrna"
    cpus 1
    label 'mchelper'
    time '2h'

    output:
    path("rRNA_Eukaryota.hmm*")

    script:
    """
    # unzip the HMM database if it's zipped, and prepare it for hmmscan
    bsdtar -xf  /opt/MCHelper/db/rRNA_Eukaryota.zip

    # Prepare the HMM database for hmmscan
    conda run -p /opt/conda/envs/MCHelper/ hmmpress rRNA_Eukaryota.hmm
    """
}
