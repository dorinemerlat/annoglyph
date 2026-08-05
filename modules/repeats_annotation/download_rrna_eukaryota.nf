process DOWNLOAD_RRNA_EUKARYOTA {
    tag "rrna"
    cpus 1
    label 'mchelper'
    time '2h'

    output:
    path "rRNA_Eukaryota.hmm*"

    script:
    """
    bsdtar -xf /opt/MCHelper/db/rRNA_Eukaryota.zip

    conda run -p /opt/conda/envs/MCHelper \
        hmmpress rRNA_Eukaryota.hmm
    """

    stub:
    """
    command -v bsdtar >/dev/null
    command -v conda >/dev/null

    touch rRNA_Eukaryota.hmm \
        rRNA_Eukaryota.hmm.h3f \
        rRNA_Eukaryota.hmm.h3i \
        rRNA_Eukaryota.hmm.h3m \
        rRNA_Eukaryota.hmm.h3p
    """
}