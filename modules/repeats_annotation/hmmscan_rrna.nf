process HMMSCAN_RRNA {
    tag "${id}"
    label 'mchelper'

    input:
    tuple val(id), val(meta), path(fasta)

    output:
    tuple val(id), val(meta), path("${id}_tes_vs_rRNA.tblout")

    script:
    """
    cp /opt/MCHelper/db/rRNA_Eukaryota.hmm .

    hmmpress rRNA_Eukaryota.hmm

    hmmscan \
        --tblout ${id}_tes_vs_rRNA.tblout \
        -E 10 \
        --noali \
        --cpu ${task.cpus} \
        rRNA_Eukaryota.hmm \
        ${fasta}
    """

    stub:
    """
    command -v hmmpress >/dev/null
    command -v hmmscan >/dev/null

    touch ${id}_tes_vs_rRNA.tblout
    """
}