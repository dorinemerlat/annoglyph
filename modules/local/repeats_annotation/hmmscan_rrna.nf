process HMMSCAN_RRNA {
    tag "${meta.id}"
    label 'mchelper'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_tes_vs_rRNA.tblout")

    script:
    """
    cp /opt/MCHelper/db/rRNA_Eukaryota.hmm .

    hmmpress rRNA_Eukaryota.hmm

    hmmscan \
        --tblout ${meta.id}_tes_vs_rRNA.tblout \
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

    touch ${meta.id}_tes_vs_rRNA.tblout
    """
}