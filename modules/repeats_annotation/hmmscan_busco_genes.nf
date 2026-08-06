process HMMSCAN_BUSCO_GENES {
    tag "${id}"
    label 'mchelper'

    input:
    tuple val(id), val(meta), path(fasta), val(busco_dataset), path(busco_dataset_dir)

    output:
    tuple val(id), val(meta), path("${id}_tes_vs_genes.tblout")

    script:
    """
    cat ${busco_dataset_dir}/hmms/*.hmm > ${busco_dataset}.hmm

    hmmpress ${busco_dataset}.hmm

    hmmscan \
        --tblout ${id}_tes_vs_genes.tblout \
        -E 10 \
        --noali \
        --cpu ${task.cpus} \
        ${busco_dataset}.hmm \
        ${fasta}
    """

    stub:
    """
    command -v hmmpress >/dev/null
    command -v hmmscan >/dev/null

    touch ${id}_tes_vs_genes.tblout
    """
}