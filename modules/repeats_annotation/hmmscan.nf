process HMMSCAN {
    tag "${id}"
    label 'mchelper'

    input:
    tuple val(id), val(molecule), path(fasta), path(hmm_files)

    output:
    tuple val(id), val(molecule), path("${id}_tes_vs_${molecule}.tblout")

    script:
    def hmm_db_name = hmm_files.find { it.name.endsWith('.hmm') }

    """
    hmmscan --tblout ${id}_tes_vs_${molecule}.tblout -E 10 --noali --cpu ${task.cpus} $hmm_db_name $fasta
    """

    stub:
    """
    touch ${id}_tes_vs_${molecule}.tblout
    """
}