process HMMSCAN {
    tag "${id}"
    cpus 5
    label 'mchelper'
    time '2h'

    input:
    tuple val(id), val(molecule), path(fasta), path(hmm_files)

    output:
    // Save the raw hmmscan table and a filtered table with the grep selection
    tuple val(id), val(molecule), path("${id}_tes_vs_${molecule}.tblout")

    script:
    def hmm_db_name = hmm_files.find { it.name.endsWith('.hmm') }
    """
    conda run -p /opt/conda/envs/MCHelper hmmscan --tblout ${id}_tes_vs_${molecule}.tblout -E 10 --noali --cpu ${task.cpus} $hmm_db_name $fasta
    """
}
