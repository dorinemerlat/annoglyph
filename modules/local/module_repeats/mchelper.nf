process MCHELPER {
    tag "${id}"
    cpus 50
    time '20d'
    label 'mchelper'
    errorStrategy 'ignore'

    input:
    tuple val(id), val(meta), path(families), path(genome)

    output:
    tuple val(id), val(meta), path("classifiedModule/denovoLibTEs_PC.classif"), path("classifiedModule/kept_seqs_classified_module.fa"), path("unclassifiedModule/kept_seqs_unclassified_module.fa"), emit: main_results
    tuple val(id), val(meta), path("classifiedModule/kept_seqs_classified_module_curated.fa"), path("classifiedModule/kept_seqs_classified_module_non_curated.fa"), emit: classfied_additional_results
    tuple val(id), val(meta), path("classifiedModule/fullLengthFrag.txt"), emit: full_length_frag
    tuple val(id), val(meta), path("classifiedModule/MSA_plots"), path("classifiedModule/MSA_seeds"), path("classifiedModule/te_aid"), emit: plots
    tuple val(id), val(meta), path("sequences_with_problems.txt"), emit: seqs_with_problems
    script:
    """
    conda run -p /opt/conda/envs/MCHelper \
        python3 MCHelper/MCHelper.py \
        -g ${genome} \
        -l ${families} \
        -o . \
        --input_type fasta \
        -r A \
        -t ${task.cpus} \
        -b arthropoda_odb10.hmm \
        -a F
    """

    stub:
    """
    touch
    """
}
