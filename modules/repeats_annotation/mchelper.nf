process MCHELPER {
    tag "${id}"
    label 'mchelper'
    scratch false
    
    input:
    tuple val(id), val(meta), path(genome), path(families), val(busco_dataset), path(busco_dataset_dir)

    output:
    tuple val(id), val(meta), path("classifiedModule/denovoLibTEs_PC.classif"), path("classifiedModule/kept_seqs_classified_module.fa"), path("unclassifiedModule/kept_seqs_unclassified_module.fa"), emit: main_results
    tuple val(id), val(meta), path("classifiedModule/kept_seqs_classified_module_curated.fa"), path("classifiedModule/kept_seqs_classified_module_non_curated.fa"), emit: classified_additional_results
    tuple val(id), val(meta), path("classifiedModule/fullLengthFrag.txt"), emit: full_length_frag
    tuple val(id), val(meta), path("classifiedModule/MSA_plots"), path("classifiedModule/MSA_seeds"), path("classifiedModule/te_aid"), emit: plots
    tuple val(id), val(meta), path("sequences_with_problems.txt"), emit: seqs_with_problems
    tuple val(id), val(meta), path("${id}_mchelper.log"), emit: log

    script:
    """
    set -o pipefail

    MCHelper.py \
        -g ${genome} \
        -l ${families} \
        -o . \
        --input_type fasta \
        -r A \
        -t ${task.cpus} \
        -b ${busco_dataset_dir} \
        -a F \
        2>&1 | tee ${id}_mchelper.log
    """

    stub:
    """
    mkdir -p classifiedModule/MSA_plots
    mkdir -p classifiedModule/MSA_seeds
    mkdir -p classifiedModule/te_aid
    mkdir -p unclassifiedModule

    touch classifiedModule/denovoLibTEs_PC.classif
    touch classifiedModule/kept_seqs_classified_module.fa
    touch unclassifiedModule/kept_seqs_unclassified_module.fa
    touch classifiedModule/kept_seqs_classified_module_curated.fa
    touch classifiedModule/kept_seqs_classified_module_non_curated.fa
    touch classifiedModule/fullLengthFrag.txt
    touch sequences_with_problems.txt
    touch ${id}_mchelper.log
    """
}