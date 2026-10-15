process MCHELPER {
    tag "${meta.id}"
    label 'mchelper'

    input:
    tuple val(meta), path(genome), path(families), val(busco_dataset), path(busco_dataset_dir)

    output:
    tuple val(meta), path("classifiedModule/denovoLibTEs_PC.classif"), path("classifiedModule/kept_seqs_classified_module.fa"), path("unclassifiedModule/kept_seqs_unclassified_module.fa"), emit: main_results
    tuple val(meta), path("classifiedModule/kept_seqs_classified_module_curated.fa"), path("classifiedModule/kept_seqs_classified_module_non_curated.fa"), emit: classified_additional_results
    tuple val(meta), path("classifiedModule/fullLengthFrag.txt"), emit: full_length_frag
    tuple val(meta), path("classifiedModule/MSA_plots"), path("classifiedModule/MSA_seeds"), path("classifiedModule/te_aid"), emit: plots
    tuple val(meta), path("sequences_with_problems.txt"), emit: seqs_with_problems
    tuple val(meta), path("${meta.id}_mchelper.log"), emit: log

    script:
    """
    set +e

    MCHelper.py \
        -g ${genome} \
        -l ${families} \
        -o . \
        --input_type fasta \
        -r A \
        -t ${task.cpus} \
        -b ${busco_dataset_dir} \
        -a F \
        2>&1 | tee ${meta.id}_mchelper.log

    mchelper_status=\${PIPESTATUS[0]}

    set -e

    mkdir -p \
        classifiedModule/MSA_plots \
        classifiedModule/MSA_seeds \
        classifiedModule/te_aid \
        unclassifiedModule

    touch \
        classifiedModule/denovoLibTEs_PC.classif \
        classifiedModule/kept_seqs_classified_module.fa \
        classifiedModule/kept_seqs_classified_module_curated.fa \
        classifiedModule/kept_seqs_classified_module_non_curated.fa \
        classifiedModule/fullLengthFrag.txt \
        unclassifiedModule/kept_seqs_unclassified_module.fa \
        sequences_with_problems.txt

    if [[ \$mchelper_status -ne 0 ]] && \
       ! grep -q "Classified module didn't find any TE" ${meta.id}_mchelper.log; then
        exit \$mchelper_status
    fi
    """

    stub:
    """
    mkdir -p \
        classifiedModule/MSA_plots \
        classifiedModule/MSA_seeds \
        classifiedModule/te_aid \
        unclassifiedModule

    touch \
        classifiedModule/denovoLibTEs_PC.classif \
        classifiedModule/kept_seqs_classified_module.fa \
        classifiedModule/kept_seqs_classified_module_curated.fa \
        classifiedModule/kept_seqs_classified_module_non_curated.fa \
        classifiedModule/fullLengthFrag.txt \
        unclassifiedModule/kept_seqs_unclassified_module.fa \
        sequences_with_problems.txt \
        ${meta.id}_mchelper.log
    """
}