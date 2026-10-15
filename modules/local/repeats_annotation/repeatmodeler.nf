process REPEATMODELER {
    tag "${meta.id}"
    label 'repeatmodeler'

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${meta.id}-families.fa"), emit: families
    tuple val(meta), path("${meta.id}_BuildDatabase.log"), emit: build_database_log
    tuple val(meta), path("${meta.id}_RepeatModeler.log"), emit: repeatmodeler_log

    script:
    """
    set -o pipefail
    
    BuildDatabase \
        -name ${meta.id} \
        -engine ncbi \
        ${genome} \
        2>&1 | tee ${meta.id}_BuildDatabase.log

    RepeatModeler \
        -threads ${task.cpus} \
        -engine ncbi \
        -database ${meta.id} \
        -LTRStruct \
        2>&1 | tee -a ${meta.id}_RepeatModeler.log
    """

    stub:
    """
    command -v BuildDatabase >/dev/null
    command -v RepeatModeler >/dev/null

    cat > ${meta.id}-families.fa << EOF
>${meta.id}_stub_family#Unknown
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
EOF

    touch ${meta.id}_BuildDatabase.log ${meta.id}_RepeatModeler.log 
    """
}