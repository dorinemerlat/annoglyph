process REPEATMODELER {
    tag "${id}"
    label 'repeatmodeler'

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path("${id}-families.fa"), emit: families
    tuple val(id), val(meta), path("${id}_BuildDatabase.log"), emit: build_database_log
    tuple val(id), val(meta), path("${id}_RepeatModeler.log"), emit: repeatmodeler_log

    script:
    """
    set -o pipefail
    
    BuildDatabase \
        -name ${id} \
        -engine ncbi \
        ${genome} \
        2>&1 | tee ${id}_BuildDatabase.log

    RepeatModeler \
        -threads ${task.cpus} \
        -engine ncbi \
        -database ${id} \
        -LTRStruct \
        2>&1 | tee -a ${id}_RepeatModeler.log
    """

    stub:
    """
    command -v BuildDatabase >/dev/null
    command -v RepeatModeler >/dev/null

    cat > ${id}-families.fa << EOF
>${id}_stub_family#Unknown
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
EOF

    touch ${id}_BuildDatabase.log ${id}_RepeatModeler.log 
    """
}