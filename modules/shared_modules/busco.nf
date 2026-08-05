process BUSCO {
    tag "${mode}/${prefix}"
    label 'busco'
    
    input:
    tuple val(id), val(meta), path(fasta), val(mode), val(prefix), val(busco_dataset), path(busco_dataset_path)

    output:
    tuple val(id), val(meta), val(prefix), val(mode), path("${prefix}_busco_${busco_dataset}.json"), emit: json
    tuple val(id), val(meta), val(prefix), val(mode), path("${prefix}_busco_${busco_dataset}.txt"), emit: txt
    tuple val(id), val(meta), val(prefix), val(mode), path("${prefix}_busco_${busco_dataset}.log"), emit: log

    script:
    def json_file = "${prefix}_busco_${busco_dataset}.json"
    def txt_file = "${prefix}_busco_${busco_dataset}.txt"
    def log_file = "${prefix}_busco_${busco_dataset}.log"
    
    """
    busco \
        --in ${fasta} \
        --out busco_${prefix} \
        --mode ${mode} \
        --lineage_dataset ./${busco_dataset_path} \
        --cpu ${task.cpus} \
        --offline \
        --force \
        2>&1 | tee ${log_file}

    mv busco_${prefix}/short_summary*.json ${json_file}
    mv busco_${prefix}/short_summary*.txt ${txt_file}
    """

    stub:
    """
    command -v busco >/dev/null

    echo '{}' > ${prefix}_busco_${busco_dataset}.json

    cat > ${prefix}_busco_${busco_dataset}.txt << EOF
BUSCO stub summary
Input: ${fasta}
Mode: ${mode}
Lineage: ${busco_dataset}
EOF

    touch ${prefix}_busco_${busco_dataset}.log
    """
}