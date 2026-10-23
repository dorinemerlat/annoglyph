process BUSCO_PROTEINS {
    tag "${meta.annotation_method}/${meta.id}"
    label 'busco'

    label 'program_busco'
    input:
    tuple val(meta), path(proteins)
    path lineage

    output:
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_busco.json"), emit: json
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_busco.txt"), emit: txt

    script:
    def prefix = "${meta.id}_${meta.annotation_method}"
    """
    busco --in ${proteins} --out busco_${prefix} --mode proteins \
        --lineage_dataset "./${lineage}" --cpu ${task.cpus} --offline --force ${params.program_args.busco ?: ''}
    cp busco_${prefix}/short_summary*.json ${prefix}_busco.json
    cp busco_${prefix}/short_summary*.txt ${prefix}_busco.txt
    """

    stub:
    """
    echo '{}' > ${meta.id}_${meta.annotation_method}_busco.json
    touch ${meta.id}_${meta.annotation_method}_busco.txt
    """
}
