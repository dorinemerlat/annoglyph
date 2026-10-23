process OMARK {
    tag "${meta.annotation_method}/${meta.id}"
    label 'omark'

    label 'program_omark'
    input:
    tuple val(meta), path(proteins), path(luca_db)
    path ete_taxonomy

    output:
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_omark_detailed_summary.txt"), emit: summary
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_omark.sum"), emit: sum
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_omark.ump"), emit: ump
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_contaminants.count"), emit: contaminants_count

    script:
    def prefix = "${meta.id}_${meta.annotation_method}_omark"

    def taxonomyArg = ete_taxonomy ? "--ete_ncbi_db '${ete_taxonomy}'" : '--ete_ncbi_db "$PWD/omark_taxonomy.sqlite"'

    """
    omamer search \
        --db ${luca_db} \
        --query ${proteins} \
        --out ${prefix}.omamer

    mkdir omark

    omark \
        --file ${prefix}.omamer \
        --database ${luca_db} \
        --outputFolder omark \
        --taxid ${meta.taxid} ${taxonomyArg} ${params.program_args.omark ?: ''}

    mv omark/* .

    awk '
        /^>/ { contaminated = (\$0 ~ /^>Contamination(_|\$)/); next }
        contaminated && NF { count++ }
        END { print count + 0 }
    ' ${prefix}.ump > ${meta.id}_${meta.annotation_method}_contaminants.count
    """

    stub:
    """
    touch \
        ${meta.id}_${meta.annotation_method}_omark_detailed_summary.txt \
        ${meta.id}_${meta.annotation_method}_omark.sum \
        ${meta.id}_${meta.annotation_method}_omark.ump \
        ${meta.id}_${meta.annotation_method}_contaminants.count
    echo 0 > ${meta.id}_${meta.annotation_method}_contaminants.count
    """
}
