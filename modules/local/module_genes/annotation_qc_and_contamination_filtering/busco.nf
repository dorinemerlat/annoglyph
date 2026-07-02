process BUSCO {

    tag "${id}/${annotation_method}"
    label 'exogap_busco'
    scratch false

    input:
    tuple val(id),
          val(meta),
          val(annotation_method),
          path(fasta)

    output:
    tuple val(id),
          val(meta),
          val(annotation_method),
          path("${id}_${annotation_method}_busco.json"),
          emit: json

    script:
    """
    busco \\
        --in ${fasta} \\
        --out busco_${id}_${annotation_method} \\
        --mode proteins \\
        --lineage_dataset ${meta.busco_dataset} \\
        --cpu ${task.cpus} \\
        --offline \\
        --force

    cp busco_${id}_${annotation_method}/short_summary.specific.*.json \\
       ${id}_${annotation_method}_busco.json

    test -s ${id}_${annotation_method}_busco.json
    """

    stub:
    """
    cat <<EOF > ${id}_${annotation_method}_busco.json
{}
EOF
    """
}