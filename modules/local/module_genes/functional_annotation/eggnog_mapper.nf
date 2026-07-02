process EGGNOG_MAPPER {

    tag "${id}/${annotation_method}"
    label 'exogap_eggnog'
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
          path("${id}_${annotation_method}.emapper.annotations"),
          emit: annotations

    script:
    """
    emapper.py \
        --cpu ${task.cpus} \
        --data_dir \$EGGNOG_DATA_DIR \
        -i ${fasta} \
        --itype proteins \
        --output ${id}_${annotation_method} \
        --output_dir .

    cp ${id}_${annotation_method}.emapper.annotations \
       ${id}_${annotation_method}.emapper.annotations

    test -s ${id}_${annotation_method}.emapper.annotations
    """

    stub:
    """
    touch ${id}_${annotation_method}.emapper.annotations
    """
}