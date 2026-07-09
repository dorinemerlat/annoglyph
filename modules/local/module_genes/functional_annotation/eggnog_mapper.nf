process EGGNOG_MAPPER {
    tag "${annotation_method}/${id}"
    cpus 8
    stageInMode 'copy'
    scratch false
    time '20h'
    memory { (16 * task.attempt).GB }

    input:
    tuple val(id), val(meta), val(annotation_method), path(proteins)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_eggnog/${id}_${annotation_method}_eggnog.emapper.annotations"), emit: annotations
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_eggnog/${id}_${annotation_method}_eggnog.emapper.seed_orthologs"), emit: seed_orthologs
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_eggnog/${id}_${annotation_method}_eggnog.emapper.hits"), emit: hits
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_eggnog/${id}_${annotation_method}_eggnog.emapper.orthologs"), emit: orthologs

    script:
    """
    module load eggnog-mapper

    mkdir -p ${id}_${annotation_method}_eggnog

    emapper.py -i ${proteins} \
        --cpu ${task.cpus} \
        --report_orthologs \
        --output ${id}_${annotation_method}_eggnog \
        --output_dir ${id}_${annotation_method}_eggnog \
        --scratch_dir /tmp \
        --temp_dir /tmp \
        --data_dir /shared/bank/emapperdb/5.0.2/
    """ 

    stub:
    """
    touch ${id}_${annotation_method}_eggnog/${id}_${annotation_method}_eggnog.emapper.annotations
    touch ${id}_${annotation_method}_eggnog/${id}_${annotation_method}_eggnog.emapper.seed_orthologs
    touch ${id}_${annotation_method}_eggnog/${id}_${annotation_method}_eggnog.emapper.hits
    touch ${id}_${annotation_method}_eggnog/${id}_${annotation_method}_eggnog.emapper.orthologs
    """
}