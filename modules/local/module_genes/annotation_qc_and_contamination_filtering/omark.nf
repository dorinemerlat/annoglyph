process OMARK {
    tag "${annotation_method}/${id}"
    // label 'omark'
    stageInMode 'copy'
    scratch false
    memory '20 GB'
    
    input:
    tuple val(id), val(meta), val(annotation_method), path(proteins), path(luca_db)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_omark_detailed_summary.txt"), emit: summary
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_omark.sum"), emit: sum
    // tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_omark.omq"), emit: omq
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_omark.ump"), emit: ump
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_contaminants.count"), emit: contaminants_count

    script:
    """
    module load omark
    omamer search --db $luca_db --query $proteins --out ${id}_${annotation_method}_omark.omamer

    mkdir -p omark
    omark -f ${id}_${annotation_method}_omark.omamer -d $luca_db -o omark -t $meta.taxid

    # calculate number of contaminants detected by OMARK
    grep '^A:' omark/*sum | head -n 1 | sed -E 's/.*C:([0-9]+).*/\\1/' > ${id}_${annotation_method}_contaminants.count

    # move important results 
    mv omark/* .
    """

    stub:
    """
    touch ${id}_${annotation_method}_omark_detailed_summary.txt ${id}_${annotation_method}_omark.sum ${id}_${annotation_method}_omark.omq ${id}_${annotation_method}_omark.ump ${id}_${annotation_method}_contaminants.count
    """
}
