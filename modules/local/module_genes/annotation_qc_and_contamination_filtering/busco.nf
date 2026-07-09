process BUSCO {
    tag "${annotation_method}/${id}"
    label 'busco'
    cpus 10
    scratch false

    input:
    tuple val(id), val(meta), val(annotation_method), path(proteins)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_busco.json"), emit: json
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_busco.txt"), emit: txt

    script:
    """
    busco -i $proteins -o busco_${id}_${annotation_method} -m prot -l arthropoda_odb10 -c ${task.cpus} -f
    mv busco_${id}_${annotation_method}/short_summary*.json ${id}_${annotation_method}_busco.json
    mv busco_${id}_${annotation_method}/short_summary*.txt ${id}_${annotation_method}_busco.txt
    """
    
    stub:
    """
    touch ${id}_${annotation_method}_busco.json
    touch ${id}_${annotation_method}_busco.txt
    """
}