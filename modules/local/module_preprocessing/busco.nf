process BUSCO_ON_ASSEMBLY {
    tag "on-assemblies/${id}"
    label 'busco'
    cpus 10
    scratch false

    input:
    tuple val(id), val(meta), path(proteins)

    output:
    tuple val(id), val(meta), path("${id}_on-assemblies_busco.json")
    tuple val(id), val(meta), path("${id}_on-assemblies_busco.txt")

    script:
    """
    busco -i $proteins -o busco_${id}_on-assemblies -m genome -l arthropoda_odb10 -c ${task.cpus} -f
    mv busco_${id}_on-assemblies/short_summary*.json ${id}_on-assemblies_busco.json
    mv busco_${id}_on-assemblies/short_summary*.txt ${id}_on-assemblies_busco.txt
    """
    
    stub:
    """
    touch ${id}_on-assemblies_busco.json
    touch ${id}_on-assemblies_busco.txt
    """
}