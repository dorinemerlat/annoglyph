process INFERNAL {

    tag "${id}"
    label 'infernal'
    cpus 30

    input:
    tuple val(id), val(meta), path(genome), path{rfam_clanin}, path{rfam_cm}

    output:
    tuple val(id), val(meta), path("${id}_infernal.gff"), emit: gff
    tuple val(id), val(meta), path("${id}_infernal.fa"), emit: fasta

    script:
    database_size = ${meta.assembly_size} * 2 / 1000000
    """
    cmscan -Z $database_size --cpu ${task.cpus} --cut_ga --rfam --nohmmonly --tblout ${id}_infernal.tblout -o ${id}_infernal.out \
            --fmt 2 --clanin ${rfam_clanin} ${rfam_cm} $genome
    """
    
    stub:
    """
    touch file
    """
}