process INFERNAL {
    tag "${meta.id}"
    label 'infernal'
    cpus 30
    label 'program_infernal'
    input:
    tuple val(meta), path(genome), path(rfam)
    output:
    tuple val(meta), path("${meta.id}_infernal.tblout"), emit: tblout
    tuple val(meta), path("${meta.id}_infernal.out"), emit: report
    script:
    def databaseSize = meta.assembly_size * 2 / 1000000.0
    """
    cmscan -Z ${databaseSize} --cpu ${task.cpus} \
        --cut_ga --rfam --nohmmonly --fmt 2 --oclan --oskip \
        --clanin ${rfam}/Rfam.clanin \
        --tblout ${meta.id}_infernal.tblout -o ${meta.id}_infernal.out \
        ${params.program_args.infernal ?: ''} ${rfam}/Rfam.cm ${genome}
    """
    stub:
    """
    echo '# cmscan fmt 2' > ${meta.id}_infernal.tblout
    touch ${meta.id}_infernal.out
    """
}
