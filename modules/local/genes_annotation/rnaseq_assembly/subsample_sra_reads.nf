process SUBSAMPLE_SRA_READS {
    tag "${meta.id}/${reads instanceof List ? reads.first().name : reads.name}/${params.sra_max_spots}"
    label 'annoglyph_python'
    input:
    tuple val(meta), path(reads, stageAs: 'input/*')
    output:
    tuple val(meta), path('sampled/*.fastq.gz'), emit: reads
    path 'subsampling.log', emit: log
    script:
    def files = reads instanceof List ? reads : [reads]
    def inputArgs = files.collect { "'${it}'" }.join(' ')
    """
    subsample_sra_reads.py --reads ${inputArgs} \
        --max-spots ${params.sra_max_spots} --output sampled > subsampling.log
    """
}
