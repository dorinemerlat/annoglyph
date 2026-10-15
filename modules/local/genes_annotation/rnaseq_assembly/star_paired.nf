process STAR_PAIRED {
    tag "${meta.id}"
    label 'varus'

    label 'program_star'
    input:
    tuple val(meta), path(genome), path(reads, stageAs: "input*/*")

    output:
    tuple val(meta), path("${meta.id}_paired.bam"), emit: bam
    tuple val(meta), path("${meta.id}_paired_star_Log.out"), emit: star_log
    tuple val(meta), path("${meta.id}_paired_star_Log.progress.out"), emit: star_progress
    tuple val(meta), path("${meta.id}_paired_star_Log.final.out"), emit: star_log_final

    script:
    def readFiles = reads instanceof List ? reads.flatten() : [reads]
    def compressedReads = readFiles.count { it.name.endsWith('.gz') }

    if (compressedReads > 0 && compressedReads < readFiles.size()) {
        error "Mixed compressed and uncompressed paired-end reads for '${meta.id}'."
    }

    if (readFiles.size() % 2 != 0) {
        error "Unbalanced paired-end reads for '${meta.id}'."
    }
    def pairs = readFiles.collate(2)
    def r1Csv = pairs.collect { it[0] }.join(',')
    def r2Csv = pairs.collect { it[1] }.join(',')

    def readCommand = compressedReads == readFiles.size() ? '--readFilesCommand zcat' : ''

    """
    mkdir genome

    # Generate the STAR genome index.
    STAR \
        --runThreadN ${task.cpus} \
        --runMode genomeGenerate \
        --outTmpDir STARtmp \
        --genomeDir genome \
        --genomeFastaFiles ${genome} \
        --limitGenomeGenerateRAM ${task.memory.toBytes()} ${params.program_args.star_index ?: ''}

    # Align paired-end reads and sort the BAM by coordinates.
    STAR \
        --runThreadN ${task.cpus} \
        --genomeDir genome \
        --readFilesIn "${r1Csv}" "${r2Csv}" \
        ${readCommand} \
        --outFileNamePrefix ${meta.id}_paired_star_ \
        --outSAMstrandField intronMotif \
        --outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM ${task.memory.toBytes()} ${params.program_args.star ?: ''}

    mv ${meta.id}_paired_star_Aligned.sortedByCoord.out.bam ${meta.id}_paired.bam
    """

    stub:
    """
    command -v STAR >/dev/null

    touch \
        ${meta.id}_paired.bam \
        ${meta.id}_paired_star_Log.out \
        ${meta.id}_paired_star_Log.progress.out \
        ${meta.id}_paired_star_Log.final.out
    """
}
