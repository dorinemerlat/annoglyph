process STAR_SINGLE {
    tag "${id}"
    label 'varus'

    input:
    tuple val(id), val(meta), path(genome), path(reads, stageAs: "input*/*")

    output:
    tuple val(id), val(meta), path("${id}_single.bam"), emit: bam
    tuple val(id), val(meta), path("${id}_single_star_Log.out"), emit: star_log
    tuple val(id), val(meta), path("${id}_single_star_Log.progress.out"), emit: star_progress
    tuple val(id), val(meta), path("${id}_single_star_Log.final.out"), emit: star_log_final

    script:
    def readFiles = reads instanceof List ? reads.flatten() : [reads]
    def compressedReads = readFiles.count { it.name.endsWith('.gz') }

    if (compressedReads > 0 && compressedReads < readFiles.size()) {
        error "Mixed compressed and uncompressed single-end reads for '${id}'."
    }

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
        --limitGenomeGenerateRAM ${task.memory.toBytes()}

    single_reads_csv=\$(find input* -maxdepth 1 -type f | sort | paste -sd, -)

    if [[ -z "\${single_reads_csv}" ]]; then
        echo "No single-end reads found for ${id}." >&2
        exit 1
    fi

    # Align single-end reads and sort the BAM by coordinates.
    STAR \
        --runThreadN ${task.cpus} \
        --genomeDir genome \
        --readFilesIn "\${single_reads_csv}" \
        ${readCommand} \
        --outFileNamePrefix ${id}_single_star_ \
        --outSAMstrandField intronMotif \
        --outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM ${task.memory.toBytes()}

    mv ${id}_single_star_Aligned.sortedByCoord.out.bam ${id}_single.bam
    """

    stub:
    """
    command -v STAR >/dev/null

    touch \
        ${id}_single.bam \
        ${id}_single_star_Log.out \
        ${id}_single_star_Log.progress.out \
        ${id}_single_star_Log.final.out
    """
}
