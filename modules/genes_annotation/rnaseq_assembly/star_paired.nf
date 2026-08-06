process STAR_PAIRED {
    tag "${id}"
    label 'varus'

    input:
    tuple val(id), val(meta), path(genome), path(reads, stageAs: "input*/*")

    output:
    tuple val(id), val(meta), path("${id}_paired.bam"), emit: bam
    tuple val(id), val(meta), path("${id}_paired_star_Log.out"), emit: star_log
    tuple val(id), val(meta), path("${id}_paired_star_Log.progress.out"), emit: star_progress
    tuple val(id), val(meta), path("${id}_paired_star_Log.final.out"), emit: star_log_final

    script:
    def readFiles = reads instanceof List ? reads.flatten() : [reads]
    def compressedReads = readFiles.count { it.name.endsWith('.gz') }

    if (compressedReads > 0 && compressedReads < readFiles.size()) {
        error "Mixed compressed and uncompressed paired-end reads for '${id}'."
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

    r1_csv=\$(find input1 -maxdepth 1 -type f | sort | paste -sd, -)
    r2_csv=\$(find input2 -maxdepth 1 -type f | sort | paste -sd, -)
    r1_count=\$(find input1 -maxdepth 1 -type f | wc -l)
    r2_count=\$(find input2 -maxdepth 1 -type f | wc -l)

    if [[ -z "\${r1_csv}" || -z "\${r2_csv}" ]]; then
        echo "Missing paired-end reads for ${id}." >&2
        exit 1
    fi

    if [[ "\${r1_count}" -ne "\${r2_count}" ]]; then
        echo "Unbalanced paired-end reads for ${id}: \${r1_count} R1 and \${r2_count} R2 files." >&2
        exit 1
    fi

    # Align paired-end reads and sort the BAM by coordinates.
    STAR \
        --runThreadN ${task.cpus} \
        --genomeDir genome \
        --readFilesIn "\${r1_csv}" "\${r2_csv}" \
        ${readCommand} \
        --outFileNamePrefix ${id}_paired_star_ \
        --outSAMstrandField intronMotif \
        --outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM ${task.memory.toBytes()}

    mv ${id}_paired_star_Aligned.sortedByCoord.out.bam ${id}_paired.bam
    """

    stub:
    """
    command -v STAR >/dev/null

    touch \
        ${id}_paired.bam \
        ${id}_paired_star_Log.out \
        ${id}_paired_star_Log.progress.out \
        ${id}_paired_star_Log.final.out
    """
}
