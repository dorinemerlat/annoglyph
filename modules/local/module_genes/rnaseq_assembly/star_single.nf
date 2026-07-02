process STAR_SINGLE {
    tag "${id}"
    cpus 20
    time '5d'
    label 'varus'
    scratch false
    memory { (100 * task.attempt).GB }
    stageInMode 'copy'
    maxRetries 10

    input:
    tuple val(id), val(meta), path(genome), path(reads, stageAs: "input*/*")

    output:
    tuple val(id), val(meta), path("${id}_single.bam"),  emit: bam
    tuple val(id), val(meta), path("${id}_single_star_Log.out"),  emit: star_log
    tuple val(id), val(meta), path("${id}_single_star_Log.progress.out"), emit: star_progress
    tuple val(id), val(meta), path("${id}_single_star_Log.final.out"), emit: star_log_final

    script:
    def readFiles = (reads instanceof List) ? reads : [reads]
    def readCmd = readFiles.any { it.name.endsWith('.gz') } ? "--readFilesCommand zcat" : ""
    """
    set -euo pipefail

    mkdir genome

    STAR \
    --runThreadN 4 \
    --runMode genomeGenerate \
    --outTmpDir STARtmp \
    --genomeDir genome \
    --genomeFastaFiles ${genome} \
    --limitGenomeGenerateRAM ${task.memory.toBytes()}

    single_reads_csv=\$(find input* -maxdepth 1 -type f | sort | paste -sd, -)

    if [ -z "\$single_reads_csv" ]; then
        echo "No single-end reads found for ${id}" >&2
        exit 1
    fi

    STAR --runThreadN ${task.cpus} \
        --genomeDir genome/ \
        --readFilesIn "\$single_reads_csv" \
        ${readCmd} \
        --outFileNamePrefix ${id}_single_star_ \
        --outSAMstrandField intronMotif \
        --outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM ${task.memory.toBytes()}

    mv ${id}_single_star_Aligned.sortedByCoord.out.bam ${id}_single.bam

    # calculate the number of mapped reads 
    samtools view -c -F 4  ${id}_single.bam > ${id}_single_mapped_read_count.txt
    """

    stub:
    """
    touch ${id}_single.bam ${id}_single_star_Log.out ${id}_single_star_Log.progress.out ${id}_single_star_Log.final.out
    """
}