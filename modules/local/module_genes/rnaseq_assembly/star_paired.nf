process STAR_PAIRED {
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
    tuple val(id), val(meta), path("${id}_paired.bam"), emit: bam
    tuple val(id), val(meta), path("${id}_paired_star_Log.out"),  emit: star_log
    tuple val(id), val(meta), path("${id}_paired_star_Log.progress.out"), emit: star_progress
    tuple val(id), val(meta), path("${id}_paired_star_Log.final.out"), emit: star_log_final
    
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

    r1_csv=\$(find input1 -maxdepth 1 -type f | sort | paste -sd, -)
    r2_csv=\$(find input2 -maxdepth 1 -type f | sort | paste -sd, -)
    r1_n=\$(find input1 -maxdepth 1 -type f  | wc -l)
    r2_n=\$(find input2 -maxdepth 1 -type f  | wc -l)

    if [ -z "\$r1_csv" ] || [ -z "\$r2_csv" ]; then
        echo "Missing paired-end reads for ${id}" >&2
        exit 1
    fi

    if [ "\$r1_n" -ne "\$r2_n" ]; then
        echo "Unbalanced paired-end read counts for ${id}: \$r1_n R1 vs \$r2_n R2" >&2
        exit 1
    fi

    STAR --runThreadN ${task.cpus} \
        --genomeDir genome/ \
        --readFilesIn "\$r1_csv" "\$r2_csv" \
        ${readCmd} \
        --outFileNamePrefix ${id}_paired_star_ \
        --outSAMstrandField intronMotif \
        --outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM ${task.memory.toBytes()}

    mv ${id}_paired_star_Aligned.sortedByCoord.out.bam ${id}_paired.bam

    # calculate the number of mapped reads 
    samtools view -c -F 4  ${id}_paired.bam > ${id}_paired_mapped_read_count.txt
    """

    stub:
    """
    touch ${id}_paired.bam ${id}_paired_star_Log.out ${id}_paired_star_Log.progress.out ${id}_paired_star_Log.final.out
    """
}