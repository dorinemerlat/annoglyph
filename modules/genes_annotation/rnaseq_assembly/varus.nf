process VARUS {
    tag "${id}"
    label 'varus'

    input:
    tuple val(id), val(meta), path(genome), path(runlist)

    output:
    tuple val(id), val(meta), path("${id}.bam"), emit: bam
    tuple val(id), val(meta), path("RunStatistics_${id}.csv"), emit: run_stats

    script:
    def scientificName = meta.lineage?.species?.name ?: meta.name
    def name = scientificName.replaceAll(/\s+/, '_').replace('.', '')
    def tokens = name.tokenize('_')
    def genus = tokens.first()
    def species = tokens.size() > 1 ? tokens[1..-1].join('_') : ''

    """
    mkdir -p ${name}
    cp ${runlist} ${name}/Runlist.txt

    # Write the VARUS configuration.
    cat > VARUSparameters.txt <<'EOF'
--batchSize 50000
--blockSize 5000
--components 1
--cost 0.001
--deleteLater 0
--estimator 2
--exportObservationsToFile 1
--exportParametersToFile 1
--fastqDumpCall fastq-dump
--genomeDir ./genome/
--lambda 10.0
--lessInfo 1
--loadAllOnce 0
--maxBatches 20
--mergeThreshold 10
--outFileNamePrefix ./
--pathToParameters ./VARUSparameters.txt
--pathToRuns ./
--pathToVARUS /opt/VARUS/Implementation/
--profitCondition 0
--pseudoCount 1
--qualityThreshold 5
--randomSeed 1
--readParametersFromFile 1
--runThreadN ${task.cpus}
--verbosityDebug 1
EOF

    # Reduce simultaneous requests to the SRA servers.
    sleep \$((RANDOM % 60 + 5))

    /opt/VARUS/runVARUS.pl \
        --aligner=STAR \
        --readFromTable=0 \
        --createindex=1 \
        --latinGenus=${genus} \
        --latinSpecies=${species} \
        --speciesGenome=${genome} \
        --logfile=varus.log \
        --nocreateRunList

    if [[ -f "${name}/VARUS.bam" ]]; then
        mv "${name}/VARUS.bam" "${id}.bam"
    else
        echo "WARNING: VARUS did not produce a BAM file for ${id}." >&2
        samtools view \
            -b \
            -T ${genome} \
            -o ${id}.bam \
            /dev/null
    fi

    if [[ -f "${name}/RunStatistics.csv" ]]; then
        mv "${name}/RunStatistics.csv" "RunStatistics_${id}.csv"
    else
        touch RunStatistics_${id}.csv
    fi
    """

    stub:
    """
    command -v samtools >/dev/null
    test -x /opt/VARUS/runVARUS.pl

    samtools view \
        -b \
        -T ${genome} \
        -o ${id}.bam \
        /dev/null

    touch RunStatistics_${id}.csv
    """
}