process VARUS {
    tag "${id}"
    cpus 30
    time '5d'
    label 'varus'
    scratch false
    stageInMode 'copy'
    memory '64 GB'

    input:
    tuple val(id), val(meta), path(genome), path(runlist)

    output:
    tuple val(id), val(meta), path("${id}.bam"), emit: bam
    tuple val(id), val(meta), path("RunStatistics_${id}.csv"), emit: run_stats

    script:
    def name    = meta.lineage.species.name.replaceAll(/\s+/, '_').replace('.','')
    def toks    = name.tokenize('_')
    def genus   = toks.size() > 0 ? toks[0] : name
    def species = toks.size() > 1 ? toks[1..-1].join('_') : ''

    """
    mkdir ${name}
    cp ${runlist} ${name}/Runlist.txt

    # generate VARUS parameters file
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

    # A random sleep to avoid multiple VARUS instances starting at the same time and hitting NCBI at the same time
    sleep \$(( (RANDOM % 60) + 5 ))

    # Run Varus
    /opt/VARUS/runVARUS.pl --aligner=STAR --readFromTable=0 --createindex=1 \\
        --latinGenus=${genus} --latinSpecies=${species} \\
        --speciesGenome=${genome} --logfile=varus.log --nocreateRunList \\
        > varus.log 

    if [[ -f "${name}/VARUS.bam" ]]; then
        mv "${name}/VARUS.bam" "${id}.bam"
        mv "${name}/RunStatistics.csv" "RunStatistics_${id}.csv"
    else
        touch "${id}.bam"
    fi
    """

    stub:
    """
    touch ${id}.bam RunStatistics_${id}.csv
    """
}