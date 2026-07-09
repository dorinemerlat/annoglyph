process GALBA {
    containerOptions "--home \$HOME"
    tag "${id}"
    cpus 25
    time '10d'
    memory '140 GB'
    label 'galba'
    errorStrategy 'retry'
    maxRetries 4
    
    input:
    tuple val(id), val(meta), path(genome), path(proteins)

    output:
    tuple val(id), val(meta), path("galba_${id}/galba.gff3"), emit: gff

    script:
    def name    = meta.lineage.species.name.replaceAll(/\s+/, '_').replace('.','')
    def toks    = name.tokenize('_')
    def genus   = toks.size() > 0 ? toks[0] : name
    def species = toks.size() > 1 ? toks[1..-1].join('_') : ''

    """
    cp -r /opt/Augustus/config/ .
    AUGUSTUS_CONFIG_PATH="\$PWD/config"
    set +e

    cp -r /opt/Augustus/config/ .
    export AUGUSTUS_CONFIG_PATH="\$PWD/config"

    galba.pl --species=${genus}_${species} \
        --genome=${genome} \
        --prot_seq=${proteins} \
        --gff3 \
        --threads ${task.cpus} \
        --workingdir galba_${id} \
        --AUGUSTUS_CONFIG_PATH=\$AUGUSTUS_CONFIG_PATH

    status=\$?

    if [ \$status -ne 0 ]; then
        echo "GALBA failed for ${id} on attempt ${task.attempt}/5"

        if [ ${task.attempt} -lt 5 ]; then
            echo "Retrying..."
            exit \$status
        else
            echo "GALBA failed after 5 attempts. Creating empty GFF3."

            mkdir -p galba_${id}
            touch galba_${id}/galba.gff3

            exit 0
        fi
    fi

    if [ ! -s galba_${id}/galba.gff3 ]; then
        echo "GALBA finished but galba_${id}/galba.gff3 is missing or empty."

        if [ ${task.attempt} -lt 5 ]; then
            echo "Retrying..."
            exit 1
        else
            echo "Creating empty GFF3 after 5 attempts."

            mkdir -p galba_${id}
            touch galba_${id}/galba.gff3

            exit 0
        fi
    fi
    """

    stub:
    """
    mkdir -p galba_${id}
    touch galba_${id}/galba.gff3
    """
}
