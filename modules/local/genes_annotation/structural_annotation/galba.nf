process GALBA {
    containerOptions "--home \$PWD"
    tag "${meta.id}"
    label 'galba'
    scratch false

    label 'program_galba'
    input:
    tuple val(meta), path(genome), path(proteins)

    output:
    tuple val(meta), path("galba_${meta.id}/galba.gff3"), emit: gff

    script:
    def geneticCode = meta.translation_table ?: params.translation_table
    def geneticCodeArg = geneticCode.toString() == '1' ? '' : "--translation_table=${geneticCode}"
    def name    = (meta.lineage?.species?.name ?: meta.name).replaceAll(/\s+/, '_').replace('.','')
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
        --workingdir galba_${meta.id} \
        --AUGUSTUS_CONFIG_PATH=\$AUGUSTUS_CONFIG_PATH ${params.program_args.galba ?: ''} ${geneticCodeArg}

    status=\$?

    if [ \$status -ne 0 ]; then
        echo "GALBA failed for ${meta.id} on attempt ${task.attempt}/5"

        if [ ${task.attempt} -lt 5 ]; then
            echo "Retrying..."
            exit \$status
        else
            echo "GALBA failed after 5 attempts. Creating empty GFF3."

            mkdir -p galba_${meta.id}
            touch galba_${meta.id}/galba.gff3

            exit 0
        fi
    fi

    if [ ! -s galba_${meta.id}/galba.gff3 ]; then
        echo "GALBA finished but galba_${meta.id}/galba.gff3 is missing or empty."

        if [ ${task.attempt} -lt 5 ]; then
            echo "Retrying..."
            exit 1
        else
            echo "Creating empty GFF3 after 5 attempts."

            mkdir -p galba_${meta.id}
            touch galba_${meta.id}/galba.gff3

            exit 0
        fi
    fi
    """

    stub:
    """
    mkdir -p galba_${meta.id}
    echo '##gff-version 3' > galba_${meta.id}/galba.gff3
    """
}
