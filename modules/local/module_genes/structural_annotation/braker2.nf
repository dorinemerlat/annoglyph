process BRAKER2 {
    tag "${id}"
    cpus 40
    time '10d'
    label 'braker'
    scratch false
    memory '170 GB'
    stageInMode 'copy'
    errorStrategy 'retry'
    maxRetries 4
    
    input:
    tuple val(id), val(meta), path(genome), path(proteins)

    output:
    tuple val(id), val(meta), path("braker_${id}/braker.gff3"), emit: gff

    script:
    def name    = meta.lineage.species.name.replaceAll(/\s+/, '_').replace('.','')
    def toks    = name.tokenize('_')
    def genus   = toks.size() > 0 ? toks[0] : name
    def species = toks.size() > 1 ? toks[1..-1].join('_') : ''
    def braker_species = species ? "${genus}_${species}" : genus

    """
    set +e

    cp -r /opt/Augustus/config/ .
    export AUGUSTUS_CONFIG_PATH="\$PWD/config"

    braker.pl --genome=${genome} --prot_seq=${proteins} --gff3 --species=${braker_species} --threads ${task.cpus} \
        --AUGUSTUS_CONFIG_PATH="\$AUGUSTUS_CONFIG_PATH" \
        --AUGUSTUS_BIN_PATH="/opt/Augustus/bin/" \
        --AUGUSTUS_SCRIPTS_PATH="/opt/Augustus/scripts/" \
        --DIAMOND_PATH="/opt/ETP/tools/" \
        --PROTHINT_PATH="/opt/ETP/bin/gmes/ProtHint/bin" \
        --GENEMARK_PATH=/opt/ETP/bin/ --workingdir=./braker_${id}

    status=\$?

    if [ \$status -ne 0 ]; then
        echo "BRAKER failed for ${id} on attempt ${task.attempt}/5"

        if [ ${task.attempt} -lt 5 ]; then
            echo "Retrying..."
            exit \$status
        else
            echo "BRAKER failed after 5 attempts. Creating empty GFF3."

            mkdir -p braker_${id}
            touch braker_${id}/braker.gff3

            exit 0
        fi
    fi

    if [ ! -s braker_${id}/braker.gff3 ]; then
        echo "BRAKER finished but braker_${id}/braker.gff3 is missing or empty."

        if [ ${task.attempt} -lt 5 ]; then
            echo "Retrying..."
            exit 1
        else
            echo "Creating empty GFF3 after 5 attempts."

            mkdir -p braker_${id}
            touch braker_${id}/braker.gff3

            exit 0
        fi
    fi
    """

    stub:
    """
    mkdir -p braker_${id}
    touch braker_${id}/braker.gff3
    """
}
