process BRAKER2 {
    tag "${id}"
    cpus 40
    time '10d'
    label 'braker'
    scratch false
    memory '170 GB'
    stageInMode 'copy'
    
    input:
    tuple val(id), val(meta), path(genome), path(proteins), val(specie_name)

    output:
    tuple val(id), val(meta), path("braker_${id}/braker.gff3"), emit: gff

    script:
    def name    = specie_name.replaceAll(/\s+/, '_').replace('.','')
    def toks    = name.tokenize('_')
    def genus   = toks.size() > 0 ? toks[0] : name
    def species = toks.size() > 1 ? toks[1..-1].join('_') : ''

    """
    cp -r /opt/Augustus/config/ .
    AUGUSTUS_CONFIG_PATH="\$PWD/config"

    braker.pl --genome=${genome} --prot_seq=${proteins} --gff3 --species=${genus}_${species} --threads ${task.cpus} \
        --AUGUSTUS_CONFIG_PATH=\$AUGUSTUS_CONFIG_PATH \
                                --AUGUSTUS_BIN_PATH="/opt/Augustus/bin/" \
                                --AUGUSTUS_SCRIPTS_PATH="/opt/Augustus/scripts/" \
                                --DIAMOND_PATH="/opt/ETP/tools/" \
                                --PROTHINT_PATH="/opt/ETP/bin/gmes/ProtHint/bin" \
                            --GENEMARK_PATH=/opt/ETP/bin/   --workingdir=./braker_${id}
    """

    stub:
    """
    touch braker.gff3
    """
}