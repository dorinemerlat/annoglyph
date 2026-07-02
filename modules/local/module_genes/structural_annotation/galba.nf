process GALBA {
    containerOptions "--home \$HOME"
    tag "${id}"
    cpus 25
    time '10d'
    memory '140 GB'
    label 'galba'
    maxRetries 5
    
    input:
    tuple val(id), val(meta), path(genome), path(proteins), val(specie_name)

    output:
    tuple val(id), val(meta), path("galba_${id}/galba.gff3"), emit: gff

    script:
    def name    = specie_name.replaceAll(/\s+/, '_').replace('.','')
    def toks    = name.tokenize('_')
    def genus   = toks.size() > 0 ? toks[0] : name
    def species = toks.size() > 1 ? toks[1..-1].join('_') : ''

    """
    cp -r /opt/Augustus/config/ .
    AUGUSTUS_CONFIG_PATH="\$PWD/config"

    galba.pl --species=${genus}_${species} \
        --genome=${genome} \
        --prot_seq=${proteins} \
        --gff3 \
        --threads ${task.cpus} \
        --workingdir galba_${id} \
        --AUGUSTUS_CONFIG_PATH=\$AUGUSTUS_CONFIG_PATH
    """

    stub:
    """
    touch galba_${id}.gff3
    """
}