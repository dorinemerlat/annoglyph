process BRAKER3 {
    tag "${meta.rnaseq_aligner}/${meta.id}"
    label 'braker'
    scratch false
    
    label 'program_braker3'
    input:
    tuple val(meta), path(genome), path(bam), path(proteins)

    output:
    tuple val(meta), path("braker_${meta.id}/braker.gff3"), emit: gff

    script:
    def geneticCode = meta.translation_table ?: params.translation_table
    def geneticCodeArg = geneticCode.toString() == '1' ? '' : "--translation_table=${geneticCode}"
    def scientificName = meta.lineage?.species?.name ?: meta.name
    def name = scientificName.replaceAll(/\s+/, '_').replace('.', '')
    def tokens = name.tokenize('_')
    def genus = tokens.first()
    def species = tokens.size() > 1 ? tokens[1..-1].join('_') : ''
    def brakerSpecies = species ? "${genus}_${species}" : genus
    def fungusOption = meta.lineage?.values()?.any { taxon -> taxon?.name == 'Fungi' } ? '--fungus' : ''
    
    """
    cp -r /opt/Augustus/config .
    export AUGUSTUS_CONFIG_PATH="\$PWD/config"
    cp -r /opt/ETP/bin genemark_etp
    python3 ${projectDir}/bin/patch_genemark_super.py genemark_etp/super.pl

    set +e

    braker.pl \
        --genome=${genome} \
        --prot_seq=${proteins} \
        --bam=${bam} \
        --gff3 \
        --species=${brakerSpecies} \
        --threads=${task.cpus} \
        ${fungusOption} \
        --AUGUSTUS_CONFIG_PATH="\${AUGUSTUS_CONFIG_PATH}" \
        --AUGUSTUS_BIN_PATH=/opt/Augustus/bin \
        --AUGUSTUS_SCRIPTS_PATH=/opt/Augustus/scripts \
        --DIAMOND_PATH=/opt/ETP/tools \
        --PROTHINT_PATH=/opt/ETP/bin/gmes/ProtHint/bin \
        --GENEMARK_PATH="\$PWD/genemark_etp" \
        --workingdir=braker_${meta.id} ${params.program_args.braker3 ?: ''} ${geneticCodeArg}

    braker_status=\$?

    set -e

    if [[ \$braker_status -ne 0 || ! -s braker_${meta.id}/braker.gff3 ]]; then
        echo "BRAKER3 failed for ${meta.id} using ${meta.rnaseq_aligner} on attempt ${task.attempt}/5." >&2

        if [[ ${task.attempt} -lt 5 ]]; then
            exit \$([[ \$braker_status -ne 0 ]] && echo \$braker_status || echo 1)
        fi

        echo "BRAKER3 failed after five attempts. Creating an empty GFF3 file." >&2
        mkdir -p braker_${meta.id}
        : > braker_${meta.id}/braker.gff3
    fi
    """

    stub:
    """
    command -v braker.pl >/dev/null
    test -d /opt/Augustus/config
    test -x /opt/Augustus/bin/augustus

    mkdir -p braker_${meta.id}
    echo '##gff-version 3' > braker_${meta.id}/braker.gff3
    """
}
