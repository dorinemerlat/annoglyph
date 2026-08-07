process BRAKER2 {
    tag "${id}"
    label 'braker'
    scratch false
    input:
    tuple val(id), val(meta), path(genome), path(proteins)

    output:
    tuple val(id), val(meta), path("braker_${id}/braker.gff3"), emit: gff

    script:
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

    set +e

    braker.pl \
        --genome=${genome} \
        --prot_seq=${proteins} \
        --gff3 \
        --species=${brakerSpecies} \
        --threads=${task.cpus} \
        ${fungusOption} \
        --AUGUSTUS_CONFIG_PATH="\${AUGUSTUS_CONFIG_PATH}" \
        --AUGUSTUS_BIN_PATH=/opt/Augustus/bin \
        --AUGUSTUS_SCRIPTS_PATH=/opt/Augustus/scripts \
        --DIAMOND_PATH=/opt/ETP/tools \
        --PROTHINT_PATH=/opt/ETP/bin/gmes/ProtHint/bin \
        --GENEMARK_PATH=/opt/ETP/bin \
        --workingdir=braker_${id}

    braker_status=\$?

    set -e

    if [[ \$braker_status -ne 0 || ! -s braker_${id}/braker.gff3 ]]; then
        echo "BRAKER2 failed for ${id} on attempt ${task.attempt}/5." >&2

        if [[ ${task.attempt} -lt 5 ]]; then
            if [[ \$braker_status -ne 0 ]]; then
                exit \$braker_status
            fi

            exit 1
        fi

        echo "BRAKER2 failed after five attempts. Creating an empty GFF3 file." >&2
        mkdir -p braker_${id}
        : > braker_${id}/braker.gff3
    fi
    """

    stub:
    """
    command -v braker.pl >/dev/null
    test -d /opt/Augustus/config
    test -x /opt/Augustus/bin/augustus

    mkdir -p braker_${id}
    touch braker_${id}/braker.gff3
    """
}