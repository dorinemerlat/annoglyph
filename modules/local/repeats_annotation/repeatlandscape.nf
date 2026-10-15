process REPEATLANDSCAPE {
    tag "${meta.id}"
    label 'repeatmasker'

    input:
    tuple val(meta), path(align)

    output:
    tuple val(meta), path("${meta.id}_TE_repeat_landscape.tsv"), emit: tsv
    tuple val(meta), path("${align}.divsum"), emit: divsum

    script:
    """
    # Allow RepeatMasker Perl modules to be found
    export PERL5LIB=/usr/local/share/RepeatMasker/

    # Recalculate divergence based on GC content
    calcDivergenceFromAlign.pl -s ${align}.divsum -a ${align}.with_div ${align}

    # Create the repeat landscape
    process_divergence_from_align.py \
        -d ${align}.divsum \
        -s ${meta.assembly_size} \
        -o ${meta.id}_TE_repeat_landscape.tsv
    """

    stub:
    """
    touch ${meta.id}_TE_repeat_landscape.tsv
    touch ${align}.divsum
    """
}