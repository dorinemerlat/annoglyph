process REPEATLANDSCAPE {
    tag "${id}"
    label 'repeatmasker'

    input:
    tuple val(id), val(meta), path(align)

    output:
    tuple val(id), val(meta), path("${id}_TE_repeat_landscape.tsv"), emit: tsv
    // tuple val(id), val(meta), path("${align}.divsum"), emit: divsum

    script:
    """
    # Allow to find calcDivergenceFromAlign.pl
    export PERL5LIB=/usr/local/share/RepeatMasker/

    # Re-calculate divergence adapted to GC content
    calcDivergenceFromAlign.pl -s ${align}.divsum -a ${align}.with_div $align

    # Create the repeat landscape
    process_divergence_from_align.py -d ${align}.divsum -s ${meta.assembly_size} -o ${id}_TE_repeat_landscape.tsv
    """

    stub:
    """
    touch ${id}_TE_repeat_landscape.tsv ${align}.divsum
    """
}
