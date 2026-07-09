process CD_HIT {
    tag "${id}"
    cpus 50
    time '1d'
    label 'exogap_tools'

    input:
    tuple val(id), val(meta), path(input), val(type), val(sequence_identity_threshold), val(length_of_throw_away_sequence), val(alignment_coverage_for_the_shorter_sequence)

    output:
    tuple val(id), val(meta), path("${id}.nr.fa"), emit: library
    tuple val(id), val(meta), path("${id}.nr.fa.clstr"), emit: clusters

    script:
    """
    if [ "$type" == "nucleotide" ]; then
        command="cd-hit-est"
        word_length="10"

    elif [ "$type" == "protein" ]; then
        command="cd-hit"
        word_length="5"

    else
        echo "Unknown type: $type"
        exit 1
    fi

    # if several files are provided as input, merge them into a single file
    if [ \$(echo $input | wc -w) -gt 1 ]; then
        merged_input="${id}_merged_input.fa"
        cat $input > \$merged_input
        input=\$merged_input
    else
        input=\$input
    fi

    \$command -i \$input \
        -o ${id}.nr.fa  \
        -c $sequence_identity_threshold \
        -n \$word_length \
        -l $length_of_throw_away_sequence \
        -aS $alignment_coverage_for_the_shorter_sequence \
        -g 1 \
        -sc 1 \
        -sf 1 \
        -d 0 \
        -M 0 \
        -T ${task.cpus}

    # c = sequence identity threshold
    # g = length difference cutoff
    # n = word length (default values: 5 for proteins, 10 for nucleotides)
    # l = length of throw_away_sequence (10 for repetitive elements, 10 for others)
    # M = memory limit (in MB) (default 800; 0 for unlimited)
    # T = number of threads
    # aS = alignment coverage for the shorter sequence (98 for repetitive elements, 0 for others)
    """

    stub:
    """
    touch ${id}.nr.fa ${id}.nr.fa.clustr
    """
}
