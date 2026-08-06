

process CD_HIT_EST {
    tag "${id}"
    cpus 50
    time '1d'
    label 'exogap_tools'

    input:
    tuple val(id), path(inputs), val(sequence_identity_threshold), val(length_of_throw_away_sequence), val(alignment_coverage_for_the_shorter_sequence)

    output:
    tuple val(id), path("${id}.nr.fa"), emit: library
    tuple val(id), path("${id}.nr.fa.clstr"), emit: clusters

    script:
    """
    # merge in case there is more than one input file
    cat ${inputs} > ${id}_merged_input.fa
    

    cd-hit-est \
        -i ${id}_merged_input.fa \
        -o ${id}.nr.fa \
        -c ${sequence_identity_threshold} \
        -n 10 \
        -l ${length_of_throw_away_sequence} \
        -aS ${alignment_coverage_for_the_shorter_sequence} \
        -g 1 \
        -sc 1 \
        -sf 1 \
        -d 0 \
        -M 0 \
        -T ${task.cpus}
    """

    stub:
    """
    touch ${id}.nr.fa
    touch ${id}.nr.fa.clstr
    """
}