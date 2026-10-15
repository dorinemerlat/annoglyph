

process CD_HIT_EST {
    tag "${meta.id}"
    label 'cdhit'

    input:
    tuple val(meta), path(inputs), val(sequence_identity_threshold), val(length_of_throw_away_sequence), val(alignment_coverage_for_the_shorter_sequence)

    output:
    tuple val(meta), path("${meta.id}.nr.fa"), emit: library
    tuple val(meta), path("${meta.id}.nr.fa.clstr"), emit: clusters

    script:
    """
    # merge in case there is more than one input file
    cat ${inputs} > ${meta.id}_merged_input.fa
    

    cd-hit-est \
        -i ${meta.id}_merged_input.fa \
        -o ${meta.id}.nr.fa \
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
    touch ${meta.id}.nr.fa
    touch ${meta.id}.nr.fa.clstr
    """
}