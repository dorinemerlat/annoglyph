process REFORMAT_GENOME {
    tag "$id"
    // publishDir "${params.outdir}/out/$id"
    label 'bioawk'
    memory '10 GB'

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path("genome_${id}.fa"),      emit: fasta
    tuple val(id), val(meta), path("genome_${id}.fa.tsv"),  emit: json

    script:
    def unmodified_genome = "${id}_unreformated.fa"
    def genome_code = meta.mnemonic.exogap
    """
    # Name modification to avoid overwriting the original file (symbolic link)
    mv $genome $unmodified_genome

    # Saving old names (with the new names) in a csv file
    bioawk  -v id=$genome_code -c fastx 'BEGIN{print "previous_IDs", "EXOGAP_IDs"} {print \$name, id"_"NR}' $unmodified_genome > genome_${id}.fa.tsv

    # Modification of genome file
    bioawk -v id=$genome_code -c fastx '{print ">"id"_"NR; print \$seq}' $unmodified_genome |fold -w 60 > genome_${id}.fa
    """

    stub:
    """
    touch genome_${id}.fa
    touch genome_${id}.fa.tsv
    """
}
