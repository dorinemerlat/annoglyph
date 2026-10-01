process BARRNAP {

    tag "${id}"
    label 'exogap_tools'
    cpus 30

    input:
    tuple val(id), val(meta), path(genome), val(kingdom)

    output:
    tuple val(id), val(meta), path("${id}_barrnap_${kingdom}.gff"), emit: gff

    script:
    """
    barrnap --quiet --kingdom ${kingdom} --threads ${task.cpus} "${genome}"  > barrnap.out

    # Normalize IDs and expand to gene + rRNA feature rows
    if [ \$(wc -l < barrnap.out) -gt 1 ]; then
        grep -v '^#' barrnap.out \
        | sed -E 's/barrnap:[0-9]+(\\.[0-9]+)?/barrnap/g' \
        | awk -v OFS='\\t' '{
            nb++;
            # gene feature
            print \$1,\$2,"gene",\$4,\$5,\$6,\$7,\$8,"ID=barrnap_" nb ";" \$9 ";family=Gene,rRNA";
            # child rRNA feature
            print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,"ID=rnammer_" nb "-rRNA;Parent=barrnap_" nb ";" \$9 ";family=Gene,rRNA";
        }' > "${id}_barrnap_${kingdom}.gff"
    else 
        echo "##gff-version 3" > "${id}_barrnap_${kingdom}.gff"
    fi
    """

    stub:
    """
    touch ${id}_barrnap_${kingdom}.gff
    """
}