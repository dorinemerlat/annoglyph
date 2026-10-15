process RNAMMER {

    tag "${id}"
    label 'rnammer'
    scratch false

    input:
    tuple val(id), val(meta), path(genome)

    output:
    tuple val(id), val(meta), path("${id}_rnammer.gff"), emit: gff
    tuple val(id), val(meta), path("${id}_rnammer.fa"), emit: fasta

    script:
    """
    set -euo pipefail

    rnammer -S euk -m tsu,ssu,lsu -multi \\
        -gff "${id}_rnammer.gff.tmp" \\
        -f   "${id}_rnammer.fa.tmp" \\
        -h   "${id}_rnammer.report" < "${genome}"

    # Fix GFF: rename 8s_rRNA -> 5s_rRNA
    grep -v '^#' "${id}_rnammer.gff.tmp" \\
    | awk -F'\\t' -v OFS='\\t' '\$9 == "8s_rRNA" { \$9 = "5s_rRNA" } { print }' \\
    | awk -v OFS='\\t' '{
        nb++;
        print \$1,\$2,"gene",\$4,\$5,\$6,\$7,\$8,"ID=rnammer_" nb ";Name=" \$9 ";family=Gene,rRNA";
        print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,"ID=rnammer_" nb "-rRNA;Name=" \$9 ";Parent=rnammer_" nb ";family=Gene,rRNA";
    }' \\
    > "${id}_rnammer.gff"

    # Fix FASTA header molecule tag
    sed "s|molecule=8s_rRNA|molecule=5s_rRNA|g" "${id}_rnammer.fa.tmp" > "${id}_rnammer.fa"

    # Normalize case in all produced files
    sed -i 's/s_rRNA/S_rRNA/g' ${id}_rnammer.*

    """
    
    stub:
    """
    touch ${id}_rnammer.gff
    touch ${id}_rnammer.fa
    """
}