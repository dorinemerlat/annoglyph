process REFORMAT_GFF3 {
    tag "${id}/${annotation_method}"

    input:
    tuple val(id), val(meta), val(annotation_method), path(gff)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_filtered.gff")

    script:
    """
    awk -v prefix="${meta.mnemonic.for_exogap}" 'BEGIN{OFS="\t"}
    {
        if (\$0 ~ /^#/) {
            print \$0
        } else {
            gsub(/ID=/, "ID=" prefix ".", \$9)
            gsub(/Parent=/, "Parent=" prefix ".", \$9)
            print \$0
        }
    }' $gff > ${id}_${annotation_method}_filtered.gff
    """

    stub:
    """
    touch ${id}_${annotation_method}_filtered.gff
    """
}