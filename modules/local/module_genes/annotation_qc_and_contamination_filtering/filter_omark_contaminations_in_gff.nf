process FILTER_OMARK_CONTAMINATIONS_IN_GFF {
    tag "${annotation_method}/${id}"
    label 'exogap_python'
    memory '20 GB'
    scratch false
    
    input:
    tuple val(id), val(meta), val(annotation_method), path(gff), path(contaminations)

    output:
    tuple val(id), val(meta), val(annotation_method), path("${id}_${annotation_method}_contaminations-filtered-c-${threshold}.gff")

    script:
    def contaminations_filtered = "${id}_${annotation_method}_contaminations-filtered-c-${threshold}"
    """
    # generation of filtered gff
    awk '
    /^-- Contaminant genes --[[:space:]]*\$/ {
        in_section=1
        skip_header=1
        next
    }

    /^-- .* --[[:space:]]*\$/ && in_section {
        in_section=0
        next
    }

    in_section {
        if (skip_header) { skip_header=0; next }
        if (\$1 != "") print \$1
    }
    ' ${contaminations} | sort -u > contaminant_genes.txt

    agat_sp_filter_feature_from_keep_list.pl --gff $gff --keep_list contaminant_genes.txt --output ${contaminations_filtered}.gff
    """

    stub:
    """
    touch ${id}_${annotation_method}_contaminations-filtered-c-${threshold}.gff
    """
}