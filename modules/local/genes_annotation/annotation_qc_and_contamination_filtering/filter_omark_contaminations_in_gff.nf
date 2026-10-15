process FILTER_OMARK_CONTAMINATIONS_IN_GFF {
    tag "${meta.annotation_method}/${meta.id}"
    label 'agat'
    label 'omark_filter'

    input:
    tuple val(meta), path(contaminations), path(gff)

    output:
    tuple val(meta), path("${meta.id}_${meta.annotation_method}_contaminations-filtered.gff")

    script:
    """
    # Extract contaminant gene identifiers from the OMArk report.
    awk '
    /^-- Contaminant genes --[[:space:]]*\$/ {
        in_section = 1
        skip_header = 1
        next
    }

    /^-- .* --[[:space:]]*\$/ && in_section {
        in_section = 0
        next
    }

    in_section {
        if (skip_header) {
            skip_header = 0
            next
        }

        if (\$1 != "")
            print \$1
    }
    ' ${contaminations} | sort -u > contaminant_genes.txt

    agat_sp_filter_feature_from_kill_list.pl \
        --gff ${gff} \
        --kill_list contaminant_genes.txt \
        --output ${meta.id}_${meta.annotation_method}_contaminations-filtered.gff
    """

    stub:
    """
    touch ${meta.id}_${meta.annotation_method}_contaminations-filtered.gff
    """
}