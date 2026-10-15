process NCRNA_STATS {
    tag 'all_genomes'
    label 'exogap_python'
    memory 8.GB

    input:
    path annotations, stageAs: 'annotations/*'
    path manifest
    path taxonomy
    path tree

    output:
    path 'tables', emit: tables
    path 'figures', emit: figures
    path 'README.txt', emit: readme

    script:
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib
    mkdir -p "\$MPLCONFIGDIR"
    ncrna_stats.py --annotations annotations --manifest ${manifest} \
        --taxonomy ${taxonomy} --tree ${tree} \
        --group-rank '${params.ncrna_stats_group_rank}' \
        --min-group-size ${params.ncrna_stats_min_group_size} --output .
    """

    stub:
    """
    mkdir tables figures
    touch tables/ncrna_counts.tsv tables/trna_isotypes.tsv tables/rrna_types.tsv
    touch figures/trna_isotypes_heatmap.png
    echo 'ncRNA statistics stub' > README.txt
    """
}
