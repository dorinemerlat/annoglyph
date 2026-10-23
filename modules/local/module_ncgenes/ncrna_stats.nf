process NCRNA_STATS {
    tag 'all_genomes'
    // Publish ncRNA summaries as soon as they are complete, independently of
    // protein-coding annotation selection. The final report uses this layout.
    publishDir "${params.outdir}/comparative/noncoding", mode: params.publish_dir_mode, overwrite: true,
        saveAs: { name -> name == 'tables' ? 'stats' : (name == 'figures' ? 'plots' : name) }
    label 'annoglyph_python'
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
