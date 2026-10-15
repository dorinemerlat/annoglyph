#!/usr/bin/env python3
"""Collect per-species summaries without rereading the genome sequences."""
import argparse
import json
import shutil
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from ncrna_stats import phylogeny, save_figure


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bundles', type=Path, required=True)
    parser.add_argument('--taxonomy', type=Path, required=True)
    parser.add_argument('--tree', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--assets', type=Path, required=True)
    args = parser.parse_args()
    for category in ['genome', 'repeats', 'protein_coding']:
        for section in ['stats', 'plots']:
            (args.output / category / section).mkdir(parents=True, exist_ok=True)
    summaries = [pd.read_csv(path, sep='\t') for path in args.bundles.glob('*/genome/stats/overview.tsv')]
    if not summaries:
        raise ValueError('No species summary found')
    summary = pd.concat(summaries, ignore_index=True).set_index('genome_id', verify_integrity=True)
    taxonomy = pd.read_csv(args.taxonomy, sep='\t', dtype={'id': str}).set_index('id', verify_integrity=True)
    tree, leaves, order, unmatched = phylogeny(args.tree, taxonomy, list(summary.index))
    summary = summary.loc[order]
    shutil.copyfile(args.taxonomy, args.output / 'taxonomy.tsv')
    shutil.copyfile(args.tree, args.output / 'taxonomy_tree.nwk')
    for category, columns in {
        'genome': ['genome_size_bp', 'contigs', 'N50_bp', 'L50', 'GC_percent', 'N_bp', 'softmasked_bp'],
        'repeats': ['repeat_count', 'repeat_covered_bp', 'repeat_coverage_percent'],
        'protein_coding': ['protein_count'],
    }.items():
        summary[columns].T.to_csv(args.output / category / 'stats/summary.tsv', sep='\t')
        fig, axes = plt.subplots(len(columns), 1, figsize=(max(8, len(order) * 0.5), len(columns) * 2.3), squeeze=False)
        for ax, column in zip(axes[:, 0], columns):
            ax.bar(order, summary[column])
            ax.set_ylabel(column)
            ax.tick_params(axis='x', labelrotation=90)
        fig.tight_layout()
        save_figure(fig, args.output / category / 'plots/overview')
    for category, relative, label in [('repeats', 'repeats/stats/composition.tsv', 'repeat_type'),
                                       ('protein_coding', 'protein_coding/stats/features.tsv', 'feature')]:
        rows = []
        for bundle in args.bundles.iterdir():
            path = bundle / relative
            if not path.exists():
                continue
            table = pd.read_csv(path, sep='\t')
            if category == 'repeats':
                table = table[table.level == 'repeat_type'].rename(columns={'category': label})
            for row in table.to_dict('records'):
                rows.append({'genome_id': bundle.name, 'type': row[label], 'count': row['count']})
        if rows:
            matrix = pd.DataFrame(rows).pivot_table(index='type', columns='genome_id', values='count', aggfunc='sum', fill_value=0).reindex(columns=order, fill_value=0)
            matrix.to_csv(args.output / category / 'stats/counts.tsv', sep='\t')
            fig, ax = plt.subplots(figsize=(max(8, len(order) * 0.5), max(4, len(matrix) * 0.3)))
            sns.heatmap(matrix, ax=ax, cmap='viridis', cbar_kws={'label': 'Feature rows'})
            save_figure(fig, args.output / category / 'plots/counts_heatmap')
    for asset in json.loads(args.assets.read_text()):
        destination = args.output / asset['destination']
        destination.parent.mkdir(parents=True, exist_ok=True)
        if Path(asset['source']).is_dir():
            shutil.copytree(asset['source'], destination, dirs_exist_ok=True)
        else:
            shutil.copyfile(asset['source'], destination)
    selections = [pd.read_csv(path, sep='\t') for path in args.bundles.glob('*/protein_coding/selection/selected_method.tsv')]
    pd.concat(selections, ignore_index=True).to_csv(args.output / 'protein_coding/stats/best_annotations.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
