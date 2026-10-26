#!/usr/bin/env python3
"""Summarize reconciled ncRNAs and plot genome comparisons in taxonomy-tree order."""
import argparse
import hashlib
import re
from collections import Counter, defaultdict
from io import StringIO
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Phylo
from scipy.stats import pearsonr, spearmanr

from annoglyph_gff import read_gff

ISOTYPES = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile',
            'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val',
            'Sec', 'Sup', 'Undet', 'Unknown']
RRNA_TYPES = ['5S', '5.8S', '12S', '16S', '18S', '23S', '25S', '26S', '28S', 'SSU', 'LSU', 'Unknown']


def attribute_text(feature):
    return ' '.join(value for key, values in feature.attrs.items()
                    if key.lower() in ('isotype', 'amino_acid', 'product', 'name', 'note', 'model')
                    for value in values)


def isotype(feature):
    text = attribute_text(feature)
    for candidate in ISOTYPES[:-1]:
        if re.search(r'(?<![A-Za-z])' + candidate + r'(?![a-z])', text, re.I):
            return candidate
    return 'Unknown'


def anticodon(feature):
    for key, values in feature.attrs.items():
        if key.lower() in ('anticodon', 'anticodon_seq', 'anticodon_sequence'):
            for value in values:
                match = re.fullmatch(r'[ACGTU]{3}', value.strip(), re.I)
                if not match:
                    match = re.search(r'(?:seq[:=])\s*([ACGTU]{3})\b', value, re.I)
                if match:
                    return (match.group(1) if match.lastindex else match.group()).upper().replace('U', 'T')
    text = attribute_text(feature)
    for expression in [r'\bPredicted\s+([ACGTU]{3})\s+anticodon\b',
                       r'\banticodon\s*[:=]?\s*([ACGTU]{3})\b',
                       r'(?:' + '|'.join(ISOTYPES[:-1]) + r')[-_]?([ACGTU]{3})\b']:
        match = re.search(expression, text, re.I)
        if match:
            return match.group(1).upper().replace('U', 'T')
    return None


def rrna_type(feature):
    text = re.sub(r'5[_ .]8S', '5.8S', attribute_text(feature), flags=re.I)
    for candidate in ['5.8S', '28S', '26S', '25S', '23S', '18S', '16S', '12S', '5S', 'SSU', 'LSU']:
        if re.search(r'(?<![\d.])' + re.escape(candidate) + r'(?![\d.])', text, re.I):
            return candidate
    return 'Unknown'


def rna_kind(feature):
    if feature.kind in ('tRNA', 'rRNA'):
        return feature.kind
    if feature.kind.endswith('RNA') and feature.kind != 'mRNA':
        if feature.attrs.get('rna_type'):
            return feature.attrs['rna_type'][0]
        text = ' '.join(feature.attrs.get('rfam_type', []))
        for kind in ['snoRNA', 'snRNA', 'miRNA', 'lncRNA', 'ribozyme', 'antisense_RNA']:
            if kind in text or feature.kind == kind:
                return kind
        if 'Cis-reg' in text:
            return 'Regulatory RNA'
        return 'Unclassified ncRNA'
    return None


def save_figure(fig, stem):
    for extension in ['png', 'pdf', 'svg']:
        fig.savefig(str(stem) + '.' + extension, dpi=200, bbox_inches='tight')
    plt.close(fig)


def phylogeny(tree_path, taxonomy, genomes):
    text = tree_path.read_text()
    # The existing taxonomy exporter writes unquoted scientific names with spaces.
    def quote_label(match):
        label = match.group(1).strip()
        if ' ' in label and not label.startswith("'"):
            return "'" + label.replace("'", "''") + "'"
        return label
    text = re.sub(r'(?<=[(),])([^():;,]+)(?=[:;,)])', quote_label, text)
    tree = Phylo.read(StringIO(text), 'newick')
    aliases = defaultdict(list)
    for genome in genomes:
        row = taxonomy.loc[genome]
        for label in {genome, str(row.get('name', '')), str(row.get('taxonomy', '')), str(row.get('species', ''))}:
            if label and label != 'nan':
                aliases[label].append(genome)
    leaf_genomes = {}
    order = []
    for leaf in tree.get_terminals():
        label = str(leaf.name or '').split('|')[0]
        matches = sorted(set(aliases.get(label, [])))
        leaf_genomes[id(leaf)] = matches
        order.extend(genome for genome in matches if genome not in order)
    unmatched = sorted(set(genomes) - set(order))
    order.extend(unmatched)
    return tree, leaf_genomes, order, unmatched


def draw_tree(ax, tree, leaf_genomes, order):
    positions = {genome: index + 0.5 for index, genome in enumerate(order)}
    depths = tree.depths(unit_branch_lengths=True)
    coordinates = {}
    for clade in tree.find_clades(order='postorder'):
        if clade.is_terminal():
            values = [positions[g] for g in leaf_genomes[id(clade)]]
            if values:
                coordinates[id(clade)] = np.mean(values)
        else:
            values = [coordinates[id(child)] for child in clade.clades if id(child) in coordinates]
            if values:
                x = float(np.mean(values))
                coordinates[id(clade)] = x
                y = depths[clade]
                ax.plot([min(values), max(values)], [y, y], color='black', lw=0.7)
                for child in clade.clades:
                    if id(child) in coordinates:
                        cx = coordinates[id(child)]
                        ax.plot([cx, cx], [y, depths[child]], color='black', lw=0.7)
    ax.set_xlim(0, len(order))
    ax.invert_yaxis()
    ax.axis('off')
    ax.set_title('Taxonomic tree order', fontsize=10)


def heatmap(matrix, stem, title, labels, tree, leaf_genomes, order):
    width = max(7, 0.65 * len(order))
    height = max(5, 0.25 * len(matrix) + 2.5)
    fig = plt.figure(figsize=(width, height))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.3, max(3, 0.25 * len(matrix))], width_ratios=[1, 0.035], hspace=0.12)
    tree_ax = fig.add_subplot(grid[0, 0])
    draw_tree(tree_ax, tree, leaf_genomes, order)
    ax, color_ax = fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])
    sns.heatmap(matrix, ax=ax, cbar_ax=color_ax, cmap='viridis', vmin=0,
                xticklabels=[labels[g] for g in order], yticklabels=matrix.index,
                cbar_kws={'label': 'Number of loci'})
    ax.set_title(title)
    ax.set_xlabel('Genome')
    ax.set_ylabel('')
    ax.tick_params(axis='x', labelrotation=90)
    ax.tick_params(axis='y', labelrotation=0, labelsize=8)
    save_figure(fig, stem)


def correlation_groups(taxonomy, order, tree, leaf_genomes, rank, minimum):
    groups = [('all', 'all', order)]
    for group, frame in taxonomy.loc[order].groupby(rank, dropna=False, sort=True):
        members = list(frame.index)
        groups.append(('taxonomy:' + rank, str(group) if pd.notna(group) else 'Unassigned', members))
    seen = {tuple(sorted(order))}
    for index, clade in enumerate(tree.find_clades()):
        members = sorted({genome for leaf in clade.get_terminals() for genome in leaf_genomes[id(leaf)]})
        key = tuple(members)
        if len(members) < minimum or key in seen:
            continue
        seen.add(key)
        groups.append(('tree', str(clade.name or f'clade_{index}'), members))
    return groups


def correlations(metrics, taxonomy, groups, rank, minimum, output):
    rows, membership = [], []
    for source, group, members in groups:
        slug = re.sub(r'[^A-Za-z0-9_-]+', '_', source + '_' + group).strip('_')[:90]
        slug += '_' + hashlib.sha256((source + group).encode()).hexdigest()[:8]
        for genome in members:
            membership.append([source, group, genome])
        frame = metrics.loc[members]
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
        for ax, metric, label in zip(axes, ['trna_count', 'trna_per_kb'], ['tRNA loci', 'tRNA loci per kb']):
            x, y = frame['genome_size_bp'] / 1e6, frame[metric]
            valid = len(frame) >= minimum and x.nunique() > 1 and y.nunique() > 1
            reason = 'ok' if valid else ('insufficient_genomes' if len(frame) < minimum else 'constant_values')
            pearson_r = pearson_p = spearman_r = spearman_p = np.nan
            if valid:
                pearson_r, pearson_p = pearsonr(x, y)
                spearman_r, spearman_p = spearmanr(x, y)
                slope, intercept = np.polyfit(x, y, 1)
                domain = np.linspace(x.min(), x.max(), 100)
                ax.plot(domain, intercept + slope * domain, color='gray', lw=1)
            categories = taxonomy.loc[members, rank].fillna('Unassigned')
            palette = dict(zip(sorted(categories.unique()), sns.color_palette('tab20', n_colors=categories.nunique())))
            for category, color in palette.items():
                subset = categories[categories == category].index
                ax.scatter(x.loc[subset], y.loc[subset], label=category, color=color)
            if len(frame) <= 20:
                for genome in members:
                    ax.annotate(str(taxonomy.loc[genome, 'name']), (x[genome], y[genome]),
                                xytext=(4, 4), textcoords='offset points', fontsize=6)
            annotation = f'n={len(frame)}\nPearson r={pearson_r:.3f}\nSpearman rho={spearman_r:.3f}' if valid else f'n={len(frame)}; {reason}'
            ax.text(0.02, 0.98, annotation, transform=ax.transAxes, va='top', fontsize=8)
            ax.set_xlabel('Genome size (Mb)')
            ax.set_ylabel(label)
            ax.legend(title=rank, fontsize=6, title_fontsize=7)
            rows.append([source, group, metric, len(frame), pearson_r, pearson_p, spearman_r, spearman_p, reason])
        fig.suptitle(f'{source} / {group}')
        fig.tight_layout()
        save_figure(fig, output / 'figures' / 'correlations' / slug)
    pd.DataFrame(rows, columns=['group_source', 'group', 'metric', 'n_genomes', 'pearson_r',
                               'pearson_p', 'spearman_rho', 'spearman_p', 'status']).to_csv(
        output / 'tables' / 'correlations.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(membership, columns=['group_source', 'group', 'genome_id']).to_csv(
        output / 'tables' / 'group_membership.tsv', sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--annotations', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--taxonomy', type=Path, required=True)
    parser.add_argument('--tree', type=Path, required=True)
    parser.add_argument('--group-rank', default='family')
    parser.add_argument('--min-group-size', type=int, default=3)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.min_group_size < 3:
        parser.error('min-group-size must be at least 3')
    manifest = pd.read_csv(args.manifest, sep='\t', dtype={'genome_id': str})
    if manifest.empty or manifest.genome_id.duplicated().any():
        raise ValueError('Expected a nonempty manifest with one annotation per genome')
    if (manifest.genome_size_bp <= 0).any():
        raise ValueError('Genome lengths must be positive')
    taxonomy = pd.read_csv(args.taxonomy, sep='\t', dtype={'id': str}).set_index('id', verify_integrity=True)
    if args.group_rank not in taxonomy.columns:
        raise ValueError(f'Taxonomy rank not found: {args.group_rank}')
    genomes = manifest.genome_id.tolist()
    missing = set(genomes) - set(taxonomy.index)
    if missing:
        raise ValueError(f'Missing taxonomy for genomes: {sorted(missing)}')
    tree, leaf_genomes, order, unmatched = phylogeny(args.tree, taxonomy, genomes)
    for directory in ['tables', 'figures', 'figures/correlations']:
        (args.output / directory).mkdir(parents=True, exist_ok=True)
    counts, trna_counts, rrna_counts, ac_counts = {}, {}, {}, {}
    loci, metric_rows = [], []
    for row in manifest.itertuples(index=False):
        genome = row.genome_id
        features = read_gff(args.annotations / row.annotation)
        broad, trna, rrna, ac = Counter(), Counter(), Counter(), Counter()
        seen = set()
        for feature in features:
            kind = rna_kind(feature)
            if kind is None:
                continue
            # Multipart GFF rows describe one locus, not several genes.
            identifier = feature.attrs.get('ID', [f'{feature.seqid}:{feature.start}-{feature.end}:{feature.strand}:{kind}'])[0]
            if identifier in seen:
                continue
            seen.add(identifier)
            broad[kind] += 1
            subtype, codon = '', None
            if kind == 'tRNA':
                subtype, codon = isotype(feature), anticodon(feature)
                trna[subtype] += 1
                if codon:
                    ac[subtype + ':' + codon] += 1
            elif kind == 'rRNA':
                subtype = rrna_type(feature)
                rrna[subtype] += 1
            loci.append([genome, identifier, kind, subtype, codon or 'NA', feature.seqid,
                         feature.start, feature.end, feature.strand, feature.source])
        counts[genome], trna_counts[genome], rrna_counts[genome], ac_counts[genome] = broad, trna, rrna, ac
        metric_rows.append([genome, row.genome_size_bp, sum(broad.values()), broad['tRNA'],
                            broad['tRNA'] / (row.genome_size_bp / 1000), broad['rRNA'],
                            sum(ac.values()), broad['tRNA'] - sum(ac.values())])
    metrics = pd.DataFrame(metric_rows, columns=['genome_id', 'genome_size_bp', 'ncrna_count',
                            'trna_count', 'trna_per_kb', 'rrna_count', 'trna_with_anticodon', 'trna_without_anticodon']).set_index('genome_id').loc[order]

    def matrix(data, rows, filename):
        frame = pd.DataFrame({genome: {key: data[genome][key] for key in rows} for genome in order}, index=rows).astype(int)
        frame.index.name = 'type'
        frame.to_csv(args.output / 'tables' / filename, sep='\t')
        return frame

    kinds = ['tRNA', 'rRNA'] + sorted({kind for counter in counts.values() for kind in counter} - {'tRNA', 'rRNA'})
    summary = matrix(counts, kinds, 'ncrna_counts.tsv')
    summary.loc['Total'] = summary.sum(axis=0)
    summary.to_csv(args.output / 'tables' / 'ncrna_counts.tsv', sep='\t')
    trna_matrix = matrix(trna_counts, ISOTYPES, 'trna_isotypes.tsv')
    rrna_matrix = matrix(rrna_counts, RRNA_TYPES, 'rrna_types.tsv')
    ac_rows = sorted({key for counter in ac_counts.values() for key in counter})
    ac_matrix = matrix(ac_counts, ac_rows, 'trna_anticodons.tsv')
    metrics.join(taxonomy, how='left').to_csv(args.output / 'tables' / 'genome_metrics.tsv', sep='\t')
    pd.DataFrame(loci, columns=['genome_id', 'locus_id', 'RNA_type', 'subtype', 'anticodon',
                              'seqid', 'start', 'end', 'strand', 'source']).to_csv(
        args.output / 'tables' / 'ncrna_loci.tsv', sep='\t', index=False)
    labels = {genome: str(taxonomy.loc[genome, 'name']) + '\n' + genome for genome in order}
    heatmap(summary.drop(index='Total'), args.output / 'figures' / 'ncrna_types_heatmap',
            'Non-coding RNA types', labels, tree, leaf_genomes, order)
    heatmap(trna_matrix, args.output / 'figures' / 'trna_isotypes_heatmap', 'tRNA isotypes', labels, tree, leaf_genomes, order)
    heatmap(rrna_matrix, args.output / 'figures' / 'rrna_types_heatmap', 'rRNA types', labels, tree, leaf_genomes, order)
    if ac_rows:
        heatmap(ac_matrix, args.output / 'figures' / 'trna_anticodons_heatmap', 'tRNA isotype / anticodon', labels, tree, leaf_genomes, order)
    groups = correlation_groups(taxonomy, order, tree, leaf_genomes, args.group_rank, args.min_group_size)
    correlations(metrics, taxonomy, groups, args.group_rank, args.min_group_size, args.output)
    (args.output / 'README.txt').write_text(
        'Counts use retained RNA loci only, excluding gene/exon parent/child rows.\n'
        'All count matrices have genomes in columns, in taxonomy-tree order.\n'
        'Anticodons use 5-prime to 3-prime DNA notation (U normalized to T); missing values are not inferred.\n'
        'Density = retained tRNA count / (total assembly length in bp / 1000).\n'
        'Correlations are ordinary descriptive Pearson/Spearman comparisons, not phylogenetically corrected.\n'
        'Density uses genome size in its denominator: correlation with genome size is mathematically coupled.\n'
        'P-values are unadjusted exploratory values; overlapping clades are not independent tests.\n'
        f'Correlations need at least {args.min_group_size} genomes and nonconstant variables.\n'
        'No anticodon heatmap is produced when no retained tRNA has an annotated anticodon.\n'
        f'Genomes unmatched in tree (appended in ID order): {", ".join(unmatched) or "none"}.\n')


if __name__ == '__main__':
    main()
