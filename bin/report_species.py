#!/usr/bin/env python3
"""Add per-genome tables/figures and arrange an export in the public results layout."""
import argparse
import csv
import json
import shutil
from collections import Counter, defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
from exogap_gff import read_gff
from ncrna_stats import isotype, anticodon, rrna_type, rna_kind, save_figure


def coverage(features):
    intervals = defaultdict(list)
    for feature in features:
        intervals[feature.seqid].append((feature.start, feature.end))
    total = 0
    for spans in intervals.values():
        start = end = None
        for left, right in sorted(spans):
            if start is None:
                start, end = left, right
            elif left <= end + 1:
                end = max(end, right)
            else:
                total += end - start + 1
                start, end = left, right
        if start is not None:
            total += end - start + 1
    return total


def table_counter(counter, path, label):
    pd.DataFrame(sorted(counter.items()), columns=[label, 'count']).to_csv(path, sep='\t', index=False)


def bar(counter, path, title, ylabel='Count'):
    fig, ax = plt.subplots(figsize=(max(6, len(counter) * 0.35), 4))
    if counter:
        ax.bar(list(counter), list(counter.values()))
        ax.tick_params(axis='x', labelrotation=90)
    else:
        ax.text(0.5, 0.5, 'No annotations', transform=ax.transAxes, ha='center')
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    save_figure(fig, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--assets', type=Path, required=True)
    args = parser.parse_args()
    root = args.bundle
    metadata = json.loads((root / 'metadata.json').read_text())
    features = read_gff(root / 'annotation.gff3')
    for category in ['genome', 'repeats', 'protein_coding', 'noncoding']:
        for section in ['sequences', 'stats', 'plots']:
            (root / category / section).mkdir(parents=True, exist_ok=True)
    (root / 'submission').mkdir()
    (root / 'protein_coding' / 'selection').mkdir()
    paths = {'annotation.gff3': 'genome/annotation.gff3', 'genome.masked.fna': 'genome/sequences/genome.softmasked.fna',
             'cds.fna': 'protein_coding/sequences/cds.fna', 'proteins.faa': 'protein_coding/sequences/proteins.faa',
             'ncrna.fna': 'noncoding/sequences/ncrna.fna', 'genome.fsa': 'submission/genome.fsa',
             'annotation.tbl': 'submission/annotation.tbl', 'metadata.json': 'submission/metadata.json',
             'id_mapping.tsv': 'submission/id_mapping.tsv', 'validation.tsv': 'submission/validation.tsv'}
    for source, destination in paths.items():
        shutil.move(root / source, root / destination)
    (root / 'README.txt').unlink(missing_ok=True)
    (root / 'checksums.sha256').unlink(missing_ok=True)
    sizes, gc, acgt, gaps, softmasked = [], 0, 0, 0, 0
    with (root / 'genome/sequences/genome.fna').open('w') as handle:
        for record in SeqIO.parse(root / 'genome/sequences/genome.softmasked.fna', 'fasta'):
            seq = str(record.seq)
            upper = seq.upper()
            sizes.append(len(seq))
            gc += upper.count('G') + upper.count('C')
            acgt += sum(upper.count(base) for base in 'ACGT')
            gaps += upper.count('N')
            softmasked += sum(base.islower() for base in seq)
            handle.write('>' + record.id + '\n')
            for i in range(0, len(seq), 60):
                handle.write(upper[i:i + 60] + '\n')
    total = sum(sizes)
    cumulative = 0
    n50 = l50 = 0
    for index, size in enumerate(sorted(sizes, reverse=True), 1):
        cumulative += size
        if cumulative >= total / 2:
            n50, l50 = size, index
            break
    genome_stats = {'genome_size_bp': total, 'contigs': len(sizes), 'N50_bp': n50, 'L50': l50,
                    'GC_percent': gc / acgt * 100 if acgt else None, 'N_bp': gaps, 'softmasked_bp': softmasked}
    pd.DataFrame([genome_stats]).to_csv(root / 'genome/stats/assembly.tsv', sep='\t', index=False)
    fig, ax = plt.subplots()
    ax.hist(sizes, bins=min(40, max(1, len(sizes))))
    ax.set(xlabel='Contig length (bp)', ylabel='Contigs', title=metadata['organism'])
    save_figure(fig, root / 'genome/plots/contig_lengths')
    groups = {category: [f for f in features if f.attrs['ID'][0].startswith(metadata['sample'] + '.' + prefix + '.')]
              for category, prefix in [('protein_coding', 'coding'), ('repeats', 'repeats'), ('noncoding', 'ncrna')]}
    for category, records in groups.items():
        with (root / category / 'annotation.gff3').open('w') as handle:
            handle.write('##gff-version 3\n')
            for feature in records:
                handle.write(feature.line() + '\n')
    repeats = groups['repeats']
    repeat_rows = []
    for level in ['repeat_type', 'repeat_class', 'repeat_order', 'repeat_superfamily', 'repeat_family']:
        classified = defaultdict(list)
        for feature in repeats:
            classified[feature.attrs.get(level, ['Unknown'])[0]].append(feature)
        for name, records in sorted(classified.items()):
            covered = coverage(records)
            repeat_rows.append([level, name, len(records), covered, covered / total * 100 if total else 0])
    pd.DataFrame(repeat_rows, columns=['level', 'category', 'count', 'covered_bp', 'coverage_percent']).to_csv(
        root / 'repeats/stats/composition.tsv', sep='\t', index=False)
    repeat_counts = Counter(f.attrs.get('repeat_type', ['Unknown'])[0] for f in repeats)
    bar(repeat_counts, root / 'repeats/plots/composition', 'Repeat annotation counts')
    coding = groups['protein_coding']
    feature_counts = Counter(f.kind for f in coding)
    table_counter(feature_counts, root / 'protein_coding/stats/features.tsv', 'feature')
    bar(feature_counts, root / 'protein_coding/plots/features', 'Protein-coding annotation features')
    lengths = [len(record.seq) for record in SeqIO.parse(root / 'protein_coding/sequences/proteins.faa', 'fasta')]
    pd.DataFrame({'protein_length_aa': lengths}).to_csv(root / 'protein_coding/stats/protein_lengths.tsv', sep='\t', index=False)
    fig, ax = plt.subplots()
    if lengths:
        ax.hist(lengths, bins=40)
    ax.set(xlabel='Protein length (aa)', ylabel='Proteins', title='Selected proteins')
    save_figure(fig, root / 'protein_coding/plots/protein_lengths')
    pd.DataFrame([{'genome_id': metadata['sample'], 'method_annotation': metadata['selected_method']}]).to_csv(
        root / 'protein_coding/selection/selected_method.tsv', sep='\t', index=False)
    kinds, trnas, rrnas, anticodons, seen = Counter(), Counter(), Counter(), Counter(), set()
    for feature in groups['noncoding']:
        kind = rna_kind(feature)
        identifier = feature.attrs['ID'][0]
        if kind is None or identifier in seen:
            continue
        seen.add(identifier)
        kinds[kind] += 1
        if kind == 'tRNA':
            trnas[isotype(feature)] += 1
            codon = anticodon(feature)
            if codon:
                anticodons[isotype(feature) + ':' + codon] += 1
        elif kind == 'rRNA':
            rrnas[rrna_type(feature)] += 1
    for counter, name, label in [(kinds, 'rna_types', 'RNA_type'), (trnas, 'trna_isotypes', 'isotype'),
                                  (rrnas, 'rrna_types', 'rRNA_type'), (anticodons, 'trna_anticodons', 'isotype_anticodon')]:
        table_counter(counter, root / 'noncoding/stats' / (name + '.tsv'), label)
        if counter or name != 'trna_anticodons':
            bar(counter, root / 'noncoding/plots' / name, label)
    summary = {**genome_stats, 'genome_id': metadata['sample'], 'organism': metadata['organism'],
               'repeat_count': len(repeats), 'repeat_covered_bp': coverage(repeats),
               'repeat_coverage_percent': coverage(repeats) / total * 100 if total else 0,
               'protein_count': len(lengths), 'ncrna_count': sum(kinds.values()),
               'trna_count': kinds['tRNA'], 'rrna_count': kinds['rRNA']}
    pd.DataFrame([summary]).to_csv(root / 'genome/stats/overview.tsv', sep='\t', index=False)
    for asset in json.loads(args.assets.read_text()):
        destination = root / asset['destination']
        destination.parent.mkdir(parents=True, exist_ok=True)
        if Path(asset['source']).is_dir():
            shutil.copytree(asset['source'], destination, dirs_exist_ok=True)
        else:
            shutil.copyfile(asset['source'], destination)
    # Keep selected protein and assembly BUSCO percentages visible per species.
    for category, directory in [('genome', 'genome/stats/busco'), ('protein_coding', 'protein_coding/stats/quality')]:
        for source in (root / directory).glob('*.json'):
            result = json.loads(source.read_text()).get('results', {})
            values = {label: result.get(key) for label, key in [('Single', 'Single copy percentage'),
                      ('Duplicated', 'Multi copy percentage'), ('Fragmented', 'Fragmented percentage'), ('Missing', 'Missing percentage')]}
            values = {key: float(value) for key, value in values.items() if value is not None}
            if values:
                bar(values, root / category / 'plots/busco', 'BUSCO assessment', 'Percentage')


if __name__ == '__main__':
    main()
