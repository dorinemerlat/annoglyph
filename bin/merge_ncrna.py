#!/usr/bin/env python3
"""Reconcile tRNAscan-SE, Barrnap and Infernal into one GFF3 per genome."""
import argparse
import csv
import re
from bisect import bisect_right
from collections import defaultdict
from annoglyph_gff import Feature, namespace, read_gff, sort_key


def interval_index(features, kind):
    grouped = defaultdict(list)
    for f in features:
        if f.kind.lower() == kind.lower():
            grouped[f.seqid, f.strand].append((f.start, f.end))
    result = {}
    for key, intervals in grouped.items():
        merged = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1] + 1:
                merged[-1][1] = max(end, merged[-1][1])
            else:
                merged.append([start, end])
        result[key] = ([v[0] for v in merged], [v[1] for v in merged])
    return result


def overlaps(feature, index):
    starts, ends = index.get((feature.seqid, feature.strand), ([], []))
    position = bisect_right(starts, feature.end) - 1
    return position >= 0 and ends[position] >= feature.start


def family_types(path):
    result = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if not fields[0].startswith('RF'):
                continue
            if len(fields) < 19:
                raise ValueError('Expected Rfam database_files/family.txt with at least 19 columns')
            result[fields[0]] = fields[18]
    return result


def infernal_features(path, types):
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.split(maxsplit=27)
            if len(fields) < 27:
                raise ValueError(f'{path}:{number}: expected cmscan --fmt 2 output')
            # --oclan --oskip suppresses lower-scoring same-clan overlaps.
            if fields[19] == '=':
                continue
            accession, model = fields[2], fields[1]
            if accession not in types:
                raise ValueError(f'Missing Rfam family metadata for {accession}')
            family_type = types[accession]
            kind = 'tRNA' if re.search(r'\btRNA\b', family_type) else (
                'rRNA' if re.search(r'\brRNA\b', family_type) else 'ncRNA')
            start, end = sorted((int(fields[9]), int(fields[10])))
            attrs = {'Name': [model], 'Dbxref': [f'Rfam:{accession}'],
                     'rfam_type': [family_type], 'evalue': [fields[17]],
                     'model': [model]}
            if kind == 'ncRNA':
                attrs['ncRNA_class'] = [family_type.split(';')[-1].strip() or 'other']
            yield Feature(fields[3], 'Infernal', kind, start, end,
                          fields[16], fields[11], '.', attrs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trnascan', required=True)
    parser.add_argument('--barrnap', nargs='+', required=True)
    parser.add_argument('--infernal', required=True)
    parser.add_argument('--families', required=True)
    parser.add_argument('--prefix', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--audit', required=True)
    args = parser.parse_args()
    trna = read_gff(args.trnascan)
    namespace(trna, args.prefix + '.trnascan')
    rrna = []
    seen = set()
    for number, path in enumerate(args.barrnap, 1):
        features = read_gff(path)
        namespace(features, f'{args.prefix}.barrnap{number}')
        for feature in features:
            # Exact duplicate calls from the two kingdom searches are redundant.
            key = (feature.seqid, feature.start, feature.end, feature.strand,
                   feature.kind, tuple(feature.attrs.get('Name', [])))
            if key not in seen:
                rrna.append(feature)
                seen.add(key)
    priority = {'tRNA': interval_index(trna, 'tRNA'), 'rRNA': interval_index(rrna, 'rRNA')}
    retained = []
    with open(args.audit, 'w') as handle:
        audit = csv.writer(handle, delimiter='\t')
        audit.writerow(['model', 'seqid', 'start', 'end', 'strand', 'decision'])
        for feature in infernal_features(args.infernal, family_types(args.families)):
            conflict = feature.kind in priority and overlaps(feature, priority[feature.kind])
            decision = ('removed_for_tRNAscan-SE' if feature.kind == 'tRNA'
                        else 'removed_for_Barrnap') if conflict else 'retained'
            audit.writerow([feature.attrs['Name'][0], feature.seqid, feature.start,
                            feature.end, feature.strand, decision])
            if not conflict:
                retained.append(feature)
    namespace(retained, args.prefix + '.infernal')
    with open(args.output, 'w') as handle:
        handle.write('##gff-version 3\n')
        for feature in sorted(trna + rrna + retained, key=sort_key):
            handle.write(feature.line() + '\n')


if __name__ == '__main__':
    main()
