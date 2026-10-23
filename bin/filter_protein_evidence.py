#!/usr/bin/env python3
"""Exclude the target genus from automatically downloaded protein evidence."""
import argparse
import csv
import gzip
import json
import re
import sqlite3
from pathlib import Path


def records(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        header, lines = None, []
        for line in handle:
            if line.startswith('>'):
                if header is not None:
                    yield header, lines
                header, lines = line[1:].strip(), []
            elif header is not None:
                lines.append(line)
        if header is not None:
            yield header, lines


class Taxonomy:
    def __init__(self, path):
        self.db = sqlite3.connect(f'file:{Path(path).resolve()}?mode=ro', uri=True)
        self.cache = {}

    def genus(self, taxid):
        taxid = int(taxid)
        if taxid not in self.cache:
            row = self.db.execute('SELECT track FROM species WHERE taxid=?', (taxid,)).fetchone()
            genus = None
            if row:
                for ancestor in row[0].split(','):
                    node = self.db.execute('SELECT spname, rank FROM species WHERE taxid=?', (int(ancestor),)).fetchone()
                    if node and node[1] == 'genus':
                        genus = node[0]
                        break
            self.cache[taxid] = genus
        return self.cache[taxid]


def protein_taxid(header):
    match = re.search(r'\bOX=(\d+)\b', header)
    if match:
        return int(match[1])
    # Partitioned OrthoDB identifiers: NCBI taxid followed by organism/gene ID.
    match = re.match(r'(\d+)[_:]', header)
    return int(match[1]) if match else None


def filter_evidence(proteins, taxonomy, target_genus, output, report):
    counts = dict(total=0, retained=0, excluded_target_genus=0, excluded_unresolved_taxonomy=0)
    with open(output, 'w') as handle:
        for header, lines in records(proteins):
            counts['total'] += 1
            taxid = protein_taxid(header)
            genus = taxonomy.genus(taxid) if taxid else None
            if genus is None:
                counts['excluded_unresolved_taxonomy'] += 1
            elif genus.casefold() == target_genus.casefold():
                counts['excluded_target_genus'] += 1
            else:
                counts['retained'] += 1
                handle.write('>' + header + '\n')
                handle.writelines(lines)
    Path(report).write_text(json.dumps(dict(excluded_genus=target_genus, **counts), indent=2) + '\n')
    if not counts['retained']:
        raise ValueError(f'No usable protein evidence remains after excluding genus {target_genus}; see {report}')
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proteins', required=True, type=Path)
    parser.add_argument('--taxonomy', required=True, type=Path)
    parser.add_argument('--target-id', required=True)
    parser.add_argument('--taxonomy-db', required=True, type=Path)
    parser.add_argument('--output', default='filtered_protein_evidence.fa')
    parser.add_argument('--report', default='protein_evidence_filter.json')
    args = parser.parse_args()
    with args.taxonomy.open() as handle:
        rows = [row for row in csv.DictReader(handle, delimiter='\t') if row['id'] == args.target_id]
    if len(rows) != 1 or not rows[0].get('genus', '').strip():
        raise ValueError(f'Expected one taxonomy row with a genus for {args.target_id}')
    filter_evidence(args.proteins, Taxonomy(args.taxonomy_db), rows[0]['genus'].strip(), args.output, args.report)


if __name__ == '__main__':
    main()
