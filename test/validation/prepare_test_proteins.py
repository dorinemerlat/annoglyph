#!/usr/bin/env python3
"""Select reproducible, reviewed yeast evidence without downloading OrthoDB."""
import argparse
import csv
import hashlib
import re
import textwrap
from pathlib import Path


def records(path):
    header, sequence = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(sequence)
                header, sequence = line[1:].strip(), []
            else:
                sequence.append(line.strip())
    if header is not None:
        yield header, ''.join(sequence)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--count', type=int, default=2000)
    args = parser.parse_args()
    unique = {}
    for header, sequence in records(args.source):
        if not header.startswith('sp|') or not 90 <= len(sequence) <= 5000:
            continue
        if not re.fullmatch('[ACDEFGHIKLMNPQRSTVWYXBZJUO]+', sequence):
            continue
        accession = header.split('|')[1]
        taxid = re.search(r'\bOX=(\d+)', header).group(1)
        checksum = hashlib.sha256(sequence.encode()).hexdigest()
        candidate = (accession, taxid, sequence, header, checksum)
        previous = unique.get(checksum)
        # Retain the test's S. paradoxus evidence when proteins are identical.
        if previous is None or (taxid == '27291', accession) > (previous[1] == '27291', previous[0]):
            unique[checksum] = candidate
    candidates = list(unique.values())
    close = sorted((x for x in candidates if x[1] == '27291'), key=lambda x: x[0])
    remainder = sorted((x for x in candidates if x[1] != '27291'), key=lambda x: hashlib.sha256(('annoglyph-test:' + x[0]).encode()).hexdigest())
    chosen = (close + remainder)[:args.count]
    if len(chosen) != args.count:
        raise ValueError(f'Only {len(chosen)} eligible reviewed, nonredundant proteins')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w') as handle:
        for accession, taxid, sequence, header, checksum in chosen:
            handle.write('>' + accession + '\n' + textwrap.fill(sequence, 60) + '\n')
    with args.output.with_suffix('.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['accession', 'taxid', 'length', 'sha256', 'uniprot_header'])
        for accession, taxid, sequence, header, checksum in chosen:
            writer.writerow([accession, taxid, len(sequence), checksum, header])
    print(f'{len(chosen)} proteins, {sum(len(x[2]) for x in chosen)} amino acids, {args.output.stat().st_size} bytes')
    print('S. paradoxus proteins:', sum(x[1] == '27291' for x in chosen))


if __name__ == '__main__':
    main()
