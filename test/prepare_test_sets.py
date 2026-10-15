#!/usr/bin/env python3
"""Create samplesheets and minimal/full evidence and assembly statistics from installed NCBI test assemblies."""
import csv
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def fasta_records(path):
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


def write_evidence(output, records):
    with output.open('w') as handle, output.with_suffix('.tsv').open('w') as provenance:
        writer = csv.writer(provenance, delimiter='\t')
        writer.writerow(['id', 'assembly', 'original_header', 'length', 'sha256'])
        for checksum, accession, header, sequence in records:
            identifier = 'test_' + checksum[:20]
            handle.write('>' + identifier + '\n')
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset:offset + 60] + '\n')
            writer.writerow([identifier, accession, header, len(sequence), checksum])


def assembly_stats(path):
    lengths, gc, canonical, ambiguous = [], 0, 0, 0
    for _, sequence in fasta_records(path):
        sequence = sequence.upper()
        lengths.append(len(sequence))
        gc += sequence.count('G') + sequence.count('C')
        acgt = sum(sequence.count(base) for base in 'ACGT')
        canonical += acgt
        ambiguous += len(sequence) - acgt
    total, cumulative, n50, l50 = sum(lengths), 0, 0, 0
    for index, length in enumerate(sorted(lengths, reverse=True), 1):
        cumulative += length
        if cumulative >= total / 2:
            n50, l50 = length, index
            break
    return [len(lengths), total, n50, l50, max(lengths, default=0),
            round(100 * gc / canonical, 3) if canonical else 0, ambiguous]


def main():
    manifest = json.loads((ROOT / 'assemblies.json').read_text())
    stats_rows = []
    for group, assemblies in manifest.items():
        directory = ROOT / ('test_' + group)
        rows, candidates = [], {}
        for assembly in assemblies:
            accession = assembly['accession']
            folder = directory / accession
            genomes = list(folder.glob('*_genomic.fna'))
            if len(genomes) != 1 or genomes[0].stat().st_size == 0:
                raise ValueError(f'Missing or ambiguous genome: {folder}')
            rows.append([assembly['name'], assembly['taxid'],
                         str(genomes[0].relative_to(ROOT.parent)), '', '', ''])
            stats_rows.append([group, assembly['name'], accession, assembly['taxid'], *assembly_stats(genomes[0])])
            proteins = folder / 'protein.faa'
            if not proteins.exists():
                print(f'{accession}: no NCBI protein annotation; using homologs from the other assemblies')
                continue
            for header, sequence in fasta_records(proteins):
                sequence = sequence.rstrip('*').upper()
                if not sequence or not re.fullmatch('[ACDEFGHIKLMNPQRSTVWYXBZJUO]+', sequence):
                    continue
                checksum = hashlib.sha256(sequence.encode()).hexdigest()
                candidates.setdefault(checksum, (accession, header, sequence))

        # Round-robin selection keeps each annotated assembly represented.
        buckets = {a['accession']: [] for a in assemblies}
        for checksum, record in sorted(candidates.items()):
            if 90 <= len(record[2]) <= 5000:
                buckets[record[0]].append((checksum, *record))
        selected = []
        for index in range(max(map(len, buckets.values()))):
            for bucket in buckets.values():
                if index < len(bucket):
                    selected.append(bucket[index])
                    if len(selected) == 2000:
                        break
            if len(selected) == 2000:
                break
        if len(selected) != 2000:
            raise ValueError(f'{group}: fewer than 2,000 eligible distinct proteins')

        with (ROOT / f'samplesheet_{group}.csv').open('w') as handle:
            writer = csv.writer(handle)
            writer.writerow(['name', 'taxid', 'fasta', 'mnemonic', 'RNASeq-dir', 'SRA'])
            writer.writerows(rows)
        output = directory / 'protein_evidence_2000.fasta'
        write_evidence(output, selected)
        full = [(checksum, *record) for checksum, record in sorted(candidates.items())]
        write_evidence(directory / 'protein_evidence_full.fasta', full)
        print(f'{group}: {len(rows)} genomes; {len(selected)} minimal / {len(full)} full proteins; {output.stat().st_size:,} bytes')

    columns = ['group', 'species', 'accession', 'taxid', 'sequences', 'total_bp',
               'N50_bp', 'L50', 'longest_bp', 'GC_percent_ACGT', 'ambiguous_bp']
    with (ROOT / 'assembly_stats.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(columns)
        writer.writerows(stats_rows)
    lines = ['# Test assembly statistics', '',
             'Computed directly from the installed genomic FASTAs by `prepare_test_sets.py`.',
             'Assemblies are used in full in both test modes. N50 uses all FASTA records;',
             'records may be contigs or scaffolds. GC excludes non-ACGT bases.', '',
             '[Download the TSV](assembly_stats.tsv)', '',
             '| Group | Species | Accession | Sequences | Size (Mb) | N50 (bp) | L50 | GC (%) | Ambiguous bases |',
             '| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for group, species, accession, taxid, count, total, n50, l50, longest, gc, ambiguous in stats_rows:
        lines.append(f'| {group} | *{species}* | {accession} | {count:,} | {total / 1e6:.3f} | {n50:,} | {l50:,} | {gc:.3f} | {ambiguous:,} |')
    (ROOT / 'assembly_stats.md').write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
