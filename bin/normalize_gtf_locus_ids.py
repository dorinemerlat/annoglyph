#!/usr/bin/env python3
"""Rebuild GTF parent rows with identifiers unique to each genomic locus."""
import argparse
from collections import defaultdict
from pathlib import Path
import re


GENE_ID = re.compile(r'\bgene_id "([^"]+)"')
TRANSCRIPT_ID = re.compile(r'\btranscript_id "([^"]+)"')


def interval_groups(records, adjacent=False):
    groups = []
    maximum = None
    for record in sorted(records, key=lambda row: (row['start'], row['end'])):
        if maximum is None or record['start'] > maximum + int(adjacent):
            groups.append([])
            maximum = record['end']
        groups[-1].append(record)
        maximum = max(maximum, record['end'])
    return groups


def parent_row(children, kind, gene, transcript=None):
    fields = children[0]['fields'][:]
    fields[2:6] = [kind, str(min(row['start'] for row in children)),
                   str(max(row['end'] for row in children)), '.']
    fields[7] = '.'
    fields[8] = f'gene_id "{gene}";'
    if transcript:
        fields[8] = f'transcript_id "{transcript}"; ' + fields[8]
    return '\t'.join(fields) + '\n'


def normalize_lines(lines):
    comments = []
    children = defaultdict(list)
    gene_bounds = defaultdict(list)
    transcript_bounds = defaultdict(list)
    for line in lines:
        if line.startswith('#') or not line.strip():
            comments.append(line)
            continue
        fields = line.rstrip('\r\n').split('\t')
        if len(fields) != 9:
            raise ValueError('Expected nine GTF columns')
        start, end = int(fields[3]), int(fields[4])
        if start < 1 or end < start:
            raise ValueError('Invalid GTF coordinates')
        gene_match = GENE_ID.search(fields[8])
        transcript_match = TRANSCRIPT_ID.search(fields[8])
        if fields[2] == 'gene':
            gene = gene_match.group(1) if gene_match else fields[8]
            gene_bounds[(gene, fields[0], fields[6])].append((start, end))
            continue
        if fields[2] in ('transcript', 'mRNA'):
            transcript = transcript_match.group(1) if transcript_match else fields[8]
            transcript_bounds[(fields[0], fields[6], transcript)].append((start, end))
            continue
        if not gene_match or not transcript_match:
            raise ValueError('Expected gene_id and transcript_id attributes')
        key = (gene_match.group(1), fields[0], fields[6], transcript_match.group(1))
        children[key].append({'start': start, 'end': end, 'fields': fields})
    if not children:
        raise ValueError('No transcript features in GTF input')

    gene_transcripts = defaultdict(list)
    for (gene, sequence, strand, transcript), records in children.items():
        bounds = transcript_bounds[(sequence, strand, transcript)]
        covered = any(start <= min(row['start'] for row in records)
                      and max(row['end'] for row in records) <= end for start, end in bounds)
        parts = [records] if covered else interval_groups(records, adjacent=True)
        for part in parts:
            gene_transcripts[(gene, sequence, strand)].append({
                'start': min(row['start'] for row in part),
                'end': max(row['end'] for row in part),
                'transcript': transcript, 'records': part,
            })
    gene_loci = defaultdict(list)
    for (gene, sequence, strand), transcripts in gene_transcripts.items():
        covered = any(start <= min(row['start'] for row in transcripts)
                      and max(row['end'] for row in transcripts) <= end
                      for start, end in gene_bounds[(gene, sequence, strand)])
        groups = [transcripts] if covered else interval_groups(transcripts)
        for group in groups:
            gene_loci[gene].append((sequence, strand, group))

    output = comments[:]
    used_transcripts = set()
    for gene, loci in gene_loci.items():
        for index, (sequence, strand, transcripts) in enumerate(sorted(
                loci, key=lambda locus: (locus[0], locus[1], min(row['start'] for row in locus[2]))), 1):
            new_gene = gene if len(loci) == 1 else f'{gene}_locus{index}'
            records = [row for transcript in transcripts for row in transcript['records']]
            output.append(parent_row(records, 'gene', new_gene))
            for transcript in transcripts:
                old_transcript = transcript['transcript']
                identifier = old_transcript
                if new_gene != gene:
                    identifier = (new_gene + old_transcript[len(gene):]
                                  if old_transcript.startswith(gene + '.')
                                  else old_transcript + '_' + new_gene)
                candidate = identifier
                suffix = 1
                while identifier in used_transcripts:
                    identifier = f'{candidate}_part{suffix}'
                    suffix += 1
                used_transcripts.add(identifier)
                output.append(parent_row(transcript['records'], 'transcript', new_gene, identifier))
                for record in transcript['records']:
                    fields = record['fields'][:]
                    fields[8] = GENE_ID.sub(lambda match: f'gene_id "{new_gene}"', fields[8])
                    fields[8] = TRANSCRIPT_ID.sub(lambda match: f'transcript_id "{identifier}"', fields[8])
                    output.append('\t'.join(fields) + '\n')
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(''.join(normalize_lines(args.input.read_text().splitlines(keepends=True))))


if __name__ == '__main__':
    main()
