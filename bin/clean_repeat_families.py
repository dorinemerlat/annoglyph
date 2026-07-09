#!/usr/bin/env python3
"""
Simplified cleaner for repeat families FASTA headers.
Usage: clean_repeat_families_simple.py -f <fasta> -c <classification.csv> -m <mnemonic>

This script renames headers to the format:
  {MNEMONIC}_TE{N}#{clean_classification}
and a short version for MCHelper:
  {MNEMONIC}_TE{N}#{clean_classification_short}

It is intentionally minimal: no `source` argument and fewer transformations.
"""

import argparse
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def load_classification(path):
    classification = {}
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            codes = [c.strip() for c in row.get('Code', '').split(';') if c.strip()]
            value = {k: v for k, v in row.items() if k != 'Code'}
            for code in codes:
                classification[code.lower()] = value
    return classification


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-f', '--fasta', required=True, help='Input FASTA file')
    p.add_argument('-c', '--classification', required=True, help='Classification CSV file')
    p.add_argument('-m', '--mnemonic', required=True, help='Mnemonic of species')
    args = p.parse_args()

    classification = load_classification(args.classification)

    base = args.fasta.rsplit('.', 1)[0]
    mchelper_out = f"{base}.cleaned_for_mchelper.fa"
    table_out = f"{base}.tsv"

    out_records_mchelper = []
    table_rows = []

    counter = 1

    for rec in SeqIO.parse(args.fasta, 'fasta'):
        seq_id = rec.id
        parts = seq_id.split('#', 1)
        if len(parts) > 1:
            repeat_classification = parts[1].strip()
        else:
            repeat_classification = ''

        if repeat_classification.lower() == 'unknown' or not repeat_classification:
            clean_classification = 'Unknown'
            clean_classification_short = 'Unknown'
            matched = None
            cls = order = superfamily = family = subfamily = subfamily2 = ''
        else:
            matched = classification.get(repeat_classification.lower()) if repeat_classification else None
            if matched is not None:
                cls = (matched.get('Class') or '').strip()
                order = (matched.get('Order') or '').strip()
                superfamily = (matched.get('SuperFamily') or '').strip()
                family = (matched.get('Family') or '').strip()
                subfamily = (matched.get('SubFamily') or '').strip()
                subfamily2 = (matched.get('Subfamily2') or '').strip()
                clean_classification = '/'.join(filter(None, [cls, order, superfamily, family, subfamily, subfamily2]))
                clean_classification_short = '/'.join(filter(None, [cls, order, superfamily]))
            else:
                clean_classification = 'Unknown'
                clean_classification_short = 'Unknown'
                cls = order = superfamily = family = subfamily = subfamily2 = ''

        new_id = f"{args.mnemonic}_TE{counter}"
        counter += 1

        # short record for MCHelper only (final FASTA)
        if matched and 'CMC' in clean_classification:
            clean_classification_short = 'ClassII/TIR/CMC'
        if matched and 'Mutator' in clean_classification:
            clean_classification_short = 'ClassII/TIR/MULE'

        repeat_header_short = f"{new_id}#{clean_classification_short}"
        mchelper_record = SeqRecord(rec.seq, id=repeat_header_short, description='')
        out_records_mchelper.append(mchelper_record)

        # add row to table: id, rm2_id, rm2_length, rm_code, cls, order, superfamily, family, subfamily, subfamily2
        rm_code = repeat_classification
        table_rows.append({
            'id': new_id,
            'rm2_id': parts[0].strip() if parts and parts[0] else '',
            'rm2_length': len(rec.seq),
            'rm_code': rm_code,
            'rm2_class': cls,
            'rm2_order': order,
            'rm2_superfamily': superfamily,
            'rm2_family': family,
            'rm2_subfamily': subfamily,
            'rm2_subfamily2': subfamily2
        })

    # write only MCHelper FASTA (overwrite)
    with open(mchelper_out, 'w') as out:
        SeqIO.write(out_records_mchelper, out, 'fasta')

    # write table TSV
    fieldnames = ['id', 'rm2_id', 'rm2_length', 'rm_code', 'rm2_class', 'rm2_order', 'rm2_superfamily', 'rm2_family', 'rm2_subfamily', 'rm2_subfamily2']
    with open(table_out, 'w', newline='', encoding='utf-8') as taf:
        writer = csv.DictWriter(taf, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for r in table_rows:
            writer.writerow(r)


if __name__ == '__main__':
    main()
