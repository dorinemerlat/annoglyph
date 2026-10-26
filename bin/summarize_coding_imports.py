#!/usr/bin/env python3
"""Summarize the audited exclusions from cached selected coding annotations."""
import argparse
import csv
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audits', nargs='+', type=Path, required=True)
    args = parser.parse_args()
    with open('coding_import_summary.tsv', 'w') as summary, open('excluded_models.tsv', 'w') as exclusions:
        summaries, details = csv.writer(summary, delimiter='\t'), csv.writer(exclusions, delimiter='\t')
        summaries.writerow(['genome_id', 'excluded_transcripts', 'excluded_feature_rows', 'retained_feature_rows', 'source_gff'])
        details.writerow(['genome_id', 'model_id', 'reason', 'source_gff'])
        records = [json.loads(path.read_text()) for path in args.audits]
        for record in sorted(records, key=lambda r: r['sample']):
            counts = record['model_exclusions']
            summaries.writerow([record['sample'], counts['excluded_transcripts'], counts['excluded_feature_rows'],
                                counts['retained_feature_rows'], record['source']])
            for identifier, reasons in sorted(counts['reasons'].items()):
                details.writerow([record['sample'], identifier, ';'.join(reasons), record['source']])


if __name__ == '__main__':
    main()
