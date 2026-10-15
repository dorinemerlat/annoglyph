#!/usr/bin/env python3
"""Run local NCBI table2asn validation; never submit files to an archive."""
import argparse
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--template', type=Path, required=True)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--sha256', required=True)
    parser.add_argument('--fail-on-errors', choices=['true', 'false'], default='true')
    args = parser.parse_args()
    digest = hashlib.sha256(args.binary.read_bytes()).hexdigest()
    if digest.lower() != args.sha256.lower():
        raise ValueError('table2asn binary SHA-256 differs from the configured value')
    submission = args.bundle / 'submission'
    metadata = json.loads((submission / 'metadata.json').read_text())
    source = metadata.get('submission', {})
    for field in ['bioproject', 'biosample', 'locus_tag_prefix']:
        if not source.get(field):
            raise ValueError(f'Official validation requires {field} in submission metadata')
    review = (submission / 'validation.tsv').read_text().splitlines()[1:]
    if review:
        raise ValueError('Resolve submission/validation.tsv issues before official validation')
    if not re.fullmatch(r'PRJ[A-Z]+[0-9]+', source['bioproject']) or not re.fullmatch(r'SAM[A-Z]+[0-9]+', source['biosample']):
        raise ValueError('Invalid BioProject/BioSample accession format')
    output = submission / 'ncbi_validation'
    output.mkdir()
    # Matching stems allow table2asn to associate the FASTA and feature table.
    shutil.copyfile(submission / 'genome.fsa', output / 'genome.fsa')
    shutil.copyfile(submission / 'annotation.tbl', output / 'genome.tbl')
    shutil.copyfile(args.template, output / 'template.sbt')
    binary = output / 'table2asn'
    shutil.copyfile(args.binary, binary)
    binary.chmod(0o755)
    command = [str(binary.resolve()), '-indir', str(output.resolve()), '-t', str((output / 'template.sbt').resolve()),
               '-M', 'n', '-Z', '-V', 'v', '-euk']
    with (output / 'execution.log').open('w') as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
    errors = []
    for path in list(output.glob('*.val')) + list(output.glob('*.dr')):
        for line in path.read_text(errors='replace').splitlines():
            if re.search(r'\b(ERROR|FATAL|REJECT)\b', line, re.I):
                errors.append({'file': path.name, 'message': line})
    sqn = list(output.glob('*.sqn'))
    validated = result.returncode == 0 and bool(sqn) and not errors and bool(list(output.glob('*.val')))
    report = {'validator': 'NCBI table2asn', 'binary_sha256': digest, 'command': command,
              'exit_status': result.returncode, 'errors': errors, 'validated': validated,
              'scope': 'local validation; archive acceptance and discrepancy review are separate',
              'archive_metadata': {'bioproject': source['bioproject'], 'biosample': source['biosample']},
              'archive_metadata_requires_portal_association': True, 'submitted': False}
    (output / 'status.json').write_text(json.dumps(report, indent=2) + '\n')
    metadata['submission_validated'] = validated
    (submission / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    binary.unlink()
    if not validated and args.fail_on_errors == 'true':
        raise SystemExit('NCBI validation did not pass; inspect the retained work directory reports')


if __name__ == '__main__':
    main()
