#!/usr/bin/env python3
"""Resolve new Rfam types automatically; keep unresolved hits unclassified."""
import argparse
import fcntl
import gzip
import hashlib
import json
import os
import sys
import tempfile
import time
from pathlib import Path
from urllib.request import urlopen
from rfam_metadata import family_types, infernal_accessions

CURRENT_URL = 'https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz'


def current_metadata(cache):
    cache.mkdir(parents=True, exist_ok=True)
    target = cache / 'family.txt'
    with (cache / '.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if target.exists() and time.time() - target.stat().st_mtime < 86400:
            return family_types(target), target, 'cached Rfam CURRENT metadata'
        try:
            with urlopen(CURRENT_URL, timeout=30) as response:
                data = gzip.decompress(response.read())
            with tempfile.NamedTemporaryFile(dir=cache, delete=False) as handle:
                temporary = Path(handle.name)
                handle.write(data)
            try:
                types = family_types(temporary)
                if not types:
                    raise ValueError('Empty Rfam metadata download')
                os.replace(temporary, target)
            finally:
                temporary.unlink(missing_ok=True)
            return types, target, 'downloaded Rfam CURRENT metadata'
        except (OSError, ValueError, EOFError) as error:
            print(f'WARNING: Rfam metadata update unavailable: {error}. Unresolved hits will be kept as Unclassified ncRNA.', file=sys.stderr)
            if target.exists():
                return family_types(target), target, 'stale cached Rfam CURRENT metadata; update unavailable'
            return {}, None, 'Rfam CURRENT metadata unavailable'


def resolve(families, tables, cache, output):
    types = family_types(families)
    required = infernal_accessions(tables)
    missing = required - types.keys()
    report = {'input_metadata': str(families.resolve()),
              'input_sha256': hashlib.sha256(families.read_bytes()).hexdigest(),
              'required_families': len(required), 'absent_from_input': sorted(missing),
              'automatically_resolved': [], 'unclassified': []}
    if missing:
        supplemental, source, status = current_metadata(cache)
        found = missing & supplemental.keys()
        types.update({accession: supplemental[accession] for accession in found})
        report.update(automatically_resolved=sorted(found), lookup=status)
        if source:
            report.update(supplemental_url=CURRENT_URL,
                          supplemental_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    report['unclassified'] = sorted(required - types.keys())
    output.mkdir(parents=True, exist_ok=True)
    # Only the official accession and type columns are needed downstream.
    with (output / 'family.txt').open('w') as handle:
        for accession, kind in sorted(types.items()):
            handle.write('\t'.join([accession] + [''] * 17 + [kind]) + '\n')
    (output / 'metadata_resolution.json').write_text(json.dumps(report, indent=2) + '\n')
    if report['unclassified']:
        print('WARNING: Keeping unresolved families as Unclassified ncRNA: ' + ', '.join(report['unclassified']), file=sys.stderr)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--families', type=Path, required=True)
    parser.add_argument('--infernal', type=Path, nargs='+', required=True)
    parser.add_argument('--cache-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(resolve(args.families, args.infernal, args.cache_dir, args.output)))


if __name__ == '__main__':
    main()
