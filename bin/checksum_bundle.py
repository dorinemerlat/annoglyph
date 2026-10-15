#!/usr/bin/env python3
import argparse
import hashlib
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--bundle', type=Path, required=True)
args = parser.parse_args()
destination = args.bundle / 'submission/checksums.sha256'
with destination.open('w') as handle:
    for path in sorted(args.bundle.rglob('*')):
        if not path.is_file() or path == destination:
            continue
        digest = hashlib.sha256()
        with path.open('rb') as source:
            for block in iter(lambda: source.read(1024 * 1024), b''):
                digest.update(block)
        handle.write(f'{digest.hexdigest()}  {path.relative_to(args.bundle)}\n')
