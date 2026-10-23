#!/usr/bin/env python3
"""Validate and translate cached final coding annotations without annotation tools."""
import argparse
import json
import shutil
from pathlib import Path
from prepare_cache_reuse import normalize_gff, reference_index


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'reference', 'mapping', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    with args.mapping.open('rb') as handle:
        next(handle)
        mapping = dict(line.split() for line in handle if line.strip())
    reference = reference_index(args.reference)
    result = normalize_gff(args.source, args.output, mapping, reference)
    # The input can be a scratch-staged symlink; publish a self-contained file.
    if args.output.is_symlink():
        original = args.output.resolve()
        args.output.unlink()
        shutil.copyfile(original, args.output)
    Path('coding_import.json').write_text(json.dumps({
        'source': str(args.source.resolve()), 'validation': result,
        'operation': 'reuse selected coding models; translate sequence IDs only'
    }, indent=2) + '\n')


if __name__ == '__main__':
    main()
