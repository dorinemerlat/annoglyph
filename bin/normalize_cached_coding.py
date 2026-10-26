#!/usr/bin/env python3
"""Validate and translate cached final coding annotations without annotation tools."""
import argparse
import json
import shutil
from pathlib import Path
from prepare_cache_reuse import normalize_gff, reference_index
from filter_cached_coding import filter_models


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample', default='unknown')
    for name in ('source', 'reference', 'mapping', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    with args.mapping.open('rb') as handle:
        next(handle)
        mapping = dict(line.split() for line in handle if line.strip())
    reference = reference_index(args.reference)
    temporary = Path('coding_filtered.gff3')
    exclusions = filter_models(args.source, temporary, Path('excluded_coding_models.gff3'),
                               Path('coding_exclusions.tsv'), mapping, reference)
    source = temporary
    result = normalize_gff(source, args.output, mapping, reference)
    # The input can be a scratch-staged symlink; publish a self-contained file.
    if args.output.is_symlink():
        original = args.output.resolve()
        args.output.unlink()
        shutil.copyfile(original, args.output)
    Path('coding_import.json').write_text(json.dumps({
        'sample': args.sample, 'source': str(args.source.resolve()),
        'validation': result, 'model_exclusions': exclusions,
        'operation': 'reuse selected coding models, translate sequence IDs, and exclude audited inconsistent models'
    }, indent=2) + '\n')


if __name__ == '__main__':
    main()
