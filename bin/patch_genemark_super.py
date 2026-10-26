#!/usr/bin/env python3
"""Make GeneMark-ETP training summaries report exon counts to parse_set.pl."""

import argparse
from pathlib import Path


def patch_super(path):
    source = path.read_text()
    original = 'print $OUT ("## ". $gcount ."\\t". $cds_length ."\\t". (scalar @introns) ."\\n");'
    corrected = original.replace('(scalar @introns)', '(scalar @cds)')
    if source.count(corrected) == 1 and original not in source:
        return
    if source.count(original) != 1:
        raise ValueError(f'Unrecognized GeneMark-ETP summary format: {path}')
    path.write_text(source.replace(original, corrected))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('super_script', type=Path)
    args = parser.parse_args()
    patch_super(args.super_script)


if __name__ == '__main__':
    main()
