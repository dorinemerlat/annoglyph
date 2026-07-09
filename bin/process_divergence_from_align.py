#!/usr/bin/env python3
"""
process_divergence_from_align.py

Convert a divsum section into % of assembly size per column.

Usage: process_divergence_from_align.py -d <divsum_file> -o <output_file> -s <assembly_size>

This replicates the behaviour of the provided bash script: find the line
"Coverage for each repeat class and divergence (Kimura)", take the next line as
header and convert subsequent numeric columns by (value / assembly_size * 100).
"""

import argparse
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Process divergence section and convert counts to percent of assembly size"
    )
    p.add_argument('-d', '--divsum', required=True, help='Input .divsum file')
    p.add_argument('-o', '--output', required=True, help='Output TSV file')
    p.add_argument('-s', '--assembly-size', required=True, type=float, help='Assembly size (number)')
    return p.parse_args()


def main():
    args = parse_args()

    divsum_path = Path(args.divsum)
    if not divsum_path.exists():
        print(f"Error: divsum file not found: {divsum_path}", file=sys.stderr)
        sys.exit(2)

    try:
        assembly_size = float(args.assembly_size)
    except Exception:
        print("Error: assembly size must be a number", file=sys.stderr)
        sys.exit(2)

    if assembly_size == 0:
        print("Error: assembly size cannot be zero", file=sys.stderr)
        sys.exit(2)

    start_marker = 'Coverage for each repeat class and divergence (Kimura)'

    with divsum_path.open('r', encoding='utf-8') as fh:
        lines = fh.readlines()

    # Find the line index of the marker
    begin_idx = None
    for i, line in enumerate(lines):
        if start_marker in line:
            begin_idx = i
            break

    if begin_idx is None:
        print(f"Error: start marker not found in {divsum_path}", file=sys.stderr)
        sys.exit(3)

    # The bash pipeline took lines from marker_line..end, then removed the first line
    # so the header is the line immediately AFTER the marker
    header_idx = begin_idx + 1
    if header_idx >= len(lines):
        print("Error: unexpected file format (no header line after marker)", file=sys.stderr)
        sys.exit(4)

    header_line = lines[header_idx].rstrip('\n')
    data_lines = lines[header_idx + 1 :]

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open('w', encoding='utf-8') as outfh:
        # write header exactly as in the file
        outfh.write(header_line + '\n')

        for raw in data_lines:
            line = raw.strip()
            if not line:
                continue

            parts = line.split()
            # keep first column as-is (class name)
            first = parts[0]
            rest = parts[1:]
            new_vals = []
            for v in rest:
                try:
                    num = float(v)
                except Exception:
                    # If not a number, keep original token
                    new_vals.append(v)
                    continue

                # convert to percent of assembly size
                converted = num / assembly_size * 100.0
                # format with sensible precision (6 decimal places, strip trailing zeros)
                s = ('{:.6f}'.format(converted)).rstrip('0').rstrip('.')
                new_vals.append(s)

            outfh.write(first + '\t' + '\t'.join(new_vals) + '\n')

    print(f"Wrote: {out_path}")


if __name__ == '__main__':
    main()
