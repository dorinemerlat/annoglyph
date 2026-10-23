#!/usr/bin/env python3
"""Keep the first N SRA spots, checking that paired FASTQ records stay synchronized."""
import argparse
from contextlib import ExitStack
import gzip
from pathlib import Path
import re


def read_record(handle):
    header = handle.readline()
    if not header:
        return None
    sequence, plus, quality = handle.readline(), handle.readline(), handle.readline()
    if not header.startswith('@') or not plus.startswith('+') or not quality:
        raise ValueError('Incomplete or invalid FASTQ record')
    if len(sequence.rstrip('\r\n')) != len(quality.rstrip('\r\n')):
        raise ValueError('FASTQ sequence and quality lengths differ')
    return header, sequence, plus, quality


def spot_id(record):
    return re.sub(r'/[12]$', '', record[0].split()[0])


def subsample(inputs, output, limit):
    if limit < 1 or len(inputs) not in (1, 2):
        raise ValueError('Provide a positive spot limit and one or two FASTQ files')
    output.mkdir(parents=True, exist_ok=True)
    if len({path.name for path in inputs}) != len(inputs):
        raise ValueError('FASTQ filenames must be distinct')
    retained = 0
    with ExitStack() as stack:
        readers = [stack.enter_context(gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) for path in inputs]
        writers = [stack.enter_context(gzip.open(output / (path.name if path.suffix == '.gz' else path.name + '.gz'), 'wt')) for path in inputs]
        for _ in range(limit):
            records = [read_record(reader) for reader in readers]
            if all(record is None for record in records):
                break
            if any(record is None for record in records):
                raise ValueError('Paired FASTQ files contain different numbers of spots')
            if len(records) == 2 and spot_id(records[0]) != spot_id(records[1]):
                raise ValueError('Paired FASTQ records do not describe the same spot')
            for writer, record in zip(writers, records):
                writer.writelines(record)
            retained += 1
    return retained


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reads', nargs='+', type=Path, required=True)
    parser.add_argument('--max-spots', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    count = subsample(args.reads, args.output, args.max_spots)
    print(f'Retained {count} spots ({len(args.reads)} FASTQ file(s))')


if __name__ == '__main__':
    main()
