#!/usr/bin/env python3
"""Prefix prediction identifiers consistently across GFF3 references."""
import argparse


def prefix_line(line, prefix):
    if line.startswith('#') or not line.strip():
        return line
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) != 9:
        raise ValueError('Expected nine tab-separated GFF3 columns')
    attributes = []
    for attribute in fields[8].split(';'):
        key, separator, value = attribute.partition('=')
        if separator and key in ('ID', 'Parent', 'Derives_from'):
            value = ','.join(prefix + '.' + item if item.startswith('g') else item
                             for item in value.split(','))
        attributes.append(key + separator + value)
    fields[8] = ';'.join(attributes)
    return '\t'.join(fields) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--prefix', required=True)
    args = parser.parse_args()
    with open(args.input) as source, open(args.output, 'w') as target:
        for line in source:
            target.write(prefix_line(line, args.prefix))


if __name__ == '__main__':
    main()
