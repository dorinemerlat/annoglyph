#!/usr/bin/env python3
"""Download Rfam models and metadata; cmpress runs in the Infernal container."""
import argparse
import gzip
import re
import shutil
from pathlib import Path
from urllib.request import urlopen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not re.fullmatch(r'CURRENT|\d+\.\d+', args.release):
        parser.error('Rfam release must be CURRENT or a numeric version')
    base = f'https://ftp.ebi.ac.uk/pub/databases/Rfam/{args.release}/'
    args.output.mkdir()
    for remote, local in [('Rfam.cm.gz', 'Rfam.cm'), ('Rfam.clanin', 'Rfam.clanin'),
                          ('database_files/family.txt.gz', 'family.txt')]:
        with urlopen(base + remote, timeout=120) as response:
            stream = gzip.GzipFile(fileobj=response) if remote.endswith('.gz') else response
            with (args.output / local).open('wb') as output:
                shutil.copyfileobj(stream, output)
    (args.output / 'release.txt').write_text(args.release + '\n')


if __name__ == '__main__':
    main()
