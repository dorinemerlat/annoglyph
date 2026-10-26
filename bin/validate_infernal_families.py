#!/usr/bin/env python3
"""Check imported Infernal tables against Rfam metadata before job submission."""
import argparse
import json
from collections import Counter
from pathlib import Path
from merge_ncrna import family_types


def validate(tables, families):
    known = family_types(families)
    missing, used = Counter(), Counter()
    for table in tables:
        with table.open() as handle:
            for number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.split(maxsplit=27)
                if len(fields) < 27:
                    raise ValueError(f'{table}:{number}: expected cmscan --fmt 2 output')
                used[fields[2]] += 1
                if fields[2] not in known:
                    missing[fields[2]] += 1
    return {'tables': len(tables), 'hits': sum(used.values()), 'families_used': len(used),
            'families_available': len(known), 'missing_families': dict(missing),
            'status': 'NEEDS_AUTOMATIC_RESOLUTION' if missing else 'PASS'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources', type=Path, required=True)
    parser.add_argument('--infernal-dir', type=Path, required=True)
    parser.add_argument('--families', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    samples = json.loads(args.sources.read_text())
    tables = [args.infernal_dir / sample / f'{sample}_infernal.tblout' for sample in samples]
    result = validate(tables, args.families)
    result['metadata'] = str(args.families.resolve())
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
