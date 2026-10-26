#!/usr/bin/env python3
"""Complete cached bundles with the exact pooled RepeatMasker input library."""
import hashlib
import json
import os
import shutil
from pathlib import Path


def main():
    project = Path(__file__).resolve().parents[1]
    manifest = json.loads((project / 'cache_save/cache_reuse_v1/manifest.json').read_text())
    sources = json.loads((project / 'outputs/cache_reuse_20261023/final_coding_sources.json').read_text())
    library = Path(manifest['repeat_library'])
    digest = hashlib.sha256()
    with library.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    expected = digest.hexdigest()
    relative = f'repeats/sequences/{library.name}'
    for genome in sources:
        for bundle in [project / 'cache_save/finalize_cache_save_v3/export_annotation' / genome / genome,
                       project / 'results/species' / genome]:
            checksums = bundle / 'submission/checksums.sha256'
            if not checksums.is_file():
                raise ValueError(f'Bundle not complete: {bundle}')
            entries = {name: value for value, name in
                       (line.split('  ', 1) for line in checksums.read_text().splitlines())}
            destination = bundle / relative
            if entries.get(relative) != expected or not destination.is_file() or destination.stat().st_size != library.stat().st_size:
                destination.parent.mkdir(parents=True, exist_ok=True)
                temporary = destination.with_name(destination.name + '.tmp')
                shutil.copyfile(library, temporary)
                os.replace(temporary, destination)
                entries[relative] = expected
                temporary = checksums.with_name(checksums.name + '.tmp')
                temporary.write_text(''.join(f'{value}  {name}\n' for name, value in sorted(entries.items())))
                os.replace(temporary, checksums)
        print(f'{genome}: pooled repeat library included', flush=True)
    report = project / 'results/pipeline_info/finalization/repeat_library.json'
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps({'source': str(library), 'sha256': expected,
                                  'bytes': library.stat().st_size, 'genomes': len(sources)}, indent=2) + '\n')


if __name__ == '__main__':
    main()
