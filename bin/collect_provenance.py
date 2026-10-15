#!/usr/bin/env python3
"""Collect input hashes, database identities, parameters and tool version evidence."""
import argparse
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def identity(path, full):
    path = Path(path)
    if not path.exists():
        return {'requested_path': str(path), 'available': False}
    resolved = path.resolve()
    entries = []
    files = sorted(p for p in resolved.rglob('*') if p.is_file()) if resolved.is_dir() else [resolved]
    for file in files:
        stat = file.stat()
        entry = {'name': str(file.relative_to(resolved)) if resolved.is_dir() else file.name,
                 'bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns}
        # Small release metadata is always hashed; large database bodies are opt-in.
        if full or (stat.st_size < 1024 * 1024 and re.search(r'(version|release|\.cfg$|\.json$)', file.name, re.I)):
            entry['sha256'] = digest(file)
        entries.append(entry)
    return {'requested_path': str(path), 'resolved_path': str(resolved), 'available': True,
            'hash_mode': 'full_sha256' if full else 'inventory_and_small_release_files', 'files': entries}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot', type=Path, required=True)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--versions', type=Path, required=True)
    parser.add_argument('--database-lock', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    snapshot = json.loads(args.snapshot.read_text())
    args.output.mkdir()
    project = Path(snapshot['workflow']['project_dir'])
    source_files = [file for directory in ['bin', 'lib', 'conf', 'workflows', 'subworkflows', 'modules', 'assets']
                    for file in (project / directory).rglob('*') if file.is_file() and file.suffix in ('.py', '.nf', '.groovy', '.config', '.json', '.yml', '.yaml', '.sh')]
    source_files += [file for file in project.glob('*.config')] + [project / 'main.nf', project / 'modules.json']
    (args.output / 'pipeline_files.sha256.json').write_text(json.dumps(
        {str(file.relative_to(project)): digest(file) for file in source_files}, indent=2) + '\n')
    try:
        status = subprocess.run(['git', '-C', str(project), 'status', '--porcelain'], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30, check=False)
        (args.output / 'workspace_status.txt').write_text(status.stdout)
    except (OSError, subprocess.TimeoutExpired):
        pass
    (args.output / 'parameters.json').write_text(json.dumps(snapshot['params'], indent=2) + '\n')
    (args.output / 'execution.json').write_text(json.dumps(snapshot['workflow'], indent=2) + '\n')
    input_records = []
    for item in json.loads(args.inputs.read_text()):
        path = Path(item['staged'])
        input_records.append({'name': item['name'], 'original_path': item['original'],
                              'bytes': path.stat().st_size, 'sha256': digest(path)})
    (args.output / 'inputs.sha256.json').write_text(json.dumps(input_records, indent=2) + '\n')
    full = bool(snapshot['params'].get('provenance_hash_databases'))
    databases = {name: identity(path, full) for name, path in snapshot['databases'].items() if path}
    for item in snapshot.get('used_database_artifacts', []):
        databases[item['name']] = identity(item['staged'], full or item['name'] == 'protein_evidence')
    if args.database_lock:
        locked = json.loads(args.database_lock.read_text())
        for name, expected in locked.items():
            actual = databases.get(name)
            if actual is None or actual.get('files') != expected.get('files'):
                raise ValueError(f'Database lock mismatch: {name}')
    for item in snapshot.get('used_database_artifacts', []):
        if item['name'] in ('protein_evidence_sources', 'uniprot_release'):
            target = args.output / 'database_metadata'
            target.mkdir(exist_ok=True)
            shutil.copyfile(item['staged'], target / (item['name'] + Path(item['staged']).suffix))
    (args.output / 'databases.lock.json').write_text(json.dumps(databases, indent=2) + '\n')
    versions_dir = args.output / 'software_versions'
    versions_dir.mkdir()
    for index, item in enumerate(json.loads(args.versions.read_text())):
        shutil.copyfile(item, versions_dir / f'nfcore_{index:04d}.yml')
    # Local modules lack standardized versions.yml; collect --version evidence
    # from cached runtime containers without running biological analyses.
    commands = {'braker3': ['braker.pl', '--version'], 'galba': ['galba.pl', '--version'],
                'varus': ['STAR', '--version'], 'trnascan': ['tRNAscan-SE', '--version'],
                'infernal': ['cmscan', '-h'], 'omark': ['omark', '--version'],
                'blast': ['blastp', '-version'], 'exogap_tools': ['barrnap', '--version']}
    cache = Path(snapshot['singularity_cache'])
    images = []
    if cache.exists():
        for image in sorted(cache.glob('*')):
            if not image.is_file() or image.suffix not in ('.img', '.sif'):
                continue
            record = {'path': str(image.resolve()), 'bytes': image.stat().st_size,
                      'mtime_ns': image.stat().st_mtime_ns}
            # Inventory, not proof that every cached image was used in this run.
            images.append(record)
        for tool, command in commands.items():
            matches = [p for p in cache.glob('*') if tool.replace('_', '-') in p.name.replace('_', '-') and p.suffix in ('.img', '.sif')]
            for image in matches:
                try:
                    result = subprocess.run(['singularity', 'exec', str(image), *command], text=True,
                                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60, check=False)
                    output = {'image': str(image), 'command': command, 'exit': result.returncode,
                              'output': result.stdout[:16000], 'scope': 'cached-image version probe; not execution evidence'}
                except (OSError, subprocess.TimeoutExpired) as error:
                    output = {'image': str(image), 'command': command, 'unavailable': str(error)}
                (versions_dir / (tool + '_' + hashlib.sha256(str(image).encode()).hexdigest()[:8] + '.json')).write_text(json.dumps(output, indent=2) + '\n')
    (args.output / 'cached_images.inventory.json').write_text(json.dumps(images, indent=2) + '\n')
    (args.output / 'limitations.json').write_text(json.dumps({
        'database_hash_mode': 'full' if full else 'metadata_only_for_large_database_files',
        'cached_images_are_not_a_list_of_executed_images': True,
        'nfcore_versions_are_task_output_evidence': True,
        'local_version_probes_do_not_establish_which_cached_image_was_used': True,
        'trace_and_runtime_container_paths_must_be_used_to_identify_executed_images': True,
        'database_lock_does_not_invalidate_analysis_caches': True}, indent=2) + '\n')


if __name__ == '__main__':
    main()
