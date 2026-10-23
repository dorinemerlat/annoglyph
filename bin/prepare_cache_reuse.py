#!/usr/bin/env python3
"""Build a validated, non-destructive import of an older AnnoGlyph cache.

Original predictions and assemblies are never edited. Only query sequence names
are translated; derived gene filtering/QC/export outputs are not imported.
"""
import argparse
import csv
import hashlib
import json
import os
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

TOKEN = re.compile(rb'(?<![A-Za-z0-9_])[A-Za-z][A-Za-z0-9]*_[0-9]+(?![A-Za-z0-9_])')


def reference_index(path):
    records = {}
    seqid = None
    digest = None
    length = 0
    with path.open('rb') as handle:
        for line in handle:
            if line.startswith(b'>'):
                if seqid is not None:
                    records[seqid] = (length, digest.hexdigest())
                seqid = line[1:].split()[0]
                if seqid in records:
                    raise ValueError(f'{path}: duplicate sequence {seqid!r}')
                digest = hashlib.sha256()
                length = 0
            elif line.strip():
                if digest is None:
                    raise ValueError(f'{path}: sequence before FASTA header')
                sequence = line.strip().upper()
                digest.update(sequence)
                length += len(sequence)
    if seqid is not None:
        records[seqid] = (length, digest.hexdigest())
    if not records:
        raise ValueError(f'{path}: empty FASTA')
    return records


def link(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.is_symlink() and target.resolve() == source.resolve():
        return
    if target.exists() or target.is_symlink():
        raise FileExistsError(target)
    target.symlink_to(os.path.relpath(source, target.parent))


def normalize_masked(source, target, mapping, reference):
    """Verify every base against the cached assembly, ignoring softmask case."""
    changed = any(k != v for k, v in mapping.items())
    temporary = target.with_name(target.name + '.tmp')
    target.parent.mkdir(parents=True, exist_ok=True)
    output = temporary.open('wb') if changed else None
    seen = set()
    seqid = None
    digest = None
    length = 0

    def finish():
        if seqid is not None:
            if seqid in seen or reference.get(seqid) != (length, digest.hexdigest()):
                raise ValueError(f'{source}: masked sequence differs from assembly: {seqid!r}')
            seen.add(seqid)

    try:
        with source.open('rb') as handle:
            for line in handle:
                if line.startswith(b'>'):
                    finish()
                    old = line[1:].split()[0]
                    seqid = mapping.get(old, old)
                    if seqid not in reference:
                        raise ValueError(f'{source}: unknown sequence {old!r}')
                    digest = hashlib.sha256()
                    length = 0
                    if output:
                        output.write(b'>' + seqid + line[1 + len(old):])
                else:
                    if line.strip():
                        if digest is None:
                            raise ValueError(f'{source}: missing FASTA header')
                        sequence = line.strip().upper()
                        digest.update(sequence)
                        length += len(sequence)
                    if output:
                        output.write(line)
        finish()
        if seen != reference.keys():
            raise ValueError(f'{source}: incomplete masked assembly ({len(seen)}/{len(reference)})')
        if output:
            output.close()
            os.replace(temporary, target)
        else:
            link(source, target)
    except BaseException:
        if output:
            output.close()
        temporary.unlink(missing_ok=True)
        raise
    return {'sequences': len(seen), 'renamed': changed, 'bases_verified': sum(v[0] for v in reference.values())}


def normalize_gff(source, target, mapping, reference, zero_length_audit=None):
    """Translate only sequence references, preserving feature IDs and attributes."""
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + '.tmp')
    rows = renamed = excluded = 0
    changed_directives = False
    audit = zero_length_audit.open('wb') if zero_length_audit else None
    try:
        with source.open('rb') as incoming, temporary.open('wb') as outgoing:
            for number, line in enumerate(incoming, 1):
                if line.startswith(b'##FASTA'):
                    raise ValueError(f'{source}:{number}: embedded FASTA requires explicit import')
                if line.startswith(b'##sequence-region '):
                    fields = line.split()
                    previous = fields[1]
                    fields[1] = mapping.get(previous, previous)
                    changed_directives |= fields[1] != previous
                    if fields[1] not in reference:
                        raise ValueError(f'{source}:{number}: unknown sequence-region')
                    line = b' '.join(fields) + b'\n'
                elif line.strip() and not line.startswith(b'#'):
                    fields = line.rstrip(b'\r\n').split(b'\t', 8)
                    if len(fields) != 9:
                        raise ValueError(f'{source}:{number}: expected 9 GFF columns')
                    old = fields[0]
                    fields[0] = mapping.get(old, old)
                    if fields[0] not in reference:
                        raise ValueError(f'{source}:{number}: unknown sequence {old!r}')
                    start, end = int(fields[3]), int(fields[4])
                    if zero_length_audit and start == end + 1 and 0 <= end <= reference[fields[0]][0]:
                        audit.write(line)
                        excluded += 1
                        continue
                    if not 1 <= start <= end <= reference[fields[0]][0]:
                        raise ValueError(f'{source}:{number}: coordinates outside cached assembly')
                    renamed += fields[0] != old
                    rows += 1
                    line = b'\t'.join(fields) + b'\n'
                outgoing.write(line)
        if audit:
            audit.close()
        if renamed or excluded or changed_directives:
            os.replace(temporary, target)
        else:
            temporary.unlink()
            link(source, target)
    except BaseException:
        if audit:
            audit.close()
        temporary.unlink(missing_ok=True)
        raise
    result = {'rows': rows, 'renamed_rows': renamed}
    if zero_length_audit:
        result.update(excluded_zero_length_rows=excluded, exclusion_audit=str(zero_length_audit))
    return result


def normalize_text(source, target, mapping):
    """Update contig references in tRNAscan sidecars in a single simultaneous pass."""
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + '.tmp')
    changed = False
    with source.open('rb') as incoming, temporary.open('wb') as outgoing:
        for line in incoming:
            converted = TOKEN.sub(lambda match: mapping.get(match[0], match[0]), line)
            changed |= converted != line
            outgoing.write(converted)
    if changed:
        os.replace(temporary, target)
    else:
        temporary.unlink()
        link(source, target)


def identity(path):
    stat = path.stat()
    return {'path': str(path.resolve()), 'bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns}


def prepare_genome(cache, output, genome):
    reference_path = cache / f'preprocessing/reformat_genome/{genome}/genome_{genome}.fa'
    mapping_path = reference_path.with_name(reference_path.name + '.tsv')
    with mapping_path.open('rb') as handle:
        next(handle)
        mapping = dict(line.split() for line in handle if line.strip())
    if len(set(mapping.values())) != len(mapping):
        raise ValueError(f'{mapping_path}: mapping is not one-to-one')
    print(f'{genome}: checking reference sequences', flush=True)
    reference = reference_index(reference_path)
    if reference.keys() != set(mapping.values()):
        raise ValueError(f'{mapping_path}: mapping does not describe cached assembly')
    prefixes = {key.rsplit(b'_', 1)[0].decode() for key in reference}
    if len(prefixes) != 1:
        raise ValueError(f'{reference_path}: multiple genome mnemonics')
    target = output / genome
    entry = {'reference': str(reference_path), 'mnemonic': prefixes.pop(),
             'mapping': identity(mapping_path), 'reference_identity': identity(reference_path),
             'predictions': {}, 'prediction_status': {}, 'validation': {}}
    masked_source = cache / f'repeats_annotation/repeatmasker/{genome}/genome_{genome}.fa.masked'
    masked = target / masked_source.name
    entry['validation']['masked'] = normalize_masked(masked_source, masked, mapping, reference)
    print(f'{genome}: masked assembly verified', flush=True)
    entry['masked'] = str(masked)
    entry['masked_source'] = identity(masked_source)
    repeat_source = cache / f'repeats_annotation/repeatmasker_out_to_gff/{genome}/{genome}_repeats.gff'
    repeat_gff = target / repeat_source.name
    entry['validation']['repeats'] = normalize_gff(repeat_source, repeat_gff, mapping, reference, target / 'excluded_zero_length_repeats.gff')
    entry['repeats_gff'] = str(repeat_gff)
    entry['repeats_source'] = identity(repeat_source)
    for key, relative in {
        'family_table': f'repeats_annotation/postprocess_mchelper/{genome}/{genome}_families_processed.tsv',
        'families': f'repeats_annotation/postprocess_mchelper/{genome}/{genome}_families_processed.fa',
        'landscape': f'repeats_annotation/reformat_repeatlandscape/{genome}/{genome}_TE_repeat_landscape_reformated.tsv',
    }.items():
        source = cache / relative
        if not source.is_file():
            raise FileNotFoundError(source)
        entry[key] = str(source)
    for method, relative in {
        'braker2': f'genes_annotation/structural_annotation/braker2/{genome}/braker_{genome}/braker.gff3',
        'galba': f'genes_annotation/structural_annotation/galba/{genome}/galba_{genome}/galba.gff3',
        'braker-with-star': f'genes_annotation/structural_annotation/braker3/star/{genome}/braker_{genome}/braker.gff3',
        'braker-with-varus': f'genes_annotation/structural_annotation/braker3/varus/{genome}/braker_{genome}/braker.gff3',
    }.items():
        source = cache / relative
        if not source.is_file() or not source.stat().st_size:
            entry['prediction_status'][method] = 'empty' if source.is_file() else 'absent'
            continue
        dest = target / 'predictions' / method / source.name
        try:
            validation = normalize_gff(source, dest, mapping, reference)
            if not validation['rows']:
                raise ValueError('GFF contains no features')
        except (ValueError, UnicodeError) as error:
            entry['prediction_status'][method] = 'invalid: ' + str(error)
            continue
        entry['predictions'][method] = str(dest)
        entry['prediction_status'][method] = 'reused'
        entry['validation'][method] = validation | {'source': identity(source)}
    if not entry['predictions']:
        raise ValueError(f'{genome}: no usable gene prediction')
    trna_source = cache / f'ncgenes_annotation/trnascanse/{genome}'
    trna_target = target / 'trnascanse'
    trna_paths = [trna_source / f'{genome}_trna.{ext}' for ext in ('gff', 'fa', 'out', 'struct', 'stats')]
    if all(p.is_file() for p in trna_paths):
        for source in trna_paths:
            dest = trna_target / source.name
            if source.suffix == '.gff':
                entry['validation']['trna'] = normalize_gff(source, dest, mapping, reference)
            else:
                normalize_text(source, dest, mapping)
        entry['trnascanse'] = str(trna_target)
    else:
        entry['trnascanse'] = None
    entry['barrnap'] = {}
    for kingdom in ('euk', 'mito'):
        source = cache / f'ncgenes_annotation/barrnap/{genome}/{genome}_barrnap_{kingdom}.gff'
        if source.is_file():
            dest = target / source.name
            entry['validation']['barrnap_' + kingdom] = normalize_gff(source, dest, mapping, reference)
            entry['barrnap'][kingdom] = str(dest)
    infernal = cache / f'ncgenes_annotation/infernal/{genome}/{genome}_infernal.tblout'
    count = 0
    with infernal.open('rb') as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(b'#'):
                continue
            fields = line.split(maxsplit=27)
            if len(fields) < 27 or fields[3] not in reference:
                raise ValueError(f'{infernal}:{number}: incompatible cached Infernal table')
            start, end = sorted((int(fields[9]), int(fields[10])))
            if not 1 <= start <= end <= reference[fields[3]][0]:
                raise ValueError(f'{infernal}:{number}: coordinates outside cached assembly')
            count += 1
    entry['validation']['infernal_rows'] = count
    return genome, entry


def build_samplesheet(source, destination, records, project):
    with source.open(newline='') as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames
        rows = {re.sub(r'\s+', '-', r['name'].lower()): r for r in reader}
    prepared = []
    for genome, entry in sorted(records.items()):
        canonical = genome.replace('sp-', 'sp.-').removesuffix('-1')
        if canonical == 'strigamia-acuminata-3':
            canonical = 'strigamia-acuminata'
        row = dict(rows[canonical])
        name = row['name'].replace('sp.', 'sp')
        if genome.endswith(('-1', '-3')):
            name += ' ' + genome[-1]
        if re.sub(r'\s+', '-', name.lower()) != genome:
            raise ValueError(f'Cannot preserve genome ID: {genome}')
        row.update(name=name, fasta=entry['reference'], mnemonic=entry['mnemonic'])
        if row.get('RNASeq-dir'):
            reads = project / 'data/rnaseq' / canonical
            if not reads.is_dir():
                raise FileNotFoundError(reads)
            row['RNASeq-dir'] = str(reads)
        prepared.append(row)
    temporary = destination.with_name(destination.name + '.tmp')
    with temporary.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator='\n')
        writer.writeheader()
        writer.writerows(prepared)
    os.replace(temporary, destination)


def install_ncrna(cache, records):
    root = cache / 'annotation_export_v1/ncgenes_annotation'
    for genome, entry in records.items():
        if entry['trnascanse']:
            link(Path(entry['trnascanse']), root / 'trnascanse' / genome)
        for kingdom, path in entry['barrnap'].items():
            link(Path(path), root / f'barrnap_{kingdom}' / kingdom / genome / Path(path).name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache-dir', type=Path, required=True)
    parser.add_argument('--samplesheet', type=Path, required=True)
    parser.add_argument('--resume-samplesheet', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    cache = args.cache_dir.resolve()
    output = cache / 'cache_reuse_v1/imported'
    manifest_path = cache / 'cache_reuse_v1/manifest.json'
    if manifest_path.exists():
        raise FileExistsError(f'{manifest_path}: import already prepared')
    genomes = sorted(p.name for p in (cache / 'preprocessing/reformat_genome').iterdir() if p.is_dir())
    records = {}
    errors = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(prepare_genome, cache, output, genome) for genome in genomes]
        for future in as_completed(futures):
            try:
                genome, entry = future.result()
                records[genome] = entry
                print(f'{genome}: assembly verified; {len(entry["predictions"])} gene predictions; repeats and ncRNA validated', flush=True)
            except Exception as error:
                errors.append(str(error))
                print('VALIDATION ERROR: ' + str(error), flush=True)
    if errors:
        raise ValueError('Import not activated: ' + '; '.join(errors))
    library = cache / 'repeats_annotation/cd_hit/repeats_families/repeats_families.nr.fa'
    if not library.is_file():
        raise FileNotFoundError(library)
    build_samplesheet(args.samplesheet.resolve(), args.resume_samplesheet.resolve(), records, cache.parent)
    install_ncrna(cache, records)
    manifest = {'schema_version': 1, 'cache_dir': str(cache), 'samplesheet': str(args.resume_samplesheet.resolve()),
                'repeat_library': str(library), 'genomes': dict(sorted(records.items())),
                'notes': ['Original cache files are preserved.',
                          'Masked assemblies verified against cached reference sequence hashes, ignoring letter case and line wrapping.',
                          'Sequence references are translated in GFF files; valid coordinates and feature relationships are preserved.',
                          'Zero-length historical RepeatMasker entries (start=end+1) are excluded and retained in per-genome audit GFFs.',
                          'Empty or invalid gene predictions are not considered reusable and can be recalculated.',
                          'Historical repeat annotations and family classifications are reused as existing experimental results.',
                          'Derived gene filtering, functional annotation, QC, selection and exports use a separate cache namespace.',
                          'Historical software versions are not inferred.']}
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Manifest: {manifest_path}', flush=True)


if __name__ == '__main__':
    main()
