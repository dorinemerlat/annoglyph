"""Exclude incoherent legacy models while retaining every original for review."""
import csv
from collections import defaultdict
from annoglyph_gff import read_gff


def filter_models(source, output, excluded_path, report_path, mapping, reference):
    # Invalid legacy coordinates must be readable here so their complete model
    # can be audited/excluded. Every downstream reader stays strict by default.
    features = read_gff(source, validate_coordinates=False)
    by_id, children = defaultdict(list), defaultdict(list)
    for f in features:
        for identifier in f.attrs.get('ID', []):
            by_id[identifier].append(f)
        for parent in f.attrs.get('Parent', []):
            children[parent].append(f)
    reasons = defaultdict(set)

    def reject(feature, reason):
        pending, seen, roots = [feature], set(), set()
        while pending:
            f = pending.pop()
            if id(f) in seen:
                continue
            seen.add(id(f))
            if f.kind in ('gene', 'mRNA', 'transcript'):
                roots.update(f.attrs.get('ID', []))
            else:
                pending.extend(p for identifier in f.attrs.get('Parent', []) for p in by_id[identifier])
        if not roots:
            roots.update(feature.attrs.get('ID', []))
        if not roots:
            raise ValueError(f'{source}: cannot isolate invalid feature without an ID')
        for identifier in roots:
            reasons[identifier].add(reason)

    for f in features:
        seqid = mapping.get(f.seqid.encode(), f.seqid.encode())
        if seqid not in reference:
            raise ValueError(f'{source}: unknown sequence {f.seqid}; check the assembly/ID mapping')
        if not 1 <= f.start <= f.end <= reference[seqid][0]:
            reject(f, 'coordinates_outside_assembly')
        if f.kind in ('gene', 'mRNA', 'transcript') and len(by_id[f.attrs['ID'][0]]) != 1:
            reject(f, 'duplicate_model_ID')
        for identifier in f.attrs.get('Parent', []):
            if not by_id[identifier]:
                reject(f, 'missing_parent')
            for parent in by_id[identifier]:
                if f.seqid != parent.seqid:
                    reject(parent, 'children_on_different_contigs')
                elif f.start < parent.start or f.end > parent.end:
                    reject(parent, 'children_outside_parent_bounds')
                if f.strand in ('+', '-') and parent.strand in ('+', '-') and f.strand != parent.strand:
                    reject(parent, 'children_on_different_strands')
    excluded_ids = set(reasons)
    excluded_rows = set()
    changed = True
    while changed:
        changed = False
        for f in features:
            own = set(f.attrs.get('ID', []))
            parents = f.attrs.get('Parent', [])
            orphaned = parents and all(p in excluded_ids for p in parents)
            empty_gene = f.kind == 'gene' and children[f.attrs['ID'][0]] and all(id(c) in excluded_rows for c in children[f.attrs['ID'][0]])
            if own & excluded_ids or orphaned or empty_gene:
                if id(f) not in excluded_rows:
                    excluded_rows.add(id(f))
                    excluded_ids.update(own)
                    changed = True
    kept, changed_parents = 0, 0
    with source.open() as incoming, output.open('w') as outgoing, excluded_path.open('w') as excluded:
        excluded.write('##gff-version 3\n')
        index = 0
        for line in incoming:
            if not line.strip() or line.startswith('#'):
                outgoing.write(line)
                continue
            feature = features[index]
            index += 1
            if id(feature) in excluded_rows:
                excluded.write(line)
                continue
            original_parents = feature.attrs.get('Parent', [])
            remaining = [p for p in original_parents if p not in excluded_ids]
            if remaining != original_parents:
                feature.attrs['Parent'] = remaining
                line = feature.line() + '\n'
                changed_parents += 1
            outgoing.write(line)
            kept += 1
    with report_path.open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['model_id', 'feature_type', 'seqid', 'reason'])
        for identifier, values in sorted(reasons.items()):
            feature = by_id[identifier][0]
            writer.writerow([identifier, feature.kind, feature.seqid, ';'.join(sorted(values))])
    return {'excluded_root_models': len(reasons), 'excluded_transcripts': sum(
                f.kind in ('mRNA', 'transcript') and id(f) in excluded_rows for f in features),
            'excluded_feature_rows': len(excluded_rows), 'retained_feature_rows': kept,
            'shared_feature_parent_references_updated': changed_parents,
            'reasons': {key: sorted(value) for key, value in reasons.items()}}
