"""Small GFF3 helpers shared by ncRNA reconciliation and final export."""
from dataclasses import dataclass
from urllib.parse import quote, unquote


@dataclass
class Feature:
    seqid: str
    source: str
    kind: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attrs: dict

    def line(self):
        attrs = ';'.join(quote(k, safe=':_-') + '=' + ','.join(
            quote(v, safe=':_-.|') for v in values) for k, values in self.attrs.items()) or '.'
        return '\t'.join(map(str, [self.seqid, self.source, self.kind, self.start,
                                   self.end, self.score, self.strand, self.phase, attrs]))


def read_gff(path):
    features = []
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith('##FASTA'):
                break
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'{path}:{number}: expected nine GFF3 columns')
            attrs = {}
            if fields[8] != '.':
                for item in fields[8].strip(';').split(';'):
                    key, sep, value = item.partition('=')
                    if not sep:
                        raise ValueError(f'{path}:{number}: malformed GFF3 attribute: {item}')
                    attrs[unquote(key)] = [unquote(v) for v in value.split(',')]
            start, end = int(fields[3]), int(fields[4])
            if start < 1 or end < start:
                raise ValueError(f'{path}:{number}: invalid coordinates')
            features.append(Feature(*fields[:3], start, end, *fields[5:8], attrs))
    return features


def namespace(features, prefix):
    """Namespace IDs and references, preserving legal multipart feature IDs."""
    mapping = {}
    for f in features:
        for old in f.attrs.get('ID', []):
            mapping.setdefault(old, f'{prefix}.{len(mapping) + 1}')
    for index, f in enumerate(features, 1):
        for key in ['ID', 'Parent', 'Derives_from']:
            if key in f.attrs:
                try:
                    f.attrs[key] = [mapping[v] for v in f.attrs[key]]
                except KeyError as error:
                    raise ValueError(f'{prefix}: unresolved {key}: {error}') from error
        f.attrs.setdefault('ID', [f'{prefix}.row{index}'])
    return mapping


def sort_key(feature):
    return (feature.seqid, feature.start, -feature.end,
            0 if feature.kind == 'gene' else 1, feature.attrs['ID'][0])
