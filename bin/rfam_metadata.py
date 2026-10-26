"""Read official Rfam types without a hand-maintained accession mapping."""
import re


def family_types(path):
    result = {}
    # Some official SQL dumps contain legacy encodings in free-text columns.
    # Accessions and the type column are ASCII; do not decode unrelated prose.
    with open(path, 'rb') as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip(b'\r\n').split(b'\t')
            if not fields[0].startswith(b'RF'):
                continue
            if len(fields) < 19:
                raise ValueError(f'{path}:{number}: expected Rfam family.txt with at least 19 columns')
            result[fields[0].decode('ascii')] = fields[18].decode('utf-8')
    return result


def classify_family(family_type):
    """Return GFF feature type, submission class and English display category."""
    tokens = [part.strip() for part in (family_type or '').split(';') if part.strip()]
    for token in tokens:
        if token.lower() in ('trna', 'rrna'):
            kind = 'tRNA' if token.lower() == 'trna' else 'rRNA'
            return kind, None, kind
    if any(token.lower() == 'cis-reg' for token in tokens):
        return 'ncRNA', 'other', 'Regulatory RNA'
    subtypes = [token for token in tokens if token.lower() not in ('gene', 'unknown', 'unclassified', '\\n')]
    if not subtypes:
        return 'ncRNA', 'other', 'Unclassified ncRNA'
    subtype = subtypes[0]
    if subtype.lower() == 'intron':
        return 'ncRNA', 'autocatalytically_spliced_intron', 'Intron'
    # Preserve future official subtypes in reports without inventing an NCBI
    # submission vocabulary value. More specific terms stay in rfam_type.
    allowed = {'antisense_RNA', 'guide_RNA', 'lncRNA', 'miRNA', 'piRNA', 'ribozyme',
               'scRNA', 'siRNA', 'snRNA', 'snoRNA', 'SRP_RNA', 'vault_RNA', 'Y_RNA'}
    return 'ncRNA', subtype if subtype in allowed else 'other', subtype


def infernal_accessions(tables):
    accessions = set()
    for table in tables:
        with open(table) as handle:
            for number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.split(maxsplit=27)
                if len(fields) < 27 or not re.fullmatch(r'RF\d+', fields[2]):
                    raise ValueError(f'{table}:{number}: expected cmscan --fmt 2 output')
                accessions.add(fields[2])
    return accessions
