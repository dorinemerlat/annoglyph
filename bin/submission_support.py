"""Explicit source/genetic-code metadata for annotation exports."""
import csv
import json
from pathlib import Path
from Bio.Data import CodonTable


def load_metadata(path, sample, assets_path=None):
    if not path:
        return {}, {}
    document = json.loads(Path(path).read_text())
    config = document.get('genomes', document).get(sample, {})
    if not config:
        raise ValueError(f'No submission metadata for genome {sample}')
    aliases = {}
    if assets_path:
        for asset in json.loads(Path(assets_path).read_text()):
            source = Path(asset['source'])
            if not source.is_file() or not source.name.endswith('.tsv'):
                continue
            with source.open() as handle:
                reader = csv.DictReader(handle, delimiter='\t')
                if reader.fieldnames == ['previous_IDs', 'AnnoGlyph_IDs']:
                    aliases.update({row['previous_IDs']: row['AnnoGlyph_IDs'] for row in reader})
    sequences = {}
    for identifier, values in config.get('sequences', {}).items():
        target = aliases.get(identifier, identifier)
        if target in sequences:
            raise ValueError(f'Duplicate sequence metadata for {target}')
        if values.get('organelle') and 'translation_table' not in values:
            raise ValueError(f'Provide an explicit genetic code for organellar sequence {identifier}')
        if 'translation_table' in values and int(values['translation_table']) not in CodonTable.unambiguous_dna_by_id:
            raise ValueError(f'Unsupported genetic code for {identifier}')
        if values.get('topology', 'linear') not in ('linear', 'circular'):
            raise ValueError(f'Invalid topology for {identifier}')
        sequences[target] = values
    return config, sequences


def modifiers(config, sequence, organism, default_code):
    values = {'organism': organism, 'moltype': 'genomic', **config.get('source', {})}
    for key in ['organelle', 'topology', 'chromosome', 'plasmid']:
        if key in sequence:
            values[{'organelle': 'location', 'plasmid': 'plasmid-name'}.get(key, key)] = sequence[key]
    code = sequence.get('translation_table', default_code)
    code_key = 'mgcode' if sequence.get('organelle') == 'mitochondrion' else ('pgcode' if sequence.get('organelle') in ('plastid', 'chloroplast') else 'gcode')
    values[code_key] = code
    for key, value in values.items():
        if any(character in str(value) for character in '[]\n\r'):
            raise ValueError(f'Invalid FASTA source modifier: {key}')
    return ' '.join(f'[{key}={value}]' for key, value in values.items())
