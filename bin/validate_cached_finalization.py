#!/usr/bin/env python3
"""Check that all requested historical annotations reached the final layout."""
import csv
import json
import sys
from pathlib import Path


def main():
    project = Path(__file__).resolve().parents[1]
    sources = json.loads((project / 'outputs/cache_reuse_20261023/final_coding_sources.json').read_text())
    out = project / 'results'
    errors = []
    required = [
        'genome/annotation.gff3', 'genome/sequences/genome.softmasked.fna',
        'genome/sequences/genome.fna', 'genome/stats/overview.tsv',
        'protein_coding/annotation.gff3', 'protein_coding/sequences/proteins.faa',
        'protein_coding/sequences/cds.fna', 'protein_coding/selection/cache_import.json',
        'noncoding/annotation.gff3', 'noncoding/sequences/ncrna.fna',
        'noncoding/stats/filtering.tsv', 'repeats/annotation.gff3',
        'repeats/plots/repeat_landscape.png', 'submission/metadata.json', 'checksums.sha256']
    for genome, source in sources.items():
        bundle = out / 'species' / genome
        for relative in required:
            path = bundle / relative
            if not path.is_file() or not path.stat().st_size:
                errors.append(f'{genome}: missing or empty {relative}')
        path = bundle / 'submission/metadata.json'
        if path.is_file():
            meta = json.loads(path.read_text())
            if meta['selected_method'] != source['method'] or meta['coding_transcripts'] == 0:
                errors.append(f'{genome}: missing/wrong selected coding annotation')
    species_dir = out / 'species'
    actual = {p.name for p in species_dir.iterdir() if p.is_dir()} if species_dir.exists() else set()
    if actual != set(sources):
        errors.append(f'Unexpected species set: missing={sorted(set(sources) - actual)}, extra={sorted(actual - set(sources))}')
    for relative in ['noncoding/stats/genome_metrics.tsv', 'protein_coding/stats/best_annotations.tsv', 'taxonomy.tsv']:
        path = out / 'comparative' / relative
        if not path.is_file():
            errors.append(f'Missing comparative/{relative}')
        else:
            with path.open() as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))
            if len(rows) != len(sources):
                errors.append(f'comparative/{relative}: {len(rows)} rows, expected {len(sources)}')
    plots = out / 'comparative/noncoding/plots'
    if not plots.is_dir() or not list(plots.glob('*.png')):
        errors.append('No comparative noncoding figures')
    result = {'status': 'FAILED' if errors else 'COMPLETE', 'genomes_expected': len(sources),
              'genomes_found': len(actual), 'errors': errors}
    destination = out / 'pipeline_info/finalization_validation.json'
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2) + '\n')
    state_path = project / 'outputs/myriapoda_run_20261023/run.json'
    state = json.loads(state_path.read_text())
    state.update(status=result['status'], validation=str(destination))
    state_path.write_text(json.dumps(state, indent=2) + '\n')
    print(json.dumps(result, indent=2))
    return bool(errors)


if __name__ == '__main__':
    sys.exit(main())
