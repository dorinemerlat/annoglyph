"""Resolve taxonomy-tree labels to annotation genome identifiers."""
import re
import pandas as pd
from ncrna_stats import phylogeny


def load_phylogeny(tree_file, genomes, taxonomy_file=None):
    genomes = sorted(genomes)
    if taxonomy_file:
        taxonomy = pd.read_csv(taxonomy_file, sep='\t', dtype={'id': str}).set_index('id', verify_integrity=True)
    else:
        # Compatibility for standalone plots using conventional species-based IDs.
        names = [re.sub(r'-\d+$', '', genome).replace('-', ' ') for genome in genomes]
        taxonomy = pd.DataFrame({'name': [name[:1].upper() + name[1:] for name in names]}, index=genomes)
    return phylogeny(tree_file, taxonomy, genomes)
